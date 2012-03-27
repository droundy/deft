\begin{code}
{-# LANGUAGE GADTs, PatternGuards #-}
module CodeGen ( RealSpace, r_var,
                 KSpace(..), k_var, kx, ky, kz, k, ksqr,
                 Scalar(..), s_var,
                 fft, ifft, integrate, grad, derive,
                 Expression, joinFFTs, (===), var,
                 Statement(..),
                 Type, 
                 makeHomogeneous, isConstant, hasexpression, factorandsum, -- for debugging only!!!
                 code, latex, setZero, codeStatement, substitute, cleanvars,
                 generateHeader, simp2, countFFT, checkDup, peakMem, reuseVar, removeAt, findToDo, latexSimp, hasFFT )
       where

import Debug.Trace

import qualified Data.Map as Map
import Data.List ( nubBy, (\\), nub, subsequences )

\end{code}

The \verb!RealSpaceField! data type describes a field in real space.

\begin{code}
data RealSpace = IFFT (Expression KSpace)
               deriving ( Eq, Ord, Show )
data KSpace = Delta | -- handy for FFT of homogeneous systems
              Kx | Ky | Kz |
              FFT (Expression RealSpace)
            deriving ( Eq, Ord, Show )
data Scalar = Integrate (Expression RealSpace)
            deriving ( Eq, Ord, Show )

kinversion :: Expression KSpace -> Expression KSpace
kinversion (Var _ _ _ (Just e)) = kinversion e
kinversion e@(Var _ _ _ Nothing) = e
kinversion (Scalar e) = Scalar e
kinversion (Cos e) = cos (kinversion e)
kinversion (Sin e) = sin (kinversion e)
kinversion (Exp e) = exp (kinversion e)
kinversion (Log e) = log (kinversion e)
kinversion (Abs e) = abs (kinversion e)
kinversion (Signum e) = signum (kinversion e)
kinversion (Product p) = product $ map ff $ product2pairs p
  where ff (e,n) = (kinversion e) ** toExpression n
kinversion (Sum s) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, kinversion y)
kinversion (Expression Kx) = -kx
kinversion (Expression Ky) = -ky
kinversion (Expression Kz) = -kz
kinversion (Expression x) = Expression x

instance Code RealSpace where
  codePrec _ (IFFT (Var _ ksp _ Nothing)) = showString ("ifft(gd, " ++ksp++ ")")
  codePrec _ (IFFT ke) = showString "ifft(gd, " . codePrec 0 ke . showString ")"
  latexPrec _ (IFFT ke) = showString "\\text{ifft}\\left(" . latexPrec 0 ke . showString "\\right)"
instance Type RealSpace where
  isRealSpace _ = Same
  derivativeHelper v ddr (IFFT ke) = derive v (fft ddr) (kinversion ke)
  zeroHelper v (IFFT ke) = ifft (setZero v ke)
  prefix e = case isifft e of Just _ -> ""
                              Nothing -> "for (int i=0; i<gd.NxNyNz; i++) {\n\t\t"
  postfix e = case isifft e of Just _ -> ""
                               Nothing -> "\t\n}\n"
  codeStatementHelper a op (Expression (IFFT (Var _ v _ Nothing))) = a ++ op ++ "ifft(gd, " ++ v ++ ");\n"
  codeStatementHelper _ _ (Expression (IFFT e)) = error ("It is a bug to generate code for a non-var input to ifft\n"++ latex e)
  codeStatementHelper a op (Var _ _ _ (Just e)) = codeStatementHelper a op e
  codeStatementHelper a op e = prefix e ++ a ++ "[i]" ++ op ++ code e ++ ";\n" ++ postfix e
  initialize (Var _ x _ Nothing) = "VectorXd " ++ x ++ "(gd.NxNyNz);"
  initialize _ = "VectorXd output(gd.NxNyNz);"
  free (Var _ x _ Nothing) = x ++ ".resize(0); // Realspace"
  free _ = error $ trace "free error" "free error"
  toScalar (IFFT ke) = makeHomogeneous ke

instance Code KSpace where
  codePrec _ Kx = showString "k_i[0]"
  codePrec _ Ky = showString "k_i[1]"
  codePrec _ Kz = showString "k_i[2]"
  codePrec _ Delta = showString "delta(k?)"
  codePrec _ (FFT r) = showString "fft(gd, " . codePrec 0 (makeHomogeneous r) . showString ")"
  latexPrec _ Kx = showString "k_{x}"
  latexPrec _ Ky = showString "k_{y}"
  latexPrec _ Kz = showString "k_{z}"
  latexPrec _ Delta = showString "\\delta(k)"
  latexPrec _ (FFT r) = showString "\\text{fft}\\left(" . latexPrec 0 r . showString "\\right)"
instance Type KSpace where
  isKSpace _ = Same
  derivativeHelper v ddk (FFT r) = derive v (ifft ddk) r
  derivativeHelper _ _ _ = 0
  zeroHelper v (FFT r) = fft (setZero v r)
  zeroHelper _ x = Expression x
  prefix e = case isfft e of Just _ -> ""
                             Nothing -> "for (int i=1; i<gd.NxNyNzOver2; i++) {\n\t\tconst int z = i % gd.NzOver2;\n\t\tconst int n = (i-z)/gd.NzOver2;\n\t\tconst int y = n % gd.Ny;\n\t\tconst int xa = (n-y)/gd.Ny;\n\t\tconst RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);\n\t\tconst Reciprocal k_i = gd.Lat.toReciprocal(rvec);\n\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3); assert(dr);\n"
  postfix e = case isfft e of Just _ -> ""
                              Nothing -> "\t}\n"
  codeStatementHelper a op (Expression (FFT (Var _ v _ Nothing))) = a ++ op ++ "fft(gd, " ++ v ++ ");\n"
  codeStatementHelper _ _ (Expression (FFT _)) = error "It is a bug to generate code for a non-var input to fft"
  codeStatementHelper a op (Var _ _ _ (Just e)) = codeStatementHelper a op e
  codeStatementHelper a op e =
          if k0code == "0"
          then a ++ "[0]" ++ op ++ "0;\n\t" ++ prefix e ++ "\t\t" ++ a ++ "[i]"  ++ op ++ code e ++ ";" ++ postfix e
          else "{\n\t\tconst int i = 0;\n\t\tconst Reciprocal k_i = Reciprocal(0,0,0);\n\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3); assert(dr);\n\t\t" ++ a ++ "[0]" ++ op ++ code (setZero kz (setZero ky (setZero kx e))) ++ ";\n\t}\n\t" ++ prefix e ++ "\t\t" ++ a ++ "[i]"  ++ op ++ code e ++ ";" ++ postfix e
      where k0code = code (setZero kz (setZero ky (setZero kx e)))
  initialize (Var _ x _ Nothing) = "VectorXcd " ++ x ++ "(gd.NxNyNzOver2);"
  initialize _ = "VectorXcd output(gd.NxNyNzOver2);"
  free (Var _ x _ Nothing) = x ++ ".resize(0); // KSpace"
  free _ = error "free error"
  toScalar Delta = 1
  toScalar Kx = s_var "_kx"
  toScalar Ky = 0
  toScalar Kz = 0
  toScalar (FFT e) = makeHomogeneous e

mapExpression :: (Type a, Type b) => (a -> Expression b) -> Expression a -> Expression b
mapExpression f (Var _ _ _ (Just e)) = mapExpression f e
mapExpression _ (Var c v t Nothing) = Var c v t Nothing
mapExpression _ (Scalar e) = Scalar e
mapExpression f (Cos e) = cos (mapExpression f e)
mapExpression f (Sin e) = sin (mapExpression f e)
mapExpression f (Exp e) = exp (mapExpression f e)
mapExpression f (Log e) = log (mapExpression f e)
mapExpression f (Abs e) = abs (mapExpression f e)
mapExpression f (Signum e) = signum (mapExpression f e)
mapExpression f (Product p) = product $ map ff $ product2pairs p
  where ff (e,n) = (mapExpression f e) ** toExpression n
mapExpression f (Sum s) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, mapExpression f y)
mapExpression f (Expression x) = f x

cleanvars :: Type a => Expression a -> Expression a
cleanvars (Var _ _ _ (Just e)) = cleanvars e
cleanvars (Var c v t Nothing) = Var c v t Nothing
cleanvars (Scalar e) = Scalar (cleanvars e)
cleanvars (Cos e) = cos (cleanvars e)
cleanvars (Sin e) = sin (cleanvars e)
cleanvars (Exp e) = exp (cleanvars e)
cleanvars (Log e) = log (cleanvars e)
cleanvars (Abs e) = abs (cleanvars e)
cleanvars (Signum e) = signum (cleanvars e)
cleanvars (Product p) = product $ map ff $ product2pairs p
  where ff (e,n) = (cleanvars e) ** toExpression n
cleanvars (Sum s) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, cleanvars y)
cleanvars (Expression e)
  | Same <- isKSpace (Expression e), FFT e' <- e = fft (cleanvars e')
  | Same <- isRealSpace (Expression e), IFFT e' <- e = ifft (cleanvars e')
  | Same <- isScalar (Expression e), Integrate e' <- e = integrate (cleanvars e')
  | otherwise = Expression e

isEven :: (Type a, Type b) => Expression b -> Expression a -> Double
isEven v e | Same <- compareExpressions v e = -1
isEven v (Var _ _ _ (Just e)) = isEven v e
isEven _ (Var _ _ _ Nothing) = 1
isEven v (Scalar e) = isEven v e
isEven _ (Cos _) = 1
isEven v (Sin e) = isEven v e
isEven v (Exp e) = if isEven v e == 1 then 1 else 0
isEven v (Log e) = if isEven v e == 1 then 1 else 0
isEven _ (Abs _) = 1
isEven v (Signum e) = isEven v e
isEven v (Product p) = product $ map ie $ product2pairs p
  where ie (x,n) = isEven v x ** n
isEven _ (Sum s) | Sum s == 0 = 1 -- ???
isEven v (Sum s) = ie (isEven v x) xs
  where (_,x):xs = sum2pairs s
        ie sofar ((_,y):ys) = if isEven v y /= sofar 
                              then 0
                              else ie sofar ys
        ie sofar [] = sofar
isEven v (Expression e) 
  | Same <- isKSpace (Expression e) =
    case e of
      FFT r -> isEven v r
      _ -> 1
isEven v (Expression e) 
  | Same <- isRealSpace (Expression e) =
    case e of
      IFFT ks -> isEven v ks
isEven v (Expression e) 
  | Same <- isScalar (Expression e) =
    case e of
      Integrate x -> isEven v x
isEven _ (Expression _) = 1 -- Technically, it might be good to recurse into this

setZero :: (Type a, Type b) => Expression b -> Expression a -> Expression a
setZero v e | Same <- compareExpressions v e = 0
setZero v (Var _ _ _ (Just e)) = setZero v e
setZero _ e@(Var _ _ _ Nothing) = e
setZero v (Scalar e) = Scalar (setZero v e)
setZero v (Cos e) = cos (setZero v e)
setZero v (Sin e) = sin (setZero v e)
setZero v (Exp e) = exp (setZero v e)
setZero v (Log e) = log (setZero v e)
setZero v (Abs e) = abs (setZero v e)
setZero v (Signum e) = signum (setZero v e)
setZero v (Product p) | product2denominator p == 1 = product $ map ff $ product2pairs p
  where ff (e,n) = (setZero v e) ** toExpression n
setZero v (Product p) = 
  if isEven v (Product p) == -1
  then 0
  else 
    if zd /= 0
    then zn / zd
    else if zn /= 0
         then error ("L'Hopital's rule failure: " ++ latex n ++ "\n /\n  " ++ latex d ++ "\n\n\n" 
                     ++ latex (Product p) ++ "\n\n\n" ++ latex zn)
         else case compareTypes v n of
                Same -> setZero v (dtop / dbot)
                  where dtop = derive v 1 n
                        dbot = derive v 1 d
                Different -> error "oops, failure in setZero that makes no sense."
  where d = product2denominator p
        n = product $ product2numerator p
        zn = setZero v n
        zd = setZero v d
setZero _ (Sum s) | Sum s == 0 = 0
setZero v (Sum s) = out
  where sz (f,x) = (f, setZero v x)
        out = pairs2sum $ map sz $ sum2pairs s
setZero v (Expression x) = zeroHelper v x

instance Code Scalar where
  codePrec _ (Integrate r) = showString "integrate(" . codePrec 0 r . showString ")"
  latexPrec _ (Integrate r) = showString "integrate(" . latexPrec 0 r . showString ")"
instance Type Scalar where
  s_var ("complex(0,1)") = Var "complex(0,1)" "complex(0,1)" "i" Nothing
  s_var v@['d',_] = Var v v v Nothing -- for differentials
  s_var vv@(a:v@(_:_)) = Var vv vv (a : '_' : '{' : v ++ "}") Nothing
  s_var v = Var v v v Nothing
  isScalar _ = Same
  derivativeHelper v dds (Integrate e) = derive v (Scalar dds*s_var "dV") e
  zeroHelper v (Integrate e) = integrate (setZero v e)
  codeStatementHelper a _ (Expression (Integrate e)) = "for (int i=0; i<gd.NxNyNz; i++) {\n\t\t" ++ a ++ " += " ++ code (e * s_var "gd.dvolume") ++ ";\n\t}\n"
  codeStatementHelper a op e = a ++ op ++ code e ++ ";\n\t"
  prefix (Expression (Integrate _))  = "for (int i=0; i<gd.NxNyNz; i++) {    "
  prefix _ = ""
  postfix (Expression (Integrate _)) = "\n}\n"
  postfix _ = ""
  initialize (Expression (Integrate _)) = "double output = 0;\n"
  initialize (Var _ x _ Nothing) = "double " ++ x ++ " = 0;\n"
  initialize _ = "double output = 0;\n"
  toScalar (Integrate r) = makeHomogeneous r
  fromScalar = id

r_var :: String -> Expression RealSpace
r_var v@([_]) = Var (v++"[i]") v v Nothing
r_var (a:v) = Var (a:v++"[i]") (a:v) (a:'_':'{':v++"}") Nothing
r_var "" = error "r_var needs non-empty string"

k_var :: String -> Expression KSpace
k_var v@([_]) = Var (v++"[i]") v v Nothing
k_var (a:v) = Var (a:v++"[i]") (a:v) (a:'_':'{':v++"}") Nothing
k_var "" = error "r_var needs non-empty string"

infix 4 ===

(===) :: Type a => String -> Expression a -> Expression a
--_ === e = e
v@(a:r@(_:_)) === e = Var c v ltx (Just e)
  where ltx = a : "_{"++r++"}"
        c = (case isScalar e of
              Same -> v
              Different -> v ++ "[i]") :: String
v === e = Var c v v (Just e)
  where c = (case isScalar e of
              Same -> v
              Different -> v ++ "[i]") :: String

var :: Type a => String -> String -> Expression a -> Expression a
var v ltx e = Var c v ltx (Just e)
  where c = (case isScalar e of
              Same -> v
              Different -> v ++ "[i]") :: String

kx :: Expression KSpace
kx = Expression Kx

ky :: Expression KSpace
ky = Expression Ky

kz :: Expression KSpace
kz = Expression Kz

k :: Expression KSpace
k = sqrt ksqr

ksqr :: Expression KSpace
ksqr = kx**2 + ky**2 + kz**2

integrate :: Type a => Expression RealSpace -> Expression a
integrate x = fromScalar $ Expression (Integrate x)

fft :: Expression RealSpace -> Expression KSpace
fft (Scalar e) = Scalar e * Expression Delta
fft r | r == 0 = 0
      | otherwise = Expression (FFT r)

ifft :: Expression KSpace -> Expression RealSpace
ifft ke | ke == 0 = 0
        | otherwise = Expression (IFFT ke)
\end{code}

The \verb!ReciprocalSpaceField! data type describes a field in real space.

\begin{code}
\end{code}

The \verb!Expression! data type describes an expression with a given
type.

\begin{code}
data Expression a = Scalar (Expression Scalar) |
                    Var String String String (Maybe (Expression a)) | -- A variable with a possible value
                    Expression a |
                    Cos (Expression a) |
                    Sin (Expression a) |
                    Exp (Expression a) |
                    Log (Expression a) |
                    Abs (Expression a) |
                    Signum (Expression a) |
                    Product (Map.Map (Expression a) Double) |
                    Sum (Map.Map (Expression a) Double)
              deriving (Eq, Ord, Show)

sum2pairs :: Type a => Map.Map (Expression a) Double -> [(Double, Expression a)]
sum2pairs s = map rev $ Map.assocs s
  where rev (a,b) = (b,a)

pairs2sum :: Type a => [(Double, Expression a)] -> Expression a
pairs2sum s = helper $ filter ((/= 0) . snd) $ filter ((/= 0) . fst) s
  where helper [] = 0
        helper [(1,e)] = e
        helper es = Sum $ fl (Map.empty) es
        fl a [] = a
        fl a ((f,Sum s'):xs) = fl a (map mulf (sum2pairs s') ++ xs)
          where mulf (ff,e) = (ff*f, e)
        fl a ((f,x):xs) = case Map.lookup x a of
                          Just f' -> if f + f' == 0
                                     then fl (Map.delete x a) xs
                                     else fl (Map.insert x (f+f') a) xs
                          Nothing -> fl (Map.insert x f a) xs
        

product2pairs :: Type a => Map.Map (Expression a) Double -> [(Expression a, Double)]
product2pairs s = Map.assocs s

map2product :: Type a => Map.Map (Expression a) Double -> Expression a
map2product p | Map.size p == 1 = case product2pairs p of [(e,1)] -> e
                                                          _ -> Product p
map2product p = helper 1 (Map.empty) $ product2pairs p
  where helper 1 a [] = Product a
        helper f a [] = Sum $ Map.singleton (Product a) f
        helper f a ((Sum x,n):xs) | [(f',x')] <- sum2pairs x = helper (f*f'**n) a ((x',n):xs)
        helper f a ((x,n):xs) = helper f (Map.insert x n a) xs

pairs2product :: Type a => [(Expression a, Double)] -> Expression a
pairs2product = map2product . fl (Map.empty)
  where fl a [] = a
        fl a ((_,0):xs) = fl a xs
        fl a ((Product p, n):xs) = fl a (map tonpower (product2pairs p) ++ xs)
            where tonpower (e,n') = (e, n'*n)
        fl a ((x,n):xs) = case Map.lookup x a of
                            Just n' -> if n + n' == 0
                                       then fl (Map.delete x a) xs
                                       else fl (Map.insert x (n+n') a) xs
                            Nothing -> fl (Map.insert x n a) xs

product2numerator :: Type a => Map.Map (Expression a) Double -> [Expression a]
product2numerator s = map f $ product2numerator_pairs s
  where f (a,b) = a ** (toExpression b)

product2numerator_pairs :: Type a => Map.Map (Expression a) Double -> [(Expression a, Double)]
product2numerator_pairs s = filter ((>=0) . snd) $ product2pairs s

product2denominator :: Type a => Map.Map (Expression a) Double -> Expression a
product2denominator s = pairs2product $ map n $ filter ((<0) . snd) $ product2pairs s
  where n (a,b) = (a, -b)

instance (Type a, Code a) => Code (Expression a) where
  codePrec = codeE
  latexPrec = latexE

codeE :: (Type a, Code a) => Int -> Expression a -> ShowS
codeE _ (Var c _ _ Nothing) = showString c
codeE p (Var _ _ _ (Just e)) = codeE p e
codeE p (Scalar x) = codePrec p x
codeE p (Expression x) = codePrec p x
codeE _ (Cos x) = showString "cos(" . codeE 0 x . showString ")"
codeE _ (Sin x) = showString "sin(" . codeE 0 x . showString ")"
codeE _ (Exp x) = showString "exp(" . codeE 0 x . showString ")"
codeE _ (Log x) = showString "log(" . codeE 0 x . showString ")"
codeE _ (Abs x) = showString "fabs(" . codeE 0 x . showString ")"
codeE _ (Signum _) = undefined
codeE _ (Product p) | Product p == 1 = showString "1"
codeE pree (Product p) = showParen (pree > 7) $
                         if den == 1 
                         then codesimple num
                         else codesimple num . showString "/" . codeE 8 den
  where codesimple [] = showString "1"
        codesimple [(a,n)] = codee a n
        codesimple [(a,n),(b,m)] = codee a n . showString "*" . codee b m
        codesimple ((a,n):es) = codee a n . showString "*" . codesimple es
        num = product2numerator_pairs p
        den = product2denominator p
        codee _ 0 = showString "1" -- this shouldn't happen...
        codee _ n | n < 0 = error "shouldn't have negative power here"
        codee x 1 = codeE 7 x
        codee x 0.5 = showString "sqrt(" . codePrec 0 x . showString ")"
        codee x nn
          | fromInteger n2 == 2*nn && odd n2 = codee x 0.5 . showString "*" . codee x (nn-0.5)
          | fromInteger n == nn && odd n = codee x 1 . showString "*" . codee x (nn-1)
          | fromInteger n == nn = 
            showParen (nn/2>1) (codee x (nn / 2)) . showString "*" . showParen (nn/2>1) (codee x (nn / 2))
          where n2 = floor (2*nn)
                n = floor nn
        codee x n = showString "pow(" . codeE 0 x . showString (", " ++ show n ++ ")")
codeE _ (Sum s) | Sum s == 0 = showString "0"
codeE p (Sum s) = showParen (p > 6) (showString me)
  where me = foldl addup "" $ sum2pairs s
        addup "" (1,e) = codeE 6 e ""
        addup "" (f,e) = if e == 1
                         then show f
                         else show f ++ "*" ++ codeE 7 e ""
        addup rest (1,e) = codeE 6 e (showString " + " $ rest)
        addup rest (f,e) = show f ++ "*" ++ codeE 7 e (showString " + " $ rest)

latexE :: (Type a, Code a) => Int -> Expression a -> ShowS
latexE p (Var _ "" "" (Just e)) = latexE p e
latexE _ (Var _ c "" _) = showString c
latexE _ (Var _ _ t _) = showString t
latexE _ x | Just xx <- isConstant x = showString (latexDouble xx)
latexE p (Scalar x) = latexPrec p x
latexE p (Expression x) = latexPrec p x
latexE _ (Cos x) = showString "\\cos(" . latexE 0 x . showString ")"
latexE _ (Sin x) = showString "\\sin(" . latexE 0 x . showString ")"
latexE _ (Exp x) = showString "\\exp(" . latexE 0 x . showString ")"
latexE _ (Log x) = showString "\\log(" . latexE 0 x . showString ")"
latexE _ (Abs x) = showString "\\left|" . latexE 0 x . showString "\\right|"
latexE _ (Signum _) = undefined
latexE p (Product x) | Map.size x == 1 && product2denominator x == 1 =
  case product2pairs x of
    [(_,0)] -> showString "1" -- this shouldn't happen...
    [(_, n)] | n < 0 -> error "shouldn't have negative power here"
    [(e, 1)] -> latexE p e
    [(e, 0.5)] -> showString "\\sqrt{" . latexE 0 e . showString "}"
    [(e, nn)]
      | fromInteger n2 == 2*nn && odd n2 -> if n2 < 10
                                            then latexE 8 e . showString ("^{\\frac" ++ show n2 ++ "2}")
                                            else latexE 8 e . showString ("^{\\frac{" ++ show n2 ++ "}2}")
      | fromInteger n == nn -> if n2 < 10
                               then latexE 8 e . showString ("^" ++ show n)
                               else latexE 8 e . showString ("^{" ++ show n ++ "}")
          where n2 = floor (2*nn)
                n = floor nn
    [(e,n)] -> latexE 8 e . showString ("^{" ++ show n ++ "}")
    _ -> error "This really cannot happen."
latexE pree (Product p) | product2denominator p == 1 = latexParen (pree > 7) $ ltexsimple $ product2numerator p
  where ltexsimple [] = showString "1"
        ltexsimple [a] = latexE 7 a
        ltexsimple [a,b] = latexE 7 a . showString " " . latexE 7 b
        ltexsimple (a:es) = latexE 7 a . showString " " . ltexsimple es
latexE pree (Product p) = latexParen (pree > 7) $ showString "\\frac{" . latexE 0 num . showString "}{" .
                                                                        latexE 0 den . showString "}"
  where num = product $ product2numerator p
        den = product2denominator p
latexE p (Sum s) = latexParen (p > 6) (showString me)
  where me = foldl addup "" $ sum2pairs s
        addup "" (1,e) = latexE 6 e ""
        addup "" (f,e) | f < 0 = "-" ++ addup "" (-f, e)
        addup "" (f,e) = if e == 1
                         then latexDouble f
                         else latexDouble f ++ " " ++ latexE 6 e ""
        addup ('-':rest) (1,e) = latexE 6 e (" - " ++ rest)
        addup ('-':rest) (f,e) = latexDouble f ++ " " ++ latexE 7 e (showString " - " $ rest)
        addup rest (1,e) = latexE 6 e (" + " ++ rest)
        addup rest (f,e) = latexDouble f ++ " " ++ latexE 7 e (showString " + " $ rest)

latexParen :: Bool -> ShowS -> ShowS
latexParen False x = x
latexParen True x = showString "\\left(" . x . showString "\\right)"

latexDouble :: Double -> String
latexDouble 0 = "0"
latexDouble x | x < 0 = "-" ++ latexDouble (-x)
latexDouble x = case double2int x of
                  Just n -> show n
                  Nothing ->
                    case double2int (x/pi) of
                      Just n -> show n ++ "\\pi"
                      Nothing ->
                        case double2int (pi/x) of
                          Just n -> "\\frac{\\pi}{" ++ show n ++ "}"
                          Nothing -> printpi [2 .. 36]
  where printpi [] = show x
        printpi (n:ns) = case double2int (fromInteger n*x/pi) of
                           Just n' -> "\\frac{" ++ show n' ++ "}{" ++ show n ++ "}\\pi"
                           Nothing ->
                             case double2int (x*fromInteger n*pi) of
                               Just n' -> "\\frac{" ++ show n' ++ "}{" ++ show n ++ "\\pi}"
                               Nothing -> printpi ns

double2int :: Double -> Maybe Int
double2int f = if abs(fromInteger (round f) - f) < 1e-13*f
               then Just (round f :: Int)
               else Nothing

class Code a  where
    codePrec  :: Int -> a -> ShowS
    codePrec _ x s = code x ++ s
    code      :: a -> String 
    code x = codePrec 0 x ""
    latexPrec :: Int -> a -> ShowS
    latexPrec _ x s = latex x ++ s
    latex     :: a -> String
    latex x = latexPrec 0 x ""

data Same a b where
    Same :: Same a a
    Different :: Same a b

toExpression :: (Type a, Real x) => x -> Expression a
toExpression 0 = Sum $ Map.empty
toExpression 1 = Product $ Map.empty
toExpression x = Sum $ Map.singleton 1 (fromRational $ toRational x)

isConstant :: Type a => Expression a -> Maybe Double
isConstant (Sum s) = case sum2pairs s of
                       [] -> Just 0
                       [(x,1)] -> Just x
                       _ -> Nothing
isConstant (Product p) = if Map.size p == 0 then Just 1 else Nothing
isConstant _ = Nothing

class (Ord a, Show a, Code a) => Type a where 
  isScalar :: Expression a -> Same a Scalar
  isScalar _ = Different
  isRealSpace :: Expression a -> Same a RealSpace
  isRealSpace _ = Different
  isKSpace :: Expression a -> Same a KSpace
  isKSpace _ = Different
  s_var :: String -> Expression a
  s_var = Scalar . s_var
  derivativeHelper :: Type b => Expression b -> Expression a -> a -> Expression b
  zeroHelper :: Type b => Expression b -> a -> Expression a
  codeStatementHelper :: String -> String -> Expression a -> String
  prefix :: Expression a -> String
  postfix :: Expression a -> String
  postfix _ = ""
  initialize :: Expression a -> String
  initialize _ = ""
  free :: Expression a -> String
  free _ = error "free nothing"
  toScalar :: a -> Expression Scalar
  fromScalar :: Expression Scalar -> Expression a
  fromScalar = Scalar



codeStatements :: [Statement] -> String
codeStatements x = unlines $ map codeStatement (reuseVar x)

codeStatement :: Statement -> String
codeStatement (AssignR x y) = "\t" ++ codeStatementHelper x " = " y
codeStatement (AssignK x y) = "\t" ++ codeStatementHelper x " = " y
codeStatement (AssignS x y) = "\t" ++ codeStatementHelper x " = " y
codeStatement (InitializeR e) = "\t" ++ initialize e
codeStatement (InitializeK e) = "\t" ++ initialize e
codeStatement (InitializeS e) = "\t" ++ initialize e
codeStatement (FreeR e) = "\t" ++ free e
codeStatement (FreeK e) = "\t" ++ free e

countFFT :: [Statement] -> Int
countFFT = sum . map helper
    where helper (AssignR _ (Expression (IFFT _))) = 1
          helper (AssignK _ (Expression (FFT _))) = 1
          helper (AssignR _ (Var _ _ _ (Just (Expression (IFFT _))))) = 1
          helper (AssignK _ (Var _ _ _ (Just (Expression (FFT _))))) = 1
          helper _ = 0

peakMem :: [Statement] -> Int
peakMem = maximum . (helper 0) . freeVectors
    where helper n (x:xs) | (InitializeR _) <- x = (n+1) : (helper (n+1) xs)
                          | (InitializeK _) <- x = (n+1) : (helper (n+1) xs)
                          | (FreeR _) <- x = (n-1) : (helper (n-1) xs)
                          | (FreeK _) <- x = (n-1) : (helper (n-1) xs)
                          | otherwise = n : helper n xs
          helper n [] = [n]

makeHomogeneous :: Type a => Expression a -> Expression Scalar
makeHomogeneous ee = 
  scalarScalar $ case isKSpace ee of
                    Same -> setZero (s_var "_kx" :: Expression Scalar) $ mapExpression toScalar ee
                    _ -> mapExpression toScalar ee
  where scalarScalar :: Expression Scalar -> Expression Scalar
        scalarScalar (Var _ _ _ (Just e)) = scalarScalar e
        scalarScalar (Var _ c _ Nothing) = s_var c
        scalarScalar (Scalar s) = s
        scalarScalar (Sum x) = pairs2sum $ map f $ sum2pairs x
          where f (a,b) = (a, scalarScalar b)
        scalarScalar (Product x) = pairs2product $ map f $ product2pairs x
          where f (a,b) = (scalarScalar a, b)
        scalarScalar (Expression e) = Expression e -- FIXME
        scalarScalar (Cos x) = cos (scalarScalar x)
        scalarScalar (Sin x) = sin (scalarScalar x)
        scalarScalar (Exp x) = exp (scalarScalar x)
        scalarScalar (Log x) = log (scalarScalar x)
        scalarScalar (Abs x) = abs (scalarScalar x)
        scalarScalar (Signum x) = signum (scalarScalar x)

instance Type a => Num (Expression a) where
  x + y | Just 0 == isConstant x = y
        | Just 0 == isConstant y = x
  Sum a + Sum b = sumup a $ sum2pairs b
      where sumup x [] = Sum x
            sumup x ((f,y):ys) = case Map.lookup y x of
                                 Nothing -> sumup (Map.insert y f x) ys
                                 Just f' -> if f + f' == 0
                                            then sumup (Map.delete y x) ys
                                            else sumup (Map.insert y (f+f') x) ys
  Sum a + b = case Map.lookup b a of
                Just fac -> if fac + 1 == 0
                            then case sum2pairs deleted of 
                                   [(1,e)] -> e
                                   [(f,e)] -> Sum $ Map.singleton e f
                                   _ -> Sum deleted
                            else Sum $ Map.insert b (fac + 1) a
                Nothing -> Sum $ Map.insert b 1 a
    where deleted = Map.delete b a
  a + Sum b = Sum b + a
  a + b = Sum (Map.singleton a 1) + b
  (-) = \x y -> x + (-y)
  -- the fromRational on the following line is needed to avoid a
  -- tricky infinite loop where -1 intepreted as (negate $ fromRational 1)
  negate = \x -> (fromRational $ -1) * x
  x * y | x == 0 = 0
        | y == 0 = 0
        | x == 1 = y
        | y == 1 = x
  Sum x * y | Just c <- isConstant y = pairs2sum $ map (f c) $ sum2pairs x
                      where f c (a,b) = (c*a, b)
  y * Sum x | Just c <- isConstant y = pairs2sum $ map (f c) $ sum2pairs x
                      where f c (a,b) = (c*a, b)
  Sum x * y | [(f,a)] <- sum2pairs x = pairs2sum [(f, a*y)]
  y * Sum x | [(f,a)] <- sum2pairs x = pairs2sum [(f, a*y)]
  Product a * Product b = puttogether a (product2pairs b)
      where puttogether x [] = map2product x
            puttogether x ((y,n):ys) =
                  case Map.lookup y x of
                    Just n' -> if n + n' == 0
                               then puttogether (Map.delete y x) ys
                               else puttogether (Map.insert y (n+n') x) ys
                    Nothing -> puttogether (Map.insert y n x) ys
  Product a * b = case Map.lookup b a of
                    Just n -> if n + 1 == 0
                             then map2product deleted
                             else map2product $ Map.insert b (n+1) a
                    Nothing -> map2product $ Map.insert b 1 a
    where deleted = Map.delete b a
  a * Product b = Product b * a
  a * b = Product (Map.singleton a 1) * b
  fromInteger = \x -> toExpression (fromInteger x :: Rational)
  abs = undefined
  signum = undefined

instance Type a => Fractional (Expression a) where
  x / y | Just yy <- isConstant y = x * toExpression (1/yy)
  x / Product y = x * Product (Map.map negate y)
  x / y = x * Product (Map.singleton y (-1))
  fromRational = toExpression

instance Type a => Floating (Expression a) where
  pi = toExpression (pi :: Double)
  exp = \x -> case x of 0 -> 1
                        _ -> Exp x
  log = \x -> case x of 1 -> 0
                        _ -> Log x
  sinh = \x -> case x of 0 -> 0
                         _ -> 0.5*(exp x - exp (-x))
  cosh = \x -> case x of 0 -> 1
                         _ -> 0.5*(exp x + exp (-x))
  sin = \x -> case x of 0 -> 0
                        _ -> Sin x
  cos = \x -> case x of
    0 -> 1
    _ -> Cos x
  a ** b | Just x <- isConstant a, Just y <- isConstant b = toExpression (x ** y)
  x ** y | y == 0 = 1
         | y == 1 = x
  (Sum x) ** c | Just n <- isConstant c,
                 [(f,y)] <- sum2pairs x = pairs2sum [(f**n, y ** c)]
  (Product x) ** c | Just n <- isConstant c = pairs2product $ map (p n) $ product2pairs x
                       where p n (e,n2) = (e,n2*n)
  x ** c | Just n <- isConstant c = pairs2product [(x,n)]
  x ** y = exp (y*log x)
  asin = undefined
  acos = undefined
  atan = undefined
  asinh = undefined
  acosh = undefined
  atanh = undefined

-- | grad takes the gradient of a scalar-valued expression with
-- respect to a particular realspace variable.

grad :: String -> Expression Scalar -> Expression RealSpace
grad v e = derive (r_var v) 1 e

data ToDo = DoK (Expression KSpace) | DoR (Expression RealSpace) | DoS (Expression Scalar) | DoNothing
            deriving (Eq, Show)

{-

How I find variables to eliminate:

1. Find a list of all variables.
2. For each variable, find a list of all (or many) subexpressions
   it occurs in.
3. For each subexpression, see if replacing that subexpression
   eliminates the variable, if so, do this.  Start with biggest
   subexpression?

-}

countVars :: Type a => Expression a -> Int
countVars s = length $ nub $ varList s
    where varList :: Type a => Expression a -> [String]
          varList (Expression e) | Same <- isKSpace (Expression e), Kx <- e = []
                                 | Same <- isKSpace (Expression e), Ky <- e = []
                                 | Same <- isKSpace (Expression e), Kz <- e = []
                                 | Same <- isRealSpace (Expression e), (IFFT e') <- e = varList e'
                                 | Same <- isKSpace (Expression e), (FFT e') <- e = varList e'
                                 | Same <- isScalar (Expression e), (Integrate e') <- e = varList e'
          varList (Var _ _ _ (Just e)) = varList e
          varList (Var _ c _ Nothing) = [c]
          varList (Sin e) = varList e
          varList (Cos e) = varList e
          varList (Log e) = varList e
          varList (Exp e) = varList e
          varList (Abs e) = varList e
          varList (Signum e) = varList e
          varList (Sum e) = concat $ map (varList . snd) (sum2pairs e)
          varList (Product e) = concat $ map (varList . fst) (product2pairs e)
          varList _ = []

removeAt :: [a] -> [Int] -> ([a], [a])
removeAt xs n = (map fst (filter (\x -> elem (snd x) n) list), map fst (filter (\x -> not $ elem (snd x) n) list))
    where list = consList 1 xs 
          consList m (p:ps) = [(p, m)] ++ (consList (m+1) ps)
          consList _ [] = []

hasFFT :: Type a => Expression a -> Bool
hasFFT (Expression e) 
    | Same <- isKSpace (Expression e), FFT _ <- e = True
    | Same <- isRealSpace (Expression e), IFFT _ <- e = True
hasFFT (Var _ _ _ (Just e)) = hasFFT e
hasFFT (Sum s) = or $ map (hasFFT . snd) (sum2pairs s)
hasFFT (Product p) = or $ map (hasFFT . fst) (product2pairs p)
hasFFT (Sin e) = hasFFT e
hasFFT (Cos e) = hasFFT e
hasFFT (Log e) = hasFFT e
hasFFT (Exp e) = hasFFT e
hasFFT (Abs e) = hasFFT e
hasFFT (Signum e) = hasFFT e
hasFFT _ = False

findToDo :: (Type a, Type b) => Expression a -> Expression b -> ToDo

findToDo everything (Expression e)
    | Same <- isKSpace (Expression e), FFT (Var _ _ _ Nothing) <- e = DoK (Expression e)
    | Same <- isKSpace (Expression e), FFT e' <- e = case findToDo everything e' of
                                                       DoNothing -> DoR $ e'
                                                       dothis -> dothis
    | Same <- isRealSpace (Expression e), IFFT (Var _ _ _ Nothing) <- e = DoR (Expression e)
    | Same <- isRealSpace (Expression e), IFFT e' <- e = case findToDo everything e' of
                                                           DoNothing -> DoK $ e'
                                                           dothis -> dothis
    | Same <- isScalar (Expression e), Integrate e' <- e = case findToDo everything e' of
                                                             DoNothing -> DoS $ Expression e
                                                             dothis -> dothis
    | otherwise = if hasFFT (Expression e)
                  then error ("FFT not detected in : " ++ show e)
                  else DoNothing -- error ("Missed Expression type in findToDo: " ++ show e)

findToDo _ (Sum s) | [] <- sum2pairs s = DoNothing
findToDo everything (Sum s) 
    | Same <- isRealSpace (Sum s), or $ map (simplifiable $ sum2pairs s) index, not $ hasFFT firstToDo = DoR firstToDo
    | Same <- isKSpace (Sum s), or $ map (simplifiable $ sum2pairs s) index, not $ hasFFT firstToDo = DoK firstToDo
    where index = reverse $ filter (\l -> length l /= 1) $ init $ take 1000 $ tail $ subsequences [1..length $ sum2pairs s] 
          simplifiable e n = (hasexpression (pairs2sum $ fst $ removeAt e n) (pairs2sum $ snd $ removeAt e n)) && ithelps e n
          oldnum = countVars everything
          ithelps :: Type a => [(Double, Expression a)] -> [Int] -> Bool
          ithelps e n | Same <- isRealSpace (pairs2sum $ fst $ removeAt e n)= countVars (substitute (pairs2sum $ fst $ removeAt e n) (r_var "WhatIfISubstituteThis") everything) < oldnum
                      | Same <- isKSpace (pairs2sum $ fst $ removeAt e n)= countVars (substitute (pairs2sum $ fst $ removeAt e n) (k_var "WhatIfISubstituteThis") everything) < oldnum
                      | otherwise = False
          firstToDo = pairs2sum $ fst $ removeAt (sum2pairs s) $ last $ filter (simplifiable $ sum2pairs s) index
findToDo everything (Sum s) = case filter (/= DoNothing) $ map sub $ sum2pairs s of
                                [] -> if hasFFT (Sum s)
                                      then error ("FFT not detected in : " ++ show (Sum s))
                                      else DoNothing
                                dothis:_ -> dothis
    where sub (_,e) = findToDo everything e

findToDo _ (Product s) | [] <- product2pairs s = DoNothing
findToDo everything (Product s) 
    | Same <- isRealSpace (Product s), or $ map (simplifiable $ product2pairs s) index, not $ hasFFT firstToDo = DoR firstToDo
    | Same <- isKSpace (Product s), or $ map (simplifiable $ product2pairs s) index, not $ hasFFT firstToDo = DoK firstToDo
    where index = reverse $ filter (\l -> length l /= 1) $ init $ take 1000 $ tail $ subsequences [1..length $ product2pairs s]
          simplifiable e n = (hasexpression (pairs2product $ fst $ removeAt e n) (pairs2product $ snd $ removeAt e n)) && ithelps e n
          oldnum = countVars everything
          ithelps :: Type a => [(Expression a, Double)] -> [Int] -> Bool
          ithelps e n | Same <- isRealSpace(pairs2product $ fst $ removeAt e n) = countVars (substitute (pairs2product $ fst $ removeAt e n) (r_var "WhatIfISubstituteThis") everything) < oldnum
                      | Same <- isKSpace(pairs2product $ fst $ removeAt e n) = countVars (substitute (pairs2product $ fst $ removeAt e n) (k_var "WhatIfISubstituteThis") everything) < oldnum
                      | otherwise = False
          firstToDo = pairs2product $ fst $ removeAt (product2pairs s) $ last $ filter (simplifiable $ product2pairs s) index

findToDo everything (Product p) = case filter (/= DoNothing) $ map sub $ product2pairs p of
                         [] -> if hasFFT (Product p)
                               then error ("FFT not detected in : " ++ show (Product p))
                               else DoNothing
                         dothis:_ -> dothis
    where sub (e, _) = findToDo everything e

findToDo x (Cos e) = findToDo x e
findToDo x (Sin e) = findToDo x e
findToDo x (Log e) = findToDo x e
findToDo x (Exp e) = findToDo x e
findToDo x (Abs e) = findToDo x e
findToDo x (Signum e) = findToDo x e
findToDo x (Var a b c (Just e)) =
  case findToDo x e of
    DoR e' -> case compareExpressions e e' of
                Same -> DoR $ Var a b c (Just e)
                Different -> DoR e'
    DoK e' -> case compareExpressions e e' of
                Same -> DoK $ Var a b c (Just e)
                Different -> DoK e'
    DoS e' -> case compareExpressions e e' of
                Same -> DoS $ Var a b c (Just e)
                Different -> DoS e'
    DoNothing -> DoNothing
--findToDo _ (Var _ _ _ (Just e)) = findToDo e
findToDo _ (Var _ _ _ Nothing) = DoNothing
findToDo everything (Scalar e) = findToDo everything e


isfft :: Expression KSpace -> Maybe (Expression RealSpace)
isfft (Expression (FFT e)) = Just e
isfft (Product p) = tofft 1 Nothing $ product2pairs p
  where tofft _ Nothing [] = Nothing
        tofft sc (Just e) [] = Just $ sc*e
        tofft sc Nothing ((Expression (FFT e),1):xs) = tofft sc (Just e) xs
        tofft sc myfft ((Scalar s,n):xs) = tofft (sc * Scalar s ** toExpression n) myfft xs
        tofft _ _ _ = Nothing
isfft (Var _ _ _ (Just e)) = isfft e
isfft _ = Nothing

isifft :: Expression RealSpace -> Maybe (Expression KSpace)
isifft (Expression (IFFT e)) = Just e
isifft (Product p) = tofft 1 Nothing $ product2pairs p
  where tofft _ Nothing [] = Nothing
        tofft sc (Just e) [] = Just $ sc*e
        tofft sc Nothing ((Expression (IFFT e),1):xs) = tofft sc (Just e) xs
        tofft sc myfft ((Scalar s,n):xs) = tofft (sc * Scalar s ** toExpression n) myfft xs
        tofft _ _ _ = Nothing
isifft (Var _ _ _ (Just e)) = isifft e
isifft _ = Nothing

joinFFTs :: Type a => Expression a -> Expression a
--joinFFTs e | cleanvars e /= e = joinFFTs $ cleanvars e
joinFFTs (Expression e)
  | Same <- isKSpace (Expression e), FFT e' <- e = fft (joinFFTs e')
  | Same <- isRealSpace (Expression e), IFFT e' <- e = ifft (joinFFTs e')
joinFFTs (Sum s) | Same <- isKSpace (Sum s) = joinup [] $ sum2pairs s
  where joinup tofft [] = fft $ factorandsum tofft
        joinup tofft ((f,x):xs) | Just e <- isfft x = joinup (toExpression f*e:tofft) xs
                                | otherwise = toExpression f*x + joinup tofft xs
joinFFTs (Sum s) | Same <- isRealSpace (Sum s) = joinup [] $ sum2pairs s
  where joinup tofft [] = ifft $ factorandsum tofft
        joinup tofft ((f,x):xs) | Just e <- isifft x = joinup (toExpression f*e:tofft) xs
                                | otherwise = toExpression f*x + joinup tofft xs
joinFFTs (Var a b c (Just e)) = Var a b c (Just (joinFFTs e))
joinFFTs e = e

simp2 :: Type a => Expression a -> ([Statement], Expression a)
simp2 = simp2helper (0 :: Int) [] -- . cleanvars
    where simp2helper n sts e = case findToDo e e of
                                  DoK ke -> simp2helper (n+1) (sts++[inike, setke]) e'
                                      where (v,vtex) = case ke of Var _ x@('k':_) t _ -> (x, t)
                                                                  Var _ x t _ -> ("k_" ++ x, t)
                                                                  _ -> ("k" ++ show n, "k_{" ++ show n ++ "}")
                                            inike = InitializeK $ Var (v++"[i]") v vtex Nothing
                                            setke = AssignK v ke
                                            e'  = substitute ke (Var (v++"[i]") v vtex Nothing) e
                                  DoR re -> simp2helper (n+1) (sts++[inire, setre]) e'
                                      where (v,vtex) = case re of Var _ x@('r':_) t _ -> (x,t)
                                                                  Var _ x t _ -> ("r_" ++ x,t)
                                                                  _ -> ("r" ++ show n, "r_{" ++ show n ++ "}")
                                            inire = InitializeR $ Var (v++"[i]") v vtex Nothing
                                            setre = AssignR v re
                                            e'  = substitute re (Var (v++"[i]") v vtex Nothing) e
                                  DoS re -> simp2helper (n+1) (sts++[inire, setre]) e'
                                      where (v,vtex) = case re of Var _ x@('s':_) t _ -> (x,t)
                                                                  Var _ x t _ -> ("s_" ++ x,t)
                                                                  _ -> ("s" ++ show n, "s_{" ++ show n ++ "}")
                                            inire = InitializeS $ Var v v vtex Nothing
                                            setre = AssignS v re
                                            e'  = substitute re (Var v v vtex Nothing) e
                                  DoNothing -> (sts, e)

factorandsum :: Type a => [Expression a] -> Expression a
factorandsum [] = 0
factorandsum (x:xs) = helper (getprodlist x) (x:xs)
  where helper :: Type a => [(Expression a, Double)] -> [Expression a] -> Expression a
        helper [] as = sum as
        helper ((f,n):fs) as = if n' /= 0
                               then (helper fs (map (/(f**(toExpression n'))) as)) * (f**(toExpression n'))
                               else helper fs as
          where n' = if n < 0 then if maximum (map (countup f) as) < 0
                                   then maximum (map (countup f) as)
                                   else 0
                              else if minimum (map (countup f) as) > 0
                                   then minimum (map (countup f) as)
                                   else 0
        countup f (Product y) = Map.findWithDefault 0 f y
        countup f (Sum a) | [(_,y)] <- sum2pairs a = countup f y
        countup f y | f == y = 1
        countup _ _ = 0
        getprodlist (Product xx) = product2pairs xx
        getprodlist (Sum a) | [(_,Product xx)] <- sum2pairs a = product2pairs xx
        getprodlist xx = [(xx,1)]

compareExpressions :: (Type a, Type b) => Expression a -> Expression b -> Same a b
compareExpressions x y | Same <- compareTypes x y, x == y = Same
compareExpressions _ _ = Different

compareTypes :: (Type a, Type b) => Expression a -> Expression b -> Same a b
compareTypes x y | Same <- isKSpace x, Same <- isKSpace y = Same
compareTypes x y | Same <- isRealSpace x, Same <- isRealSpace y = Same
compareTypes x y | Same <- isScalar x, Same <- isScalar y = Same
compareTypes _ _ = Different

derive :: (Type a, Type b) => Expression b -> Expression a -> Expression a -> Expression b
derive v dda e | Same <- compareExpressions v e = dda
derive v dda (Var _ _ _ (Just e)) = derive v dda e
derive _ _ (Var _ _ _ Nothing) = 0
derive v dda (Sum s) = factorandsum $ map dbythis $ sum2pairs s
  where dbythis (f,x) = toExpression f * derive v dda x
derive v dda (Product p) = factorandsum (map dbythis $ product2pairs p)
  where dbythis (x,n) = derive v (Product p*toExpression n*dda/x) x
derive _ _ (Scalar _) = 0 -- FIXME
derive v dda (Cos e) = derive v (-dda*sin e) e
derive v dda (Sin e) = derive v (dda*cos e) e
derive v dda (Exp e) = derive v (dda*exp e) e
derive v dda (Log e) = derive v (dda/e) e
derive _ _ (Abs _) = error "I didn't think we'd need abs"
derive _ _ (Signum _) = error "I didn't think we'd need signum"
derive v dda (Expression e) = derivativeHelper v dda e

hasexpression :: (Type a, Type b) => Expression a -> Expression b -> Bool
hasexpression x e | Same <- compareExpressions x e = True
hasexpression x (Expression v)
  | Same <- isKSpace (Expression v), FFT e <- v = hasexpression x e
  | Same <- isRealSpace (Expression v), IFFT e <- v = hasexpression x e
  | Same <- isScalar (Expression v), Integrate e <- v = hasexpression x e
  | Same <- isKSpace (Expression v), v == Kx || v == Ky || v == Kz || v == Delta = False
  | otherwise = error "Unhandled case in hasexpression"
hasexpression x@(Sum xs) e@(Sum es) -- check if x is a subexpression of e
  | Same <- compareTypes e x,
    filter (\(_, ex) -> ex == (snd $ head $ xspairs)) espairs /= [], 
    factorX <- fst $ head $ xspairs,
    factorE <- fst $ head $ filter (\(_, ex) -> ex == (snd $ head $ xspairs)) espairs,
    ratio <- factorE / factorX = 
         (and $ map (\term -> elem term $ espairs) (map (\(f, s) -> (f*ratio, s)) xspairs))
         || (or $ map (hasexpression x) (map snd espairs))
            where espairs = sum2pairs es
                  xspairs = sum2pairs xs
hasexpression x@(Product xs) e@(Product es)
  | Same <- compareTypes e x, filter (\(ex, _) -> ex == (fst $ head $ product2pairs xs)) (product2pairs es) /= [] = 
       (and $ map (\term -> elem term $ product2pairs es) (map (\(f, s) -> (f, s*((snd $ head $ filter (\(ex, _) -> ex == (fst $ head $ product2pairs xs)) (product2pairs es)) / factorX))) (product2pairs xs)))
       || (or $ map (hasexpression x) (map fst (product2pairs es)))
    where factorX = snd $ head $ product2pairs xs

hasexpression x (Sum s) = or $ map sub $ sum2pairs s
    where sub (_,e) = hasexpression x e
hasexpression x (Product p) = or $ map sub $ product2pairs p
    where sub (e, _) = hasexpression x e
hasexpression x (Cos e)   = hasexpression x e
hasexpression x (Sin e)   = hasexpression x e
hasexpression x (Log e)   = hasexpression x e
hasexpression x (Exp e)   = hasexpression x e
hasexpression x (Abs e)   = hasexpression x e
hasexpression x (Signum e) = hasexpression x e
hasexpression x (Var _ _ _ (Just e)) = hasexpression x e
hasexpression _ (Var _ _ _ Nothing) = False
hasexpression x (Scalar e) = hasexpression x e

substitute :: (Type a, Type b) => Expression a -> Expression a -> Expression b -> Expression b
substitute x y e | Same <- compareExpressions x e = y


substitute x@(Sum xs) y e@(Sum es) 
  | Same <- compareTypes x e, hasexpression x e = 
    let factorX = fst $ head $ xspairs
        factorE = fst $ head $ filter (\(_, ex) -> ex == (snd $ head $ xspairs)) espairs
        ratio = factorE / factorX
        in if (filter (\(_, ex) -> ex == (snd $ head $ xspairs)) espairs /= []) && (and $ map (\term -> elem term $ espairs) (multiplyFst ratio xspairs))
           then (pairs2sum $ map (subSnd x y) $ filter (\term -> not $ elem term (multiplyFst ratio xspairs)) espairs) + (y * (pairs2sum [(ratio, pairs2product [])]))
           else pairs2sum $ map (subSnd x y) espairs
  where subSnd m n (f, ex) = (f, substitute m n ex)
        espairs = sum2pairs es
        xspairs = sum2pairs xs
        multiplyFst factor pairs = map (\(a, b) -> (a*factor, b)) pairs

{-
substitute x@(Sum xs) y e@(Sum es) 
  | Same <- compareTypes x e, hasexpression x e =
    if (and $ map (\term -> elem term $ espairs) xspairs)
    then (pairs2sum $ map (subSnd x y) $ filter (\term -> not $ elem term xspairs) espairs) + y
    else pairs2sum $ map (subSnd x y) espairs
  where subSnd m n (f, ex) = (f, substitute m n ex)
        espairs = sum2pairs es
        xspairs = sum2pairs xs
-}
substitute x@(Product xs) y e@(Product es) 
  | Same <- compareTypes x e, hasexpression x e = 
    if (and $ map (\term -> elem term $ product2pairs es) (product2pairs xs))
    then (pairs2product $ map (subFst x y) $ filter (\term -> not $ elem term $ product2pairs xs) (product2pairs es)) * y
    else pairs2product $ map (subFst x y) (product2pairs es)
  where subFst m n (ex, f) = (substitute m n ex, f)
substitute x y (Expression v)
  | Same <- isKSpace (Expression v), FFT e <- v = Expression $ FFT (substitute x y e)
  | Same <- isRealSpace (Expression v), IFFT e <- v = Expression $ IFFT (substitute x y e)
  | Same <- isScalar (Expression v), Integrate e <- v = Expression $ Integrate (substitute x y e)
  | Same <- isKSpace (Expression v), v == Kx || v == Ky || v == Kz || v == Delta = Expression v
  | otherwise = error $ "unhandled case in substitute: " ++ show v
substitute x y (Sum s) = pairs2sum $ map sub $ sum2pairs s
    where sub (f,e) = (f, substitute x y e)
substitute x y (Product p) = pairs2product $ map sub $ product2pairs p
    where sub (e, n) = (substitute x y e, n)
substitute x y (Cos e)   = cos $ substitute x y e
substitute x y (Sin e)   = sin $ substitute x y e
substitute x y (Log e)   = log $ substitute x y e
substitute x y (Exp e)   = exp $ substitute x y e
substitute x y (Abs e)   = abs $ substitute x y e
substitute x y (Signum e) = signum $ substitute x y e
substitute x y (Var a b c (Just e)) = Var a b c (Just $ substitute x y e)
substitute _ _ v@(Var _ _ _ Nothing) = v
substitute x y (Scalar e) = Scalar (substitute x y e)

substituteS :: Type a => Expression a -> Expression a -> Statement -> Statement
substituteS x y (AssignR s e) = AssignR s (substitute x y e)
substituteS x y (AssignK s e) = AssignK s (substitute x y e)
substituteS x y (AssignS s e) = AssignS s (substitute x y e)
substituteS x y (InitializeR e) = InitializeR (substitute x y e)
substituteS x y (InitializeK e) = InitializeK (substitute x y e)
substituteS x y (InitializeS e) = InitializeS (substitute x y e)
substituteS x y (FreeR e) = FreeR (substitute x y e)
substituteS x y (FreeK e) = FreeK (substitute x y e)

\end{code}

The statement type allows us to string together a sequence of
definitions and assignments.

\begin{code}
data Statement = AssignR String (Expression RealSpace)
               | AssignK String (Expression KSpace)
               | AssignS String (Expression Scalar)
               | InitializeR (Expression RealSpace)
               | InitializeK (Expression KSpace)
               | InitializeS (Expression Scalar)
               | FreeR (Expression RealSpace)
               | FreeK (Expression KSpace)

instance Show Statement where
  showsPrec = showsS
showsS :: Int -> Statement -> ShowS
showsS _ (AssignR x y) = showString x . showString " := " . showsPrec 0 y
showsS _ (AssignK x y) = showString x . showString " := " . showsPrec 0 y
showsS _ (AssignS x y) = showString x . showString " := " . showsPrec 0 y
showsS _ (InitializeR x) = showString "Initialize " . showsPrec 0 x
showsS _ (InitializeK x) = showString "Initialize " . showsPrec 0 x
showsS _ (InitializeS x) = showString "Initialize " . showsPrec 0 x
showsS _ (FreeR x) = showString "Free " . showsPrec 0 x
showsS _ (FreeK x) = showString "Free " . showsPrec 0 x

instance Code Statement where
  codePrec = codeS
  latexPrec = latexS
codeS :: Int -> Statement -> ShowS
codeS _ (AssignR x y) = showString (prefix y) . codePrec 0 (r_var x) . showString " = " . codePrec 0 y . showString ";" . showString (postfix y)
codeS _ (AssignK x y) = showString (prefix y) . codePrec 0 (k_var x) . showString " = " . codePrec 0 y . showString ";" . showString (postfix y)
codeS _ (AssignS x y) = showString (prefix y) . codePrec 0 (s_var x :: Expression Scalar) . showString " = " . codePrec 0 y . showString ";" . showString (postfix y)
codeS _ (InitializeR e) = showString (initialize e)
codeS _ (InitializeK e) = showString (initialize e)
codeS _ (InitializeS e) = showString (initialize e)
codeS _ (FreeR e) = showString (free e)
codeS _ (FreeK e) = showString (free e)

latexS :: Int -> Statement -> ShowS
latexS _ (AssignR x y) = latexPrec 0 (r_var x) . showString " = " . latexPrec 0 (cleanvars y)
latexS _ (AssignK x y) = latexPrec 0 (k_var x) . showString " = " . latexPrec 0 (cleanvars y)
latexS _ (AssignS x y) = latexPrec 0 (s_var x :: Expression Scalar) . showString " = " . latexPrec 0 (cleanvars y)
latexS _ (InitializeR e) = showString (initialize e)
latexS _ (InitializeK e) = showString (initialize e)
latexS _ (InitializeS e) = showString (initialize e)
latexS _ (FreeR e) = showString (free e)
latexS _ (FreeK e) = showString (free e)

latexStatements :: [Statement] -> String
latexStatements x = unlines $ map (\e -> "\n\\begin{dmath}\n" ++ latex e ++ "\n\\end{dmath}") x

latexSimp :: (Type a) => Expression a -> String
latexSimp e = "\\documentclass{article}\n\\usepackage{amsmath}\n\\usepackage{breqn}\n\\begin{document}\n\n" ++
              latexStatements (f sts) ++ "\n\\begin{dmath}\n" ++ latex e' {-(cleanvars e')-} ++ "\n\\end{dmath}" ++
              "\n\\end{document}"
    where f ((InitializeR _):xs) = f xs
          f ((InitializeK _):xs) = f xs
          f ((InitializeS _):xs) = f xs
          f ((FreeR _):xs) = f xs
          f ((FreeK _):xs) = f xs
          f (x:xs) = x:(f xs)
          f [] = []
          (sts,e') = simp2 e


checkDup :: [Statement] -> [Statement]
checkDup = nubBy sameInit
    where sameInit (InitializeR x) (InitializeR y) = x == y
          sameInit (InitializeK x) (InitializeK y) = x == y
--          sameInit (AssignR x y) (AssignR x' y') = x == x' && y == y'
--          sameInit (AssignK x y) (AssignK x' y') = x == x' && y == y'
          sameInit _ _ = False

freeVectors :: [Statement] -> [Statement]
freeVectors s = reverse $ freeHelper (vecInMem s) (reverse s)
    where freeHelper v (x@(AssignR _ xn):xs) | or $ map check v = ns ++ [x] ++ freeHelper (v \\ fe) xs
              where fe = filter check v
                    ns = map freeme fe
                    check n = (hasexpression (r_var n) xn) || (hasexpression (k_var n) xn)
          freeHelper v (x@(AssignK _ xn):xs) | or $ map check v = ns ++ [x] ++ freeHelper (v \\ fe) xs
              where fe = filter check v
                    ns = map freeme fe
                    check n = (hasexpression (r_var n) xn) || (hasexpression (k_var n) xn)
          freeHelper v (x:xs) = [x] ++ freeHelper v xs
          freeHelper _ [] = []
          freeme v@('r':_) = FreeR $ r_var v
          freeme v@('k':_) = FreeK $ k_var v
          freeme _ = error "bad variable name!?!@"

vecInMem :: [Statement] -> [String]
vecInMem s = filter isFreeVec $ map ini s
     where ini (InitializeR (Var _ x _ Nothing)) = x
           ini (InitializeK (Var _ x _ Nothing)) = x
           ini _ = ""
           fre (FreeR (Var _ x _ Nothing)) = x
           fre (FreeK (Var _ x _ Nothing)) = x
           fre _ = ""
           isFreeVec x = not $ elem x $ map fre s

reuseVar :: [Statement] -> [Statement]
reuseVar ((InitializeR (Var _ ivar _ Nothing)):(AssignR n e):(FreeR (Var _ fvar _ Nothing)):xs)
    | ivar == n = (AssignR fvar e) : (reuseVar $ map (substituteS (r_var ivar) (r_var fvar)) xs)
    | otherwise = error "RS initialize error: "
reuseVar ((InitializeK (Var _ ivar _ Nothing)):(AssignK n e):(FreeK (Var _ fvar _ Nothing)):xs)
    | ivar == n = (AssignK fvar e) : (reuseVar $ map (substituteS (k_var ivar) (k_var fvar)) xs)
    | otherwise = error "initialize error"
reuseVar (x:xs) = x : (reuseVar xs)
reuseVar [] = []

\end{code}

\begin{code}

functionCode :: String -> String -> [(String, String)] -> String -> String
functionCode "" "" [] "" = ""
functionCode "" "" (x:xs) "" = if xs == [] 
                               then fst x ++ " " ++ snd x
                               else fst x ++ " " ++ snd x ++ ", " ++ functionCode "" "" xs ""
functionCode n t a b = t ++ " " ++ n ++ "(" ++ functionCode "" "" a "" ++ ") const {\n" ++ b ++ "\n}\n"



classCode :: Expression RealSpace -> [String] -> String -> String
classCode e arg n = "class " ++ n ++ " : public FunctionalInterface {\npublic:\n" ++ n ++ codeA arg ++ "  {\n\thave_integral = true;\n}\n" ++
                functionCode "I_have_analytic_grad" "bool" [] "\treturn false;" ++
                functionCode "integral" "double" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n\tdouble output=0;\n\tconst double dr = pow(gd.fineLat.volume(), 1.0/3); assert(dr);\n" ++ codeStatements codeIntegrate ++ "\t// " ++ show (countFFT codeIntegrate) ++ " Fourier transform used.\n\t// " ++ show (peakMem codeIntegrate) ++ "\n\treturn output;\n" ) ++
                functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3); assert(dr);\n" ++ codeStatements codeVTransform  ++ "\t// " ++ show (countFFT codeVTransform) ++ " Fourier transform used.\n\t// " ++ show (peakMem codeVTransform)  ++ "\n\treturn output;\n")  ++
                functionCode "transform" "double" [("double", "kT"), ("double", "x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(x==x); // to avoid an unused parameter error\n\t" ++ codeStatements codeDTransform ++ "\n\treturn output;\n") ++
                functionCode "append_to_name" "bool" [("", "std::string")] "\treturn false;" ++
                functionCode "derive" "double" [("double", "kT"), ("double", "x")] ("\tassert(kT==kT);\n\tassert(x==x);\n" ++ codeStatements codeDerive ++ "\n\treturn output;\n") ++
                functionCode "d_by_dT" "double" [("double", ""), ("double", "")] "\tassert(0); // fail\n\treturn 0;\n" ++
                functionCode "derive_homogeneous" "Expression" [("const Expression &", "")] "\tassert(0); // fail\n\treturn Expression(0);\n" ++
                functionCode "grad" "Functional" [("const Functional", "&ingrad"), ("const Functional", "&x"), ("bool", "")] "\tassert(&ingrad==&ingrad);\n\tassert(&x==&x);\n\treturn ingrad;" ++
                functionCode "grad_T" "Functional" [("const Functional", "&ingradT")] "\tassert(&ingradT==&ingradT);\n\treturn ingradT;" ++
                functionCode "grad" "void" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x"), ("const VectorXd", "&ingrad"), ("VectorXd", "*outgrad"), ("VectorXd", "*outpgrad")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n\tassert(outpgrad==outpgrad);\n\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3); assert(dr);\n" ++ codeStatements codeGrad ++ "\t// " ++ show (countFFT codeGrad) ++ " Fourier transform used.\n\t// " ++ show (peakMem codeGrad) ++ "\n") ++
                functionCode "printme" "Expression" [("const Expression", "&x")] ("\treturn funexpr(\"" ++ n ++ "()\")(x);") ++
                functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")] "\tFunctionalInterface::print_summary(prefix, energy, name);" ++
                "private:\n"++ codeArgInit arg  ++"}; // End of " ++ n ++ " class\n\t// Total " ++ (show $ (countFFT codeIntegrate + countFFT codeVTransform + countFFT codeGrad)) ++ " Fourier transform used.\n\t// peak memory used: " ++ (show $ maximum $ map peakMem [codeIntegrate, codeVTransform, codeGrad])
    where
      codeIntegrate = freeVectors (st ++ [AssignS "output" e'])
          where (st, e') = simp2 (joinFFTs $ integrate e)
      codeVTransform = freeVectors (InitializeR e: st ++ [AssignR "output" e'])
          where (st, e') = simp2 $ joinFFTs e
      codeDTransform = freeVectors [InitializeS (makeHomogeneous e), AssignS "output" (makeHomogeneous e)]
      codeDerive = freeVectors [InitializeS (makeHomogeneous (derive (r_var "x") 1 e)), AssignS "output" (makeHomogeneous (derive (r_var "x") 1 e))]
      codeGrad = freeVectors (st ++ [AssignR "(*outgrad)" (r_var "(*outgrad)" + e')])
          where (st, e') = simp2 (joinFFTs $ derive (r_var "x") (r_var "ingrad") $ cleanvars e ) -- why cleanvars here ?
      codeA [] = "()"
      codeA a = "(" ++ foldl1 (\x y -> x ++ ", " ++ y ) (map (\x -> "double " ++ x ++ "_arg") a) ++ ") : " ++ foldl1 (\x y -> x ++ ", " ++ y) (map (\x -> x ++ "(" ++ x ++ "_arg)") a)
      codeArgInit [] = ""
      codeArgInit a = unlines $ map (\x -> "\tdouble " ++ x ++ ";") a

generateHeader :: Expression RealSpace -> [String] -> String -> String
generateHeader e arg n = "// -*- mode: C++; -*-\n\n#include \"MinimalFunctionals.h\"\n#include \"utilities.h\"\n#include \"handymath.h\"\n\n" ++ 
                     classCode e arg (n ++ "_type") ++
                     "\n\nFunctional " ++ n ++"(" ++ codeA arg ++ ") {\n\treturn Functional(new " ++ n ++ "_type(" ++ codeA' arg ++ "), \"" ++ n ++ "\");\n}\n"
    where codeA [] = ""
          codeA a = foldl1 (\x y -> x ++ ", " ++ y ) (map ("double " ++) a)
          codeA' [] = ""
          codeA' a = foldl1 (\x y -> x ++ ", " ++ y ) a



\end{code}
