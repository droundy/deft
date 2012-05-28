\begin{code}
{-# LANGUAGE GADTs, PatternGuards #-}

module Expression (
                   RealSpace(..), r_var,
                   KSpace(..), k_var, kx, ky, kz, k, ksqr, setkzero,
                   Scalar(..),
                   fft, ifft, integrate, grad, derive,
                   Expression(..), joinFFTs, (===), var,
                   Type(..), Code(..), Same(..), IsTemp(..),
                   makeHomogeneous, isConstant,
                   setZero, cleanvars, factorandsum,
                   sum2pairs, pairs2sum,
                   product2pairs, pairs2product, product2denominator,
                   hasK, hasFFT, hasexpression,
                   countexpression, substitute, 
                   compareExpressions, countVars)
    where

import Debug.Trace

import qualified Data.Map as Map
import qualified Data.Set as Set
import Data.List ( nub )
\end{code}

The \verb!RealSpaceField! data type describes a field in real space.

\begin{code}
data RealSpace = IFFT (Expression KSpace)
               deriving ( Eq, Ord, Show )
data KSpace = Delta | -- handy for FFT of homogeneous systems
              Kx | Ky | Kz |
              SetKZeroValue (Expression KSpace) (Expression KSpace) |
              FFT (Expression RealSpace)
            deriving ( Eq, Ord, Show )
data Scalar = Integrate (Expression RealSpace)
            deriving ( Eq, Ord, Show )

kinversion :: Expression KSpace -> Expression KSpace
kinversion (Var _ _ _ _ (Just e)) = kinversion e
kinversion e@(Var _ _ _ _ Nothing) = e
kinversion (Scalar e) = Scalar e
kinversion (Cos e) = cos (kinversion e)
kinversion (Sin e) = sin (kinversion e)
kinversion (Exp e) = exp (kinversion e)
kinversion (Log e) = log (kinversion e)
kinversion (Abs e) = abs (kinversion e)
kinversion (Signum e) = signum (kinversion e)
kinversion (Product p _) = product $ map ff $ product2pairs p
  where ff (e,n) = (kinversion e) ** toExpression n
kinversion (Sum s _) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, kinversion y)
kinversion (Expression Kx) = -kx
kinversion (Expression Ky) = -ky
kinversion (Expression Kz) = -kz
kinversion (Expression x) = Expression x

instance Code RealSpace where
  codePrec _ (IFFT (Var _ _ ksp _ Nothing)) = showString ("ifft(gd, " ++ksp++ ")")
  codePrec _ (IFFT ke) = showString "ifft(gd, " . codePrec 0 ke . showString ")"
  latexPrec _ (IFFT ke) = showString "\\text{ifft}\\left(" . latexPrec 0 ke . showString "\\right)"
instance Type RealSpace where
  isRealSpace _ = Same
  derivativeHelper v ddr (IFFT ke) = derive v (fft ddr) (kinversion ke)
  zeroHelper v (IFFT ke) = ifft (setZero v ke)
  prefix e = case isifft e of Just _ -> ""
                              Nothing -> "for (int i=0; i<gd.NxNyNz; i++) {\n\t\t"
  postfix e = case isifft e of Just _ -> ""
                               Nothing -> "\t}\n"
  codeStatementHelper (Var _ _ a _ _) op (Expression (IFFT (Var _ _ v _ Nothing))) = a ++ op ++ "ifft(gd, " ++ v ++ ");\n"
  codeStatementHelper _ _ (Expression (IFFT e)) = error ("It is a bug to generate code for a non-var input to ifft\n"++ latex e)
  codeStatementHelper a op (Var _ _ _ _ (Just e)) = codeStatementHelper a op e
  codeStatementHelper a op e = prefix e ++ code a ++ op ++ code e ++ ";\n" ++ postfix e
  initialize (Var IsTemp _ x _ Nothing) = "VectorXd " ++ x ++ "(gd.NxNyNz);"
  initialize _ = error "VectorXd output(gd.NxNyNz);"
  free (Var IsTemp _ x _ Nothing) = x ++ ".resize(0); // Realspace"
  free e = error $ trace "free error" ("free error " ++ show e)
  toScalar (IFFT ke) = makeHomogeneous ke

instance Code KSpace where
  codePrec _ Kx = showString "k_i[0]"
  codePrec _ Ky = showString "k_i[1]"
  codePrec _ Kz = showString "k_i[2]"
  codePrec _ Delta = showString "delta(k?)"
  codePrec p (SetKZeroValue _ val) = codePrec p val
  codePrec _ (FFT r) = showString "fft(gd, " . codePrec 0 (makeHomogeneous r) . showString ")"
  latexPrec _ Kx = showString "k_{x}"
  latexPrec _ Ky = showString "k_{y}"
  latexPrec _ Kz = showString "k_{z}"
  latexPrec _ Delta = showString "\\delta(k)"
  latexPrec p (SetKZeroValue _ val) = latexPrec p val
  latexPrec _ (FFT r) = showString "\\text{fft}\\left(" . latexPrec 0 r . showString "\\right)"
instance Type KSpace where
  isKSpace _ = Same
  derivativeHelper v ddk (FFT r) = derive v (ifft ddk) r
  derivativeHelper v ddk (SetKZeroValue _ e) = derive v (setkzero 0 ddk) e -- FIXME: how best to handle k=0 derivative?
  derivativeHelper _ _ Kx = 0
  derivativeHelper _ _ Ky = 0
  derivativeHelper _ _ Kz = 0
  derivativeHelper _ _ Delta = 0
  zeroHelper v (FFT r) = fft (setZero v r)
  zeroHelper _ Kx = Expression Kx
  zeroHelper _ Ky = Expression Ky
  zeroHelper _ Kz = Expression Kz
  zeroHelper _ Delta = Expression Delta
  zeroHelper v e@(SetKZeroValue val _) | Same <- isKSpace v, v == kz = val -- slightly hokey... assuming that if we set kz = 0 then we set kx and ky = 0
                                       | otherwise = Expression e
  prefix e = case isfft e of Just _ -> ""
                             Nothing -> "for (int i=1; i<gd.NxNyNzOver2; i++) {\n\t\tconst int z = i % gd.NzOver2;\n\t\tconst int n = (i-z)/gd.NzOver2;\n\t\tconst int y = n % gd.Ny;\n\t\tconst int xa = (n-y)/gd.Ny;\n\t\tconst RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);\n\t\tconst Reciprocal k_i = gd.Lat.toReciprocal(rvec);\n\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3); assert(dr);\n"
  postfix e = case isfft e of Just _ -> ""
                              Nothing -> "\n\t}\n"
  codeStatementHelper (Var _ _ a _ _) op (Expression (FFT (Var _ _ v _ Nothing))) = a ++ op ++ "fft(gd, " ++ v ++ ");\n"
  codeStatementHelper _ _ (Expression (FFT _)) = error "It is a bug to generate code for a non-var input to fft"
  codeStatementHelper a op (Var _ _ _ _ (Just e)) = codeStatementHelper a op e
  codeStatementHelper (Var _ _ a _ _) op e =
          if k0code == "0"
          then a ++ "[0]" ++ op ++ "0;\n\t" ++ prefix e ++ "\t\t" ++ a ++ "[i]"  ++ op ++ code e ++ ";" ++ postfix e
          else "{\n\t\tconst int i = 0;\n\t\tconst Reciprocal k_i = Reciprocal(0,0,0);\n\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3); assert(dr);\n\t\t" ++ a ++ "[0]" ++ op ++ code (setZero kz (setZero ky (setZero kx e))) ++ ";\n\t}\n\t" ++ prefix e ++ "\t\t" ++ a ++ "[i]"  ++ op ++ code e ++ ";" ++ postfix e
      where k0code = code (setZero kz (setZero ky (setZero kx e)))
  codeStatementHelper _ _ _ = error "Illegal input to codeStatementHelper for kspace"
  initialize (Var IsTemp _ x _ Nothing) = "VectorXcd " ++ x ++ "(gd.NxNyNzOver2);"
  initialize _ = error "VectorXcd output(gd.NxNyNzOver2);"
  free (Var IsTemp _ x _ Nothing) = x ++ ".resize(0); // KSpace"
  free _ = error "free error"
  toScalar Delta = 1
  toScalar Kx = s_var "_kx"
  toScalar Ky = 0
  toScalar Kz = 0
  toScalar (SetKZeroValue val _) = makeHomogeneous val
  toScalar (FFT e) = makeHomogeneous e

mapExpression :: (Type a, Type b) => (a -> Expression b) -> Expression a -> Expression b
mapExpression f (Var t a b c (Just e)) = Var t a b c $ Just $ mapExpression f e
mapExpression _ (Var tt c v t Nothing) = Var tt c v t Nothing
mapExpression _ (Scalar e) = Scalar e
mapExpression f (Cos e) = cos (mapExpression f e)
mapExpression f (Sin e) = sin (mapExpression f e)
mapExpression f (Exp e) = exp (mapExpression f e)
mapExpression f (Log e) = log (mapExpression f e)
mapExpression f (Abs e) = abs (mapExpression f e)
mapExpression f (Signum e) = signum (mapExpression f e)
mapExpression f (Product p _) = product $ map ff $ product2pairs p
  where ff (e,n) = (mapExpression f e) ** toExpression n
mapExpression f (Sum s _) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, mapExpression f y)
mapExpression f (Expression x) = f x

cleanvars :: Type a => Expression a -> Expression a
cleanvars (Var _ _ _ _ (Just e)) = cleanvars e
cleanvars (Var tt c v t Nothing) = Var tt c v t Nothing
cleanvars (Scalar e) = Scalar (cleanvars e)
cleanvars (Cos e) = cos (cleanvars e)
cleanvars (Sin e) = sin (cleanvars e)
cleanvars (Exp e) = exp (cleanvars e)
cleanvars (Log e) = log (cleanvars e)
cleanvars (Abs e) = abs (cleanvars e)
cleanvars (Signum e) = signum (cleanvars e)
cleanvars (Product p _) = product $ map ff $ product2pairs p
  where ff (e,n) = (cleanvars e) ** toExpression n
cleanvars (Sum s _) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, cleanvars y)
cleanvars (Expression e)
  | Same <- isKSpace (Expression e), FFT e' <- e = fft (cleanvars e')
  | Same <- isRealSpace (Expression e), IFFT e' <- e = ifft (cleanvars e')
  | Same <- isScalar (Expression e), Integrate e' <- e = integrate (cleanvars e')
  | otherwise = Expression e

isEven :: (Type a, Type b) => Expression b -> Expression a -> Double
isEven v e | Same <- compareExpressions v e = -1
isEven v (Var _ _ _ _ (Just e)) = isEven v e
isEven _ (Var _ _ _ _ Nothing) = 1
isEven v (Scalar e) = isEven v e
isEven _ (Cos _) = 1
isEven v (Sin e) = isEven v e
isEven v (Exp e) = if isEven v e == 1 then 1 else 0
isEven v (Log e) = if isEven v e == 1 then 1 else 0
isEven _ (Abs _) = 1
isEven v (Signum e) = isEven v e
isEven v (Product p _) = product $ map ie $ product2pairs p
  where ie (x,n) = isEven v x ** n
isEven _ (Sum s i) | Sum s i == 0 = 1 -- ???
isEven v (Sum s _) = ie (isEven v x) xs
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
setZero v (Var t a b c (Just e)) = Var t a b c (Just $ setZero v e)
setZero _ e@(Var _ _ _ _ Nothing) = e
setZero v (Scalar e) = Scalar (setZero v e)
setZero v (Cos e) = cos (setZero v e)
setZero v (Sin e) = sin (setZero v e)
setZero v (Exp e) = exp (setZero v e)
setZero v (Log e) = log (setZero v e)
setZero v (Abs e) = abs (setZero v e)
setZero v (Signum e) = signum (setZero v e)
setZero v (Product p _) | product2denominator p == 1 = product $ map ff $ product2pairs p
  where ff (e,n) = (setZero v e) ** toExpression n
setZero v (Product p i) =
  if isEven v (Product p i) == -1
  then 0
  else
    if zd /= 0
    then zn / zd
    else if zn /= 0
         then error ("L'Hopital's rule failure: " ++ latex n ++ "\n /\n  " ++ latex d ++ "\n\n\n" 
                     ++ latex (Product p i) ++ "\n\n\n" ++ latex zn)
         else case compareTypes v n of
                Same -> setZero v (dtop / dbot)
                  where dtop = derive v 1 n
                        dbot = derive v 1 d
                Different -> error "oops, failure in setZero that makes no sense."
  where d = product2denominator p
        n = product $ product2numerator p
        zn = setZero v n
        zd = setZero v d
setZero _ (Sum s i) | Sum s i == 0 = 0
setZero v (Sum s _) = out
  where sz (f,x) = (f, setZero v x)
        out = pairs2sum $ map sz $ sum2pairs s
setZero v (Expression x) = zeroHelper v x

instance Code Scalar where
  codePrec _ (Integrate r) = showString "integrate(" . codePrec 0 r . showString ")"
  latexPrec _ (Integrate r) = showString "integrate(" . latexPrec 0 r . showString ")"
instance Type Scalar where
  s_var ("complex(0,1)") = Var CannotBeFreed "complex(0,1)" "complex(0,1)" "i" Nothing
  s_var v@['d',_] = Var CannotBeFreed v v v Nothing -- for differentials
  s_var vv@(a:v@(_:_)) = Var CannotBeFreed vv vv (a : '_' : '{' : v ++ "}") Nothing
  s_var v = Var CannotBeFreed v v v Nothing
  isScalar _ = Same
  derivativeHelper v dds (Integrate e) = derive v (Scalar dds*s_var "dV") e
  zeroHelper v (Integrate e) = integrate (setZero v e)
  codeStatementHelper a " = " (Expression (Integrate e)) =
    code a ++ " = 0;\n\tfor (int i=0; i<gd.NxNyNz; i++) {\n\t\t" ++
    code a ++ " += " ++ code (e * s_var "gd.dvolume") ++
    ";\n\t}\n"
  codeStatementHelper _ op (Expression (Integrate _)) = error ("Haven't implemented "++op++" for integrate...")
  codeStatementHelper a op e = code a ++ op ++ code e ++ ";"
  prefix (Expression (Integrate _))  = "for (int i=0; i<gd.NxNyNz; i++) {    "
  prefix _ = ""
  postfix (Expression (Integrate _)) = "\n}\n"
  postfix _ = ""
  initialize (Expression (Integrate _)) = "double output = 0;\n"
  initialize (Var IsTemp _ x _ Nothing) = "double " ++ x ++ " = 0;\n"
  initialize _ = error "double output = 0;\n"
  toScalar (Integrate r) = makeHomogeneous r
  fromScalar = id

r_var :: String -> Expression RealSpace
r_var v@([_]) = Var CannotBeFreed (v++"[i]") v v Nothing
r_var v | take 5 v == "rtemp" = Var IsTemp (v++"[i]") ('r':v) v Nothing
r_var v | take 4 v == "temp" = Var IsTemp ("r"++v++"[i]") ('r':v) v Nothing
r_var (a:v) = Var CannotBeFreed (a:v++"[i]") (a:v) (a:'_':'{':v++"}") Nothing
r_var "" = error "r_var needs non-empty string"

k_var :: String -> Expression KSpace
k_var v@([_]) = Var CannotBeFreed (v++"[i]") v v Nothing
k_var v | take 5 v == "ktemp" = Var IsTemp (v++"[i]") ('r':v) v Nothing
k_var v | take 4 v == "temp" = Var IsTemp ("k"++v++"[i]") ('r':v) v Nothing
k_var (a:v) = Var CannotBeFreed (a:v++"[i]") (a:v) (a:'_':'{':v++"}") Nothing
k_var "" = error "r_var needs non-empty string"

infix 4 ===

(===) :: Type a => String -> Expression a -> Expression a
--_ === e = e
v@(a:r@(_:_)) === e = Var CannotBeFreed c v ltx (Just e)
  where ltx = a : "_{"++r++"}"
        c = (case isScalar e of
              Same -> v
              Different -> v ++ "[i]") :: String
v === e = Var CannotBeFreed c v v (Just e)
  where c = (case isScalar e of
              Same -> v
              Different -> v ++ "[i]") :: String

var :: Type a => String -> String -> Expression a -> Expression a
var v ltx e = Var CannotBeFreed c v ltx (Just e)
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

setkzero :: Expression KSpace -> Expression KSpace -> Expression KSpace
setkzero zeroval otherval = Expression $ SetKZeroValue zeroval otherval

integrate :: Type a => Expression RealSpace -> Expression a
integrate (Sum s _) = pairs2sum $ map i $ sum2pairs s
  where i (f, e) = (f, integrate e)
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
-- IsTemp is a boolean type for telling whether a variable is temporary or not.
data IsTemp = IsTemp | CannotBeFreed
            deriving (Eq, Ord, Show)

data Expression a = Scalar (Expression Scalar) |
                    Var IsTemp String String String (Maybe (Expression a)) | -- A variable with a possible value
                    Expression a |
                    Cos (Expression a) |
                    Sin (Expression a) |
                    Exp (Expression a) |
                    Log (Expression a) |
                    Abs (Expression a) |
                    Signum (Expression a) |
                    Product (Map.Map (Expression a) Double) (Set.Set String) |
                    Sum (Map.Map (Expression a) Double) (Set.Set String)
              deriving (Eq, Ord, Show)

sum2pairs :: Type a => Map.Map (Expression a) Double -> [(Double, Expression a)]
sum2pairs s = map rev $ Map.assocs s
  where rev (a,b) = (b,a)

pairs2sum :: Type a => [(Double, Expression a)] -> Expression a
pairs2sum s = helper $ filter ((/= 0) . snd) $ filter ((/= 0) . fst) s
  where helper [] = 0
        helper [(1,e)] = e
        helper es = map2sum $ fl (Map.empty) es
        fl a [] = a
        fl a ((f,Sum s' _):xs) = fl a (map mulf (sum2pairs s') ++ xs)
          where mulf (ff,e) = (ff*f, e)
        fl a ((f,x):xs) = case Map.lookup x a of
                          Just f' -> if f + f' == 0
                                     then fl (Map.delete x a) xs
                                     else fl (Map.insert x (f+f') a) xs
                          Nothing -> fl (Map.insert x f a) xs

map2sum :: Type a => Map.Map (Expression a) Double -> Expression a
map2sum s | Map.size s == 1 =
  case sum2pairs s of [(1,e)] -> e
                      _ -> Sum s (Set.unions $ map varSet $ Map.keys s)
map2sum s = Sum s (Set.unions $ map varSet $ Map.keys s)

product2pairs :: Type a => Map.Map (Expression a) Double -> [(Expression a, Double)]
product2pairs s = Map.assocs s

map2product :: Type a => Map.Map (Expression a) Double -> Expression a
map2product p | Map.size p == 1 =
  case product2pairs p of [(e,1)] -> e
                          _ -> Product p (Set.unions $ map varSet $ Map.keys p)
map2product p = helper 1 (Map.empty) $ product2pairs p
  where helper 1 a [] = Product a (vl a)
        helper f a [] = Sum (Map.singleton (Product a i) f) i
          where i = vl a
        helper f a ((Sum x _,n):xs) | [(f',x')] <- sum2pairs x = helper (f*f'**n) a ((x',n):xs)
        helper f a ((x,n):xs) = helper f (Map.insert x n a) xs
        vl = Set.unions . map varSet . Map.keys

pairs2product :: Type a => [(Expression a, Double)] -> Expression a
pairs2product = map2product . fl (Map.empty)
  where fl a [] = a
        fl a ((_,0):xs) = fl a xs
        fl a ((Product p _, n):xs) = fl a (map tonpower (product2pairs p) ++ xs)
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
codeE _ (Var _ c _ _ Nothing) = showString c
codeE p (Var _ _ _ _ (Just e)) = codeE p e
codeE p (Scalar x) = codePrec p x
codeE p (Expression x) = codePrec p x
codeE _ (Cos x) = showString "cos(" . codeE 0 x . showString ")"
codeE _ (Sin x) = showString "sin(" . codeE 0 x . showString ")"
codeE _ (Exp x) = showString "exp(" . codeE 0 x . showString ")"
codeE _ (Log x) = showString "log(" . codeE 0 x . showString ")"
codeE _ (Abs x) = showString "fabs(" . codeE 0 x . showString ")"
codeE _ (Signum _) = undefined
codeE _ (Product p i) | Product p i == 1 = showString "1"
codeE pree (Product p _) = showParen (pree > 7) $
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
codeE _ (Sum s i) | Sum s i == 0 = showString "0"
codeE p (Sum s _) = showParen (p > 6) (showString me)
  where me = foldl addup "" $ sum2pairs s
        addup "" (1,e) = codeE 6 e ""
        addup "" (f,e) = if e == 1
                         then show f
                         else show f ++ "*" ++ codeE 7 e ""
        addup rest (1,e) = codeE 6 e (showString " + " $ rest)
        addup rest (f,e) = show f ++ "*" ++ codeE 7 e (showString " + " $ rest)

latexE :: (Type a, Code a) => Int -> Expression a -> ShowS
latexE p (Var _ _ "" "" (Just e)) = latexE p e
latexE _ (Var _ _ c "" _) = showString c
latexE _ (Var _ _ _ t _) = showString t
latexE _ x | Just xx <- isConstant x = showString (latexDouble xx)
latexE p (Scalar x) = latexPrec p x
latexE p (Expression x) = latexPrec p x
latexE _ (Cos x) = showString "\\cos(" . latexE 0 x . showString ")"
latexE _ (Sin x) = showString "\\sin(" . latexE 0 x . showString ")"
latexE _ (Exp x) = showString "\\exp(" . latexE 0 x . showString ")"
latexE _ (Log x) = showString "\\log(" . latexE 0 x . showString ")"
latexE _ (Abs x) = showString "\\left|" . latexE 0 x . showString "\\right|"
latexE _ (Signum _) = undefined
latexE p (Product x _) | Map.size x == 1 && product2denominator x == 1 =
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
latexE pree (Product p _) | product2denominator p == 1 = latexParen (pree > 7) $ ltexsimple $ product2numerator p
  where ltexsimple [] = showString "1"
        ltexsimple [a] = latexE 7 a
        ltexsimple [a,b] = latexE 7 a . showString " " . latexE 7 b
        ltexsimple (a:es) = latexE 7 a . showString " " . ltexsimple es
latexE pree (Product p _) = latexParen (pree > 7) $ showString "\\frac{" . latexE 0 num . showString "}{" .
                                                                        latexE 0 den . showString "}"
  where num = product $ product2numerator p
        den = product2denominator p
latexE p (Sum s _) = latexParen (p > 6) (showString me)
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
toExpression 0 = Sum Map.empty Set.empty
toExpression 1 = Product Map.empty Set.empty
toExpression x = Sum (Map.singleton 1 (fromRational $ toRational x)) Set.empty

isConstant :: Type a => Expression a -> Maybe Double
isConstant (Sum s _) = case sum2pairs s of
                       [] -> Just 0
                       [(x,1)] -> Just x
                       _ -> Nothing
isConstant (Product p _) = if Map.size p == 0 then Just 1 else Nothing
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
  codeStatementHelper :: Expression a -> String -> Expression a -> String
  prefix :: Expression a -> String
  postfix :: Expression a -> String
  postfix _ = ""
  initialize :: Expression a -> String
  free :: Expression a -> String
  free _ = error "free nothing"
  toScalar :: a -> Expression Scalar
  fromScalar :: Expression Scalar -> Expression a
  fromScalar = Scalar




makeHomogeneous :: Type a => Expression a -> Expression Scalar
makeHomogeneous ee = 
  scalarScalar $ case isKSpace ee of
                    Same -> setZero (s_var "_kx" :: Expression Scalar) $ mapExpression toScalar ee
                    _ -> mapExpression toScalar ee
  where scalarScalar :: Expression Scalar -> Expression Scalar
        scalarScalar (Var t a b c (Just e)) = Var t a b c (Just $ scalarScalar e)
        scalarScalar (Var _ _ c l Nothing) = Var CannotBeFreed c c l Nothing
        scalarScalar (Scalar s) = s
        scalarScalar (Sum x _) = pairs2sum $ map f $ sum2pairs x
          where f (a,b) = (a, scalarScalar b)
        scalarScalar (Product x _) = pairs2product $ map f $ product2pairs x
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
  Sum a _ + Sum b _ = sumup a $ sum2pairs b
      where sumup x [] = map2sum x
            sumup x ((f,y):ys) = case Map.lookup y x of
                                 Nothing -> sumup (Map.insert y f x) ys
                                 Just f' -> if f + f' == 0
                                            then sumup (Map.delete y x) ys
                                            else sumup (Map.insert y (f+f') x) ys
  Sum a _ + b = case Map.lookup b a of
                Just fac -> if fac + 1 == 0
                            then case sum2pairs deleted of 
                                   [(1,e)] -> e
                                   [(f,e)] -> pairs2sum [(f,e)]
                                   _ -> map2sum deleted
                            else map2sum $ Map.insert b (fac + 1) a
                Nothing -> map2sum $ Map.insert b 1 a
    where deleted = Map.delete b a
  a + Sum i b = Sum i b + a
  a + b = Sum (Map.singleton a 1) Set.empty + b
  (-) = \x y -> x + (-y)
  -- the fromRational on the following line is needed to avoid a
  -- tricky infinite loop where -1 intepreted as (negate $ fromRational 1)
  negate = \x -> (fromRational $ -1) * x
  x * y | x == 0 = 0
        | y == 0 = 0
        | x == 1 = y
        | y == 1 = x
  Sum x i * y | Just c <- isConstant y = Sum (Map.map (\f -> c*f) x) i
  y * Sum x i | Just c <- isConstant y = Sum (Map.map (\f -> c*f) x) i
  Sum x _ * y | [(f,a)] <- sum2pairs x = pairs2sum [(f, a*y)]
  y * Sum x _ | [(f,a)] <- sum2pairs x = pairs2sum [(f, a*y)]
  Product a _ * Product b _ = puttogether a (product2pairs b)
      where puttogether x [] = map2product x
            puttogether x ((y,n):ys) =
                  case Map.lookup y x of
                    Just n' -> if n + n' == 0
                               then puttogether (Map.delete y x) ys
                               else puttogether (Map.insert y (n+n') x) ys
                    Nothing -> puttogether (Map.insert y n x) ys
  Product a _ * b = case Map.lookup b a of
                    Just n -> if n + 1 == 0
                             then map2product deleted
                             else map2product $ Map.insert b (n+1) a
                    Nothing -> map2product $ Map.insert b 1 a
    where deleted = Map.delete b a
  a * Product b i = Product b i * a
  a * b = Product (Map.singleton a 1) (varSet a) * b
  fromInteger = \x -> toExpression (fromInteger x :: Rational)
  abs = undefined
  signum = undefined

instance Type a => Fractional (Expression a) where
  x / y | Just yy <- isConstant y = x * toExpression (1/yy)
  x / Product y i = x * Product (Map.map negate y) i
  x / y = x * pairs2product [(y, -1)]
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
  (Sum x _) ** c | Just n <- isConstant c,
                   [(f,y)] <- sum2pairs x = pairs2sum [(f**n, y ** c)]
  (Product x _) ** c | Just n <- isConstant c = pairs2product $ map (p n) $ product2pairs x
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

countVars :: Type a => Expression a -> Int
countVars s = if newval == oldval then newval
                                  else error "bug in countVars!"
  where newval = Set.size $ varSet s
        oldval = length $ filter (/= "") $ nub $ varList s

varSet :: Type a => Expression a -> Set.Set String
varSet (Expression e) | Same <- isKSpace (Expression e), Kx <- e = Set.empty
                       | Same <- isKSpace (Expression e), Ky <- e = Set.empty
                       | Same <- isKSpace (Expression e), Kz <- e = Set.empty
                       | Same <- isRealSpace (Expression e), (IFFT e') <- e = varSet e'
                       | Same <- isKSpace (Expression e), (FFT e') <- e = varSet e'
                       | Same <- isScalar (Expression e), (Integrate e') <- e = varSet e'
                       | Same <- isKSpace (Expression e), SetKZeroValue _ e' <- e = varSet e'
                       | otherwise = error "There is no other possible expression"
varSet (Var _ _ _ _ (Just e)) = varSet e
varSet (Var IsTemp _ c _ Nothing) = Set.singleton c
varSet (Var CannotBeFreed _ _ _ Nothing) = Set.empty
varSet (Sin e) = varSet e
varSet (Cos e) = varSet e
varSet (Log e) = varSet e
varSet (Exp e) = varSet e
varSet (Abs e) = varSet e
varSet (Signum e) = varSet e
varSet (Sum _ i) = i
varSet (Product _ i) = i
varSet (Scalar _) = Set.empty

varList :: Type a => Expression a -> [String]
varList (Expression e) | Same <- isKSpace (Expression e), Kx <- e = []
                       | Same <- isKSpace (Expression e), Ky <- e = []
                       | Same <- isKSpace (Expression e), Kz <- e = []
                       | Same <- isRealSpace (Expression e), (IFFT e') <- e = varList e'
                       | Same <- isKSpace (Expression e), (FFT e') <- e = varList e'
                       | Same <- isScalar (Expression e), (Integrate e') <- e = varList e'
                       | Same <- isKSpace (Expression e), SetKZeroValue _ e' <- e = varList e'
                       | otherwise = error "There is no other possible expression"
varList (Var _ _ _ _ (Just e)) = varList e
varList (Var t _ c _ Nothing) | t == IsTemp = [c]
                                 | otherwise = []
varList (Sin e) = varList e
varList (Cos e) = varList e
varList (Log e) = varList e
varList (Exp e) = varList e
varList (Abs e) = varList e
varList (Signum e) = varList e
varList (Sum e _) = concat $ map (varList . snd) (sum2pairs e)
varList (Product e _) = concat $ map (varList . fst) (product2pairs e)
varList (Scalar _) = []

hasFFT :: Type a => Expression a -> Bool
hasFFT (Expression e) 
    | Same <- isKSpace (Expression e), FFT _ <- e = True
    | Same <- isRealSpace (Expression e), IFFT _ <- e = True
hasFFT (Var _ _ _ _ (Just e)) = hasFFT e
hasFFT (Sum s _) = or $ map (hasFFT . snd) (sum2pairs s)
hasFFT (Product p _) = or $ map (hasFFT . fst) (product2pairs p)
hasFFT (Sin e) = hasFFT e
hasFFT (Cos e) = hasFFT e
hasFFT (Log e) = hasFFT e
hasFFT (Exp e) = hasFFT e
hasFFT (Abs e) = hasFFT e
hasFFT (Signum e) = hasFFT e
hasFFT _ = False

hasK :: Type a => Expression a -> Bool
hasK (Expression e) 
    | Same <- isKSpace (Expression e), Kx <- e = True
    | Same <- isKSpace (Expression e), Ky <- e = True
    | Same <- isKSpace (Expression e), Kz <- e = True
    | otherwise = False
hasK (Var _ _ _ _ (Just e)) = hasK e
hasK (Sum s _) = or $ map (hasK . snd) (sum2pairs s)
hasK (Product p _) = or $ map (hasK . fst) (product2pairs p)
hasK (Sin e) = hasK e
hasK (Cos e) = hasK e
hasK (Log e) = hasK e
hasK (Exp e) = hasK e
hasK (Abs e) = hasK e
hasK (Signum e) = hasK e
hasK (Var _ _ _ _ Nothing) = False
hasK (Scalar _) = False

isfft :: Expression KSpace -> Maybe (Expression RealSpace)
isfft (Expression (FFT e)) = Just e
isfft (Product p _) = tofft 1 Nothing $ product2pairs p
  where tofft _ Nothing [] = Nothing
        tofft sc (Just e) [] = Just $ sc*e
        tofft sc Nothing ((Expression (FFT e),1):xs) = tofft sc (Just e) xs
        tofft sc myfft ((Scalar s,n):xs) = tofft (sc * Scalar s ** toExpression n) myfft xs
        tofft _ _ _ = Nothing
isfft (Var _ _ _ _ (Just e)) = isfft e
isfft _ = Nothing

isifft :: Expression RealSpace -> Maybe (Expression KSpace)
isifft (Expression (IFFT e)) = Just e
isifft (Product p _) = tofft 1 Nothing $ product2pairs p
  where tofft _ Nothing [] = Nothing
        tofft sc (Just e) [] = Just $ sc*e
        tofft sc Nothing ((Expression (IFFT e),1):xs) = tofft sc (Just e) xs
        tofft sc myfft ((Scalar s,n):xs) = tofft (sc * Scalar s ** toExpression n) myfft xs
        tofft _ _ _ = Nothing
isifft (Var _ _ _ _ (Just e)) = isifft e
isifft _ = Nothing

joinFFTs :: Type a => Expression a -> Expression a
--joinFFTs e | cleanvars e /= e = joinFFTs $ cleanvars e
joinFFTs (Expression e)
  | Same <- isKSpace (Expression e), FFT e' <- e = fft (joinFFTs e')
  | Same <- isRealSpace (Expression e), IFFT e' <- e = ifft (joinFFTs e')
joinFFTs (Sum s i) | Same <- isKSpace (Sum s i) = joinup [] $ sum2pairs s
  where joinup tofft [] = fft $ factorandsum tofft
        joinup tofft ((f,x):xs) | Just e <- isfft x = joinup (toExpression f*e:tofft) xs
                                | otherwise = toExpression f*x + joinup tofft xs
joinFFTs (Sum s i) | Same <- isRealSpace (Sum s i) = joinup [] $ sum2pairs s
  where joinup tofft [] = ifft $ factorandsum tofft
        joinup tofft ((f,x):xs) | Just e <- isifft x = joinup (toExpression f*e:tofft) xs
                                | otherwise = toExpression f*x + joinup tofft xs
joinFFTs (Var t a b c (Just e)) = Var t a b c (Just (joinFFTs e))
joinFFTs e = e

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
        countup f (Product y _) = Map.findWithDefault 0 f y
        countup f (Sum a _) | [(_,y)] <- sum2pairs a = countup f y
        countup f y | f == y = 1
        countup _ _ = 0
        getprodlist (Product xx _) = product2pairs xx
        getprodlist (Sum a _) | [(_,Product xx _)] <- sum2pairs a = product2pairs xx
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
derive v@(Var _ a b c _) dda (Var t aa bb cc (Just e))
  | Same <- compareTypes v dda =
    case isConstant $ derive v dda e of
      Just x -> toExpression x
      Nothing ->
        case isConstant $ derive v 1 e of
          Just x -> toExpression x * dda
          Nothing -> dda*(Var t ("d" ++ aa ++ "_by_d" ++ a) ("d" ++ bb ++ "_by_d" ++ b)
                          ("\\frac{\\partial "++cc ++"}{\\partial "++c++"}") $ Just $
                          (derive v 1 e))
derive v dda (Var _ _ _ _ (Just e)) = derive v dda e
derive _ _ (Var _ _ _ _ Nothing) = 0
derive v dda (Sum s _) = factorandsum $ map dbythis $ sum2pairs s
  where dbythis (f,x) = toExpression f * derive v dda x
derive v dda (Product p i) = factorandsum (map dbythis $ product2pairs p)
  where dbythis (x,n) = derive v (Product p i*toExpression n*dda/x) x
derive _ _ (Scalar _) = 0 -- FIXME
derive v dda (Cos e) = derive v (-dda*sin e) e
derive v dda (Sin e) = derive v (dda*cos e) e
derive v dda (Exp e) = derive v (dda*exp e) e
derive v dda (Log e) = derive v (dda/e) e
derive _ _ (Abs _) = error "I didn't think we'd need abs"
derive _ _ (Signum _) = error "I didn't think we'd need signum"
derive v dda (Expression e) = derivativeHelper v dda e

hasexpression :: (Type a, Type b) => Expression a -> Expression b -> Bool
hasexpression x e = countexpression x e > 0

countexpression :: (Type a, Type b) => Expression a -> Expression b -> Int
countexpression x e = snd $ subAndCount x (s_var "WeAreCounting") e

substitute :: (Type a, Type b) => Expression a -> Expression a -> Expression b -> Expression b
substitute x y e = fst $ subAndCount x y e

subAndCount :: (Type a, Type b) => Expression a -> Expression a -> Expression b -> (Expression b, Int)
subAndCount x y e | Same <- compareExpressions x e = (y, 1)
                  | not (Set.isSubsetOf (varSet x) (varSet e)) = (e, 0) -- quick check
subAndCount x@(Sum xs _) y e@(Sum es _)
  | Same <- compareTypes x e,
    ((factorX, termX):_) <- xspairs,
    Just factorE <- Map.lookup termX es,
    ratio <- factorE / factorX,
    Just es' <- filterout ratio xspairs es,
    (e'',n) <- subAndCount x y (mksum es') = (e'' + toExpression ratio*y, n+1)
  where xspairs = sum2pairs xs
        filterout ratio ((f,x1):rest) emap =
          case Map.lookup x1 emap of
            Just f' -> if ratio*f == f' then filterout ratio rest (Map.delete x1 emap)
                                        else Nothing
            Nothing -> Nothing
        filterout _ [] emap = Just emap
        mksum = pairs2sum . sum2pairs
subAndCount x@(Product xs _) y e@(Product es _)
  | Same <- compareTypes x e,
    ((termX, factorX):_) <- xspairs,
    Just factorE <- Map.lookup termX es,
    ratio <- factorE / factorX,
    Just es' <- filterout ratio xspairs es,
    (e'',n) <- subAndCount x y (mkproduct es') = (e'' * y**(toExpression ratio), n+1)
  where xspairs = product2pairs xs
        filterout ratio ((x1,f):rest) emap =
          case Map.lookup x1 emap of
            Just f' -> if ratio*f == f' then filterout ratio rest (Map.delete x1 emap)
                                        else Nothing
            Nothing -> Nothing
        filterout _ [] emap = Just emap
        mkproduct = pairs2product . product2pairs
subAndCount x y (Expression v)
  | Same <- isKSpace (Expression v), FFT e <- v = let (e', n) = subAndCount x y e in (Expression $ FFT e', n)
  | Same <- isRealSpace (Expression v), IFFT e <- v = let (e', n) = subAndCount x y e in (Expression $ IFFT e', n)
  | Same <- isScalar (Expression v), Integrate e <- v = let (e', n) = subAndCount x y e in (Expression $ Integrate e', n)
  | Same <- isKSpace (Expression v), v == Kx || v == Ky || v == Kz || v == Delta = (Expression v, 0)
  | Same <- isKSpace (Expression v), SetKZeroValue z e <- v = 
    let (e', n) = subAndCount x y e in (Expression $ SetKZeroValue z e', n)
  | otherwise = error $ "unhandled case in subAndCount: " ++ show v
subAndCount x y (Sum s i) = if n > 0 then (pairs2sum $ map justfe results, n)
                                     else (Sum s i, 0)
    where n = sum $ map (snd . snd) results
          results = map sub $ sum2pairs s
          justfe (f, (e, _)) = (f,e)
          sub (f, e) = (f, subAndCount x y e)
subAndCount x y (Product p i) = if num > 0 then (pairs2product $ map justen results, num)
                                           else (Product p i, 0)
    where num = sum $ map (snd . fst) results
          results = map sub $ product2pairs p
          justen ((e, _), n) = (e, n)
          sub (e, n) = (subAndCount x y e, n)
subAndCount x y (Cos e)   = (cos e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Sin e)   = (sin e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Log e)   = (log e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Exp e)   = (exp e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Abs e)   = (abs e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Signum e) = (signum e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Var t a b c (Just e)) = (Var t a b c (Just e'), n)
    where (e', n) = subAndCount x y e
subAndCount _ _ v@(Var _ _ _ _ Nothing) = (v, 0)
subAndCount x y (Scalar e) = (Scalar e', n)
    where (e', n) = subAndCount x y e



-- factorize :: Type a => Expression a -> Expression a
-- factorize (Expression e)
--   | Same <- isKSpace (Expression e), FFT e' <- e = fft (factorize e')
--   | Same <- isRealSpace (Expression e), IFFT e' <- e = ifft (factorize e')
-- factorize (Sum s) | Same <- isKSpace (Sum s) = joinup [] $ sum2pairs s
--   where joinup tofft [] = fft $ factorandsum tofft
--         joinup tofft ((f,x):xs) | Just e <- isfft x = joinup (toExpression f*e:tofft) xs
--                                 | otherwise = toExpression f*x + joinup tofft xs
-- factorize (Sum s) | Same <- isRealSpace (Sum s) = joinup [] $ sum2pairs s
--   where joinup tofft [] = ifft $ factorandsum tofft
--         joinup tofft ((f,x):xs) | Just e <- isifft x = joinup (toExpression f*e:tofft) xs
--                                 | otherwise = toExpression f*x + joinup tofft xs
-- factorize (Var a b c (Just e)) = Var a b c (Just (factorize e))
-- factorize e = e


\end{code}
