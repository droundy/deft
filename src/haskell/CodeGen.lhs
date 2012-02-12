\begin{code}
{-# LANGUAGE GADTs, PatternGuards #-}
module CodeGen ( RealSpace(..), r_var,
                 KSpace(..), k_var, kx, ky, kz, k, ksqr,
                 Scalar(..), s_var,
                 fft, ifft, integrate, grad, derive,
                 Expression,                 
                 Statement(..),
                 Type, 
                 makeHomogeneous, isConstant, hasexpression, -- for debugging only!!!
                 code, latex, setZero, codeStatement, substitute,
                 generateHeader, simp, simp2, countFFT, checkDup)
       where

--import Debug.Trace

import qualified Data.Map as Map
import Hash ( hash )
import Data.List ( nubBy, (\\) )

\end{code}

The \verb!RealSpaceField! data type describes a field in real space.

\begin{code}
data RealSpace = R String |
                 IFFT (Expression KSpace)
               deriving ( Eq, Ord, Show )
data KSpace = K String |
              Delta | -- handy for FFT of homogeneous systems
              Kx | Ky | Kz |
              FFT (Expression RealSpace)
            deriving ( Eq, Ord, Show )
data Scalar = S String |
              Integrate (Expression RealSpace)
            deriving ( Eq, Ord, Show )

instance Code RealSpace where
  codePrec _ (R v) = showString (v ++ "[i]")
  codePrec _ (IFFT ksp@(Expression (K _))) = showString "ifft(gd, " . codePrec 0 (makeHomogeneous ksp) . showString ")"
  codePrec _ (IFFT ke) = showString "ifft(gd, " . codePrec 0 ke . showString ")"
  latexPrec _ (R v@('{':_)) = showString v
  latexPrec _ (R (a:v@(_:_))) = showString (a : '_' : '{' : v ++ "}")
  latexPrec _ (R v) = showString v
  latexPrec _ (IFFT ke) = showString "\\text{ifft}\\left(" . latexPrec 0 ke . showString "\\right)"
instance Type RealSpace where
  isRealSpace _ = Same
  derivativeHelper v ddr r | Same <- isRealSpace (Expression v), v == r = ddr
  derivativeHelper v ddr (IFFT ke) = derive v (fft ddr) ke
  derivativeHelper _ _ _ = 0
  zeroHelper v x | Same <- isRealSpace (Expression v), v == x = 0
  zeroHelper v (IFFT ke) = ifft (setZero v ke)
  zeroHelper _ x = Expression x
  simpHelper (R x) = ([], r_var x)
  simpHelper (IFFT ksp@(Expression (K _))) = (st++[InitializeR (r_var temp), AssignR temp (ifft ksp')], r_var temp) 
      where (st, ksp') = simp ksp
            Expression (R temp) = r_var ("tempIFFT" ++ "_" ++ hash (show (ifft ksp')))
  simpHelper (IFFT ksp) = (st++[InitializeK (k_var tempk), AssignK tempk ksp', InitializeR (r_var temp), AssignR temp (ifft ksp'')], r_var temp)
      where (st, ksp') = simp ksp
            ksp''@(Expression (K tempk)) = k_var ("temp" ++ "_" ++ hash (show ksp'))
            Expression (R temp) = r_var ("tempIFFT" ++ "_" ++ hash (show (ifft ksp'')))
  var v (Scalar _) = s_var v
  var v _ = r_var v
  prefix "" (Expression (IFFT _)) = ""
  prefix "" _ = "for (int i=0; i<gd.NxNyNz; i++) {\n\t\t"
  prefix v _ = "Grid " ++ v ++ "(gd);\nfor (int i=0; i<gd.NxNyNz; i++) {\n\t\t"
  postfix (Expression (IFFT _)) = ""
  postfix _ = "\t\n}\n"
  codeStatementHelper a op (Expression (IFFT (Expression (K v)))) = a ++ op ++ "ifft(gd, " ++ v ++ ");\n"
  codeStatementHelper _ _ (Expression (IFFT e)) = error ("It is a bug to generate code for a non-var input to ifft\n"++ latex e)
  codeStatementHelper a op e = prefix "" e ++ a ++ "[i]" ++ op ++ code e ++ ";\n" ++ postfix e
  initialize (Expression (R x)) = "VectorXd " ++ x ++ "(gd.NxNyNz);\n\t//printf(\"Memory use alloc " ++ x ++ " is %g with peak %g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n"
  initialize _ = "VectorXd output(gd.NxNyNz);\n\t//printf(\"Memory use output is %g with peak %g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n"
  free (Expression (R x)) = x ++ ".resize(0);\n\t//printf(\"Memory use free " ++ x ++ " is %g with peak %g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n"
  free _ = error "free error"
  toScalar (R v) = s_var v
  toScalar (IFFT ke) = makeHomogeneous ke

instance Code KSpace where
  codePrec _ (K v) = showString (v ++ "[i]")
  codePrec _ Kx = showString "k_i[0]"
  codePrec _ Ky = showString "k_i[1]"
  codePrec _ Kz = showString "k_i[2]"
  codePrec _ Delta = showString "delta(k?)"
  codePrec _ (FFT r) = showString "fft(gd, " . codePrec 0 (makeHomogeneous r) . showString ")"
  latexPrec _ (K v@('{':_)) = showString v
  latexPrec _ (K (a:v@(_:_))) = showString (a : '_' : '{' : v ++ "}")
  latexPrec _ (K v) = showString v
  latexPrec _ Kx = showString "kx"
  latexPrec _ Ky = showString "ky"
  latexPrec _ Kz = showString "kz"
  latexPrec _ Delta = showString "delta(k)"
  latexPrec _ (FFT r) = showString "\\text{fft}(" . latexPrec 0 r . showString ")"
instance Type KSpace where
  isKSpace _ = Same
  derivativeHelper v ddk kk | Same <- isKSpace (Expression v), kk == v = ddk
  derivativeHelper v ddk (FFT r) = derive v (ifft ddk) r
  derivativeHelper _ _ _ = 0
  zeroHelper v x | Same <- isKSpace (Expression v), x == v = 0
  zeroHelper v (FFT r) = fft (setZero v r)
  zeroHelper _ x = Expression x
  simpHelper (FFT rsp@(Expression (R _))) = (st ++ [InitializeK (k_var temp), AssignK temp (fft rsp')], k_var temp)
      where (st, rsp') = simp rsp
            Expression (K temp) = k_var ("tempFFT" ++ "_" ++ hash (show (fft rsp')))
  simpHelper (FFT rsp) = (st++[InitializeR (r_var tempr), AssignR tempr rsp', InitializeK (k_var temp), AssignK temp (fft rsp'')], k_var temp)
      where (st, rsp') = simp rsp
            rsp''@(Expression (R tempr)) = r_var ("tempR" ++ "_" ++ hash (show rsp'))
            Expression (K temp) = k_var ("tempFFT" ++ "_" ++ hash (show (fft rsp'')))
  simpHelper ksp = ([], Expression ksp)
  var v (Scalar _) = s_var v
  var v _ = k_var v
  prefix "" (Expression (FFT _)) = ""
  prefix "" _ = "for (int i=1; i<gd.NxNyNzOver2; i++) {\n\t\tconst int z = i % gd.NzOver2;\n\t\tconst int n = (i-z)/gd.NzOver2;\n\t\tconst int y = n % gd.Ny;\n\t\tconst int xa = (n-y)/gd.Ny;\n\t\tconst RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);\n\t\tconst Reciprocal k_i = gd.Lat.toReciprocal(rvec);\n\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3); assert(dr);\n"
  prefix v _ = "ReciprocalGrid " ++ v ++ "(gd);\nfor (int i=0; i<gd.NxNyNz; i++) {\n\t\t"
  postfix (Expression (FFT _)) = ""
  postfix _ = "\t}\n"
  codeStatementHelper a op (Expression (FFT (Expression (R v)))) = a ++ op ++ "fft(gd, " ++ v ++ ");\n"
  codeStatementHelper _ _ (Expression (FFT _)) = error "It is a bug to generate code for a non-var input to fft"
  codeStatementHelper a op e = "{\n\t\tconst int i = 0;\n\t\tconst int z = i % gd.NzOver2;\n\t\tconst int n = (i-z)/gd.NzOver2;\n\t\tconst int y = n % gd.Ny;\n\t\tconst int xa = (n-y)/gd.Ny;\n\t\tconst RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);\n\t\tconst Reciprocal k_i = gd.Lat.toReciprocal(rvec);\n\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3); assert(dr);\n\t\t" ++ a ++ "[0]" ++ op ++ code (setZero Kz (setZero Ky (setZero Kx e))) ++ ";\n\t}\n\t" ++ prefix "" e ++ "\t\t" ++ a ++ "[i]"  ++ op ++ code e ++ ";" ++ postfix e
  initialize (Expression (K x)) = "VectorXcd " ++ x ++ "(gd.NxNyNzOver2);\n\t//printf(\"Memory use k alloc " ++ x ++ " is %g with peak %g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n"
  initialize _ = "VectorXcd output(gd.NxNyNzOver2);\n\t//printf(\"Memory use k output is %g with peak %g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n"
  free (Expression (K x)) = x ++ ".resize(0);\n\t//printf(\"Memory use free " ++ x ++ " is %g with peak %g\\n\", current_memory()/1024.0/1024, peak_memory()/1024.0/1024);\n"
  free _ = error "free error"
  toScalar (K v) = s_var v
  toScalar Delta = 1
  toScalar Kx = s_var "_kx"
  toScalar Ky = 0
  toScalar Kz = 0
  toScalar (FFT e) = makeHomogeneous e

mapExpression :: (Type a, Type b) => (a -> Expression b) -> Expression a -> Expression b
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

isEven :: (Type a, Type b) => b -> Expression a -> Double
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
isEven v (Expression e) | Same <- isKSpace (Expression v), 
                          Same <- isKSpace (Expression e),
                          v == e = -1
isEven v (Expression e) | Same <- isRealSpace (Expression v), 
                          Same <- isRealSpace (Expression e),
                          v == e = -1
isEven v (Expression e) | Same <- isScalar (Expression v), 
                          Same <- isScalar (Expression e),
                          v == e = -1
isEven v (Expression e) 
  | Same <- isKSpace (Expression e) =
    case e of
      FFT r -> isEven v r
      _ -> 1
isEven v (Expression e) 
  | Same <- isRealSpace (Expression e) =
    case e of
      IFFT ks -> isEven v ks
      _ -> 1
isEven v (Expression e) 
  | Same <- isScalar (Expression e) =
    case e of
      Integrate x -> isEven v x
      _ -> 1
isEven _ (Expression _) = 1 -- Technically, it might be good to recurse into this

setZero :: (Type a, Type b) => b -> Expression a -> Expression a
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
  if --trace ("isEven " ++ latex (Product p) ++ " = " ++ show (isEven v (Product p)))
     isEven v (Product p) == -1
  then 0
  else 
    if --trace ("setZero Product " ++ latex (Product p)) 
       zd /= 0
    then zn / zd
    else if --trace ("L'Hopital's rule: " ++ latex (Product p)) 
            zn /= 0
         then error ("L'Hopital's rule failure: " ++ latex n ++ "\n /\n  " ++ latex d ++ "\n\n\n" 
                     ++ latex (Product p) ++ "\n\n\n" ++ latex zn)
         else case isKSpace (Expression v) of
              Same -> 
                case isKSpace n of
                  Same -> --trace ("Need to derive: dtop is " ++ latex dtop ++ " dbot is " ++ latex dbot) 
                          setZero v (dtop / dbot)
                    where dtop = derive v 1 n
                          dbot = derive v 1 d
                  _ -> error "oopsies"
              _ -> 
                case isScalar (Expression v) of
                  Same -> 
                    case isScalar n of
                      Same -> setZero v (derive v 1 n / derive v 1 d)
                      _ -> error "oopsies"
                  _ -> 
                    case isRealSpace (Expression v) of
                      Same -> 
                        case isRealSpace n  of
                          Same -> --trace ("Need to derive: dtop is " ++ latex dtop ++ " dbot is " ++ latex dbot
                                  --       ++
                                  --       "\n\ttop is " ++ latex n ++ " bot is " ++ latex d) 
                                  setZero v (dtop / dbot)
                            where dtop = derive v 1 n
                                  dbot = derive v 1 d
                          _ -> error "oopsies"
                      _ -> error "oops" -- setZero v (derive v 1 a / derive v 1 b)
  where d = product2denominator p
        n = product $ product2numerator p
        zn = setZero v n
        zd = setZero v d
setZero _ (Sum s) | Sum s == 0 = 0
setZero v (Sum s) = --trace ("out " ++ latex (Sum s) ++ " = " ++ show (map show $ map sz $ sum2pairs s) 
                    --      ++ " = " ++ latex out)
                    out
  where sz (f,x) = (f, setZero v x)
        out = pairs2sum $ map sz $ sum2pairs s
setZero v (Expression x) = zeroHelper v x

instance Code Scalar where
  codePrec _ (S v) = showString v
  codePrec _ (Integrate r) = showString "(" . codePrec 0 r . showString ") * gd.dvolume"
  latexPrec _ (S v@('{':_)) = showString v
  latexPrec _ (S (a:v@(_:_))) = showString (a : '_' : '{' : v ++ "}")
  latexPrec _ (S v) = showString v
  latexPrec _ (Integrate r) = showString "integrate(" . latexPrec 0 r . showString ")"
instance Type Scalar where
  s_var = Expression . S
  isScalar _ = Same
  derivativeHelper v dds (Integrate e) = derive v (Scalar dds*s_var "dV") e
  derivativeHelper v dds s | Same <- isScalar (Expression v), v == s = dds
  derivativeHelper _ _ _ = 0
  zeroHelper v x | Same <- isScalar (Expression v), v == x = 0
  zeroHelper v (Integrate e) = integrate (setZero v e)
  zeroHelper _ x = Expression x
  simpHelper (Integrate r) = (st, integrate r')
      where (st, r') = simp r
  simpHelper sc = ([], Expression sc)
  var v _ = s_var v
  codeStatementHelper a _ (Expression (Integrate e)) = "for (int i=0; i<gd.NxNyNz; i++) {\n\t\t" ++ a ++ " += " ++ code e ++ ";\n\t}\n"
  codeStatementHelper a op e = a ++ op ++ code e ++ ";\n\t"
  prefix "" (Expression (Integrate _))  = "for (int i=0; i<gd.NxNyNz; i++) {    "
  prefix _ _ = ""
  postfix (Expression (Integrate _)) = "\n}\n"
  postfix _ = ""
  initialize (Expression (Integrate _)) = "double output = 0;\n"
  initialize (Expression (S x)) = "double " ++ x ++ ";\n"
  initialize _ = "double output = 0;\n"
  toScalar (Integrate r) = makeHomogeneous r
  toScalar v = Expression v

r_var :: String -> Expression RealSpace
r_var v = Expression (R v)

k_var :: String -> Expression KSpace
k_var v = Expression (K v)

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

integrate :: Expression RealSpace -> Expression Scalar
integrate x = Expression (Integrate x)

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
        helper [(x,Sum y)] | [(x2,y')] <- sum2pairs y = helper [(x*x2, y')]
        helper es = Sum $ fl (Map.empty) es
        fl a [] = a
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
                         else show f ++ "*" ++ codeE 6 e ""
        addup rest (1,e) = codeE 6 e (showString " + " $ rest)
        addup rest (f,e) = show f ++ "*" ++ codeE 6 e (showString " + " $ rest)

latexE :: (Type a, Code a) => Int -> Expression a -> ShowS
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
                                            then latexE 8 e . showString ("^\\frac" ++ show n2 ++ "2")
                                            else latexE 8 e . showString ("^\\frac{" ++ show n2 ++ "}2")
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
latexE pree (Product p) = latexParen (pree > 7) $ showString "\\frac{" . latexE 7 num . showString "}{" .
                                                                        latexE 7 den . showString "}"
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
        addup ('-':rest) (f,e) = latexDouble f ++ " " ++ latexE 6 e (showString " - " $ rest)
        addup rest (1,e) = latexE 6 e (" + " ++ rest)
        addup rest (f,e) = latexDouble f ++ " " ++ latexE 6 e (showString " + " $ rest)

latexParen :: Bool -> ShowS -> ShowS
latexParen False x = x
latexParen True x = showString "\\left(" . x . showString "\\right)"

latexDouble :: Double -> String
latexDouble 0 = "0"
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
  derivativeHelper :: Type b => b -> Expression a -> a -> Expression b
  zeroHelper :: Type b => b -> a -> Expression a
  simpHelper :: a -> ([Statement], Expression a)
  var :: String -> Expression a -> Expression a
  codeStatementHelper :: String -> String -> Expression a -> String
  prefix :: String -> Expression a -> String
  postfix :: Expression a -> String
  postfix _ = ""
  initialize :: Expression a -> String
  initialize _ = ""
  free :: Expression a -> String
  free _ = error "free nothing"
  toScalar :: a -> Expression Scalar

codeStatements :: [Statement] -> String
codeStatements = unlines . map codeStatement

codeStatement :: Statement -> String
codeStatement (AssignR x y) = "\t" ++ codeStatementHelper x " = " y
codeStatement (AssignK x y) = "\t" ++ codeStatementHelper x " = " y
codeStatement (AssignS x y) = "\t" ++ codeStatementHelper x " = " y
codeStatement (InitializeR e) = "\t" ++ initialize e
codeStatement (InitializeK e) = "\t" ++ initialize e
codeStatement (InitializeS e) = "\t" ++ initialize e
codeStatement (FreeR e) = "\t" ++ free e
codeStatement (FreeK e) = "\t" ++ free e

{-
codeStatement :: Statement a -> String
codeStatement (a := e) = "\t" ++ codeStatementHelper a " = " e
codeStatement (a :?= e) = "\t" ++ codeStatementHelper (a ++ "_" ++ hash (show e)) " = " e
codeStatement (a :+= e) = "\t" ++ codeStatementHelper a " += " e
codeStatement (Initialize e) = "\t" ++ initialize e
codeStatement (Return _) = ""
codeStatement (a :>> b) = codeStatement a ++ codeStatement b
-}

countFFT :: [Statement] -> Int
countFFT = sum . map helper
    where helper (AssignR _ (Expression (IFFT _))) = 1
          helper (AssignK _ (Expression (FFT _))) = 1
          helper _ = 0

makeHomogeneous :: Type a => Expression a -> Expression Scalar
makeHomogeneous ee = 
  scalarScalar $ case isKSpace ee of
                    Same -> setZero (S "_kx") $ mapExpression toScalar ee
                    _ -> mapExpression toScalar ee
  where scalarScalar :: Expression Scalar -> Expression Scalar
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
grad v e = derive (R v) 1 e

data ToDo = DoK (Expression KSpace) | DoR (Expression RealSpace) | DoNothing
            deriving Eq

findToDo :: Type a => Expression a -> ToDo
findToDo (Expression e)
    | Same <- isKSpace (Expression e), FFT (Expression (R _)) <- e = DoK (Expression e)
    | Same <- isKSpace (Expression e), FFT e' <- e = case findToDo e' of
                                                       DoNothing -> DoR $ e'
                                                       dothis -> dothis
    | Same <- isRealSpace (Expression e), IFFT (Expression (K _)) <- e = DoR (Expression e)
    | Same <- isRealSpace (Expression e), IFFT e' <- e = case findToDo e' of
                                                           DoNothing -> DoK $ e'
                                                           dothis -> dothis
    | Same <- isScalar (Expression e), Integrate e' <- e = findToDo e'
findToDo (Sum s) = case filter (/= DoNothing) $ map sub $ sum2pairs s of
                     [] -> DoNothing
                     dothis:_ -> dothis
    where sub (_,e) = findToDo e
findToDo (Product p) = case filter (/= DoNothing) $ map sub $ product2pairs p of
                         [] -> DoNothing
                         dothis:_ -> dothis
    where sub (e, _) = findToDo e
findToDo (Cos e) = findToDo e
findToDo (Sin e) = findToDo e
findToDo (Log e) = findToDo e
findToDo (Exp e) = findToDo e
findToDo (Abs e) = findToDo e
findToDo (Signum e) = findToDo e
findToDo _ = DoNothing

simp2 :: Type a => Expression a -> ([Statement], Expression a)
simp2 = simp2helper (0 :: Int) []
    where simp2helper n sts e = case findToDo e of
                                  DoK ke -> simp2helper (n+1) (sts++[inike, setke]) e'
                                      where v = "var" ++ show n
                                            inike = InitializeK $ k_var v
                                            setke = AssignK v ke
                                            e'  = substitute ke (k_var v) e
                                  DoR re -> simp2helper (n+1) (sts++[inire, setre]) e'
                                      where v = "var" ++ show n
                                            inire = InitializeR $ r_var v
                                            setre = AssignR v re
                                            e'  = substitute re (r_var v) e
                                  DoNothing -> (sts, e)

simp :: Type a => Expression a -> ([Statement], Expression a)
simp (Expression e) = simpHelper e
simp (Sum s) = (concatMap fst simped, pairs2sum $ map snd simped)
  where es = sum2pairs s
        simped = map simpme es
        simpme (f,x) = (st, (f,x'))
          where (st, x') = simp x
simp (Cos e) = (st, Cos e')
    where (st, e') = simp e
simp (Sin e) = (st, Sin e')
    where (st, e') = simp e
simp (Log e) = (st, Log e')
    where (st, e') = simp e
simp (Exp e) = (st, Exp e')
    where (st, e') = simp e
simp (Scalar s) = (st, Scalar s')
    where (st, s') = simp s
simp (Product p) = (concatMap fst simped, product $ map snd simped)
  where es = product2pairs p
        simped = map simpme es
        simpme (x,n) = (st, x' ** toExpression n)
          where (st,x') = simp x
simp _ = error "simp incomplete"

factorandsum :: Type a => [Expression a] -> Expression a
factorandsum [] = 0
factorandsum (x:xs) = helper (getprodlist x) (x:xs)
  where helper :: Type a => [(Expression a, Double)] -> [Expression a] -> Expression a
        helper [] as = sum as
        helper ((f,n):fs) as = if --trace ("factoring " ++ show n' ++ " copies of "++latex f) 
                                  n' /= 0
                               then (helper fs (map (/(f**(toExpression n'))) as)) * (f**(toExpression n'))
                               else helper fs as
          where n' = if n < 0 then if --trace (show f ++ show (map (countup f) as)) 
                                      maximum (map (countup f) as) < 0
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

derive :: (Type a, Type b) => b -> Expression a -> Expression a -> Expression b
derive v dda (Sum s) = sum $ map dbythis $ sum2pairs s
  where dbythis (f,x) = toExpression f * derive v dda x
derive v dda (Product p) = factorandsum (map dbythis $ product2pairs p)
  where dbythis (x,n) = derive v (Product p*toExpression n*dda/x) x
derive v _ (Scalar x) = derive v 1 x -- FIXME
derive v dda (Cos e) = derive v (-dda*sin e) e
derive v dda (Sin e) = derive v (dda*cos e) e
derive v dda (Exp e) = derive v (dda*exp e) e
derive v dda (Log e) = derive v (dda/e) e
derive _ _ (Abs _) = error "I didn't think we'd need abs"
derive _ _ (Signum _) = error "I didn't think we'd need signum"
derive v dda (Expression e) = derivativeHelper v dda e

hasexpression :: (Type a, Type b) => Expression a -> Expression b -> Bool
hasexpression x e | Same <- isKSpace e, Same <- isKSpace x, x == e = True
                  | Same <- isRealSpace e, Same <- isRealSpace x, x == e = True
                  | Same <- isScalar e, Same <- isScalar x, x == e = True
hasexpression x (Expression v) | Same <- isKSpace (Expression v), FFT e <- v = hasexpression x e
                               | Same <- isRealSpace (Expression v), IFFT e <- v = hasexpression x e
                               | Same <- isScalar (Expression v), Integrate e <- v = hasexpression x e
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
hasexpression _ _ = False

substitute :: (Type a, Type b) => Expression a -> Expression a -> Expression b -> Expression b
substitute x y e | Same <- isKSpace e, Same <- isKSpace x, x == e = y
                 | Same <- isRealSpace e, Same <- isRealSpace x, x == e = y
                 | Same <- isScalar e, Same <- isScalar x, x == e = y
substitute x y (Expression v) | Same <- isKSpace (Expression v), FFT e <- v = Expression $ FFT (substitute x y e)
                              | Same <- isRealSpace (Expression v), IFFT e <- v = Expression $ IFFT (substitute x y e)
                              | Same <- isScalar (Expression v), Integrate e <- v = Expression $ Integrate (substitute x y e)
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
substitute _ _ e = e

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
codeS _ (AssignR x y) = showString (prefix "" y) . codePrec 0 (var x y) . showString " = " . codePrec 0 y . showString ";" . showString (postfix y)
codeS _ (AssignK x y) = showString (prefix "" y) . codePrec 0 (var x y) . showString " = " . codePrec 0 y . showString ";" . showString (postfix y)
codeS _ (AssignS x y) = showString (prefix "" y) . codePrec 0 (var x y) . showString " = " . codePrec 0 y . showString ";" . showString (postfix y)
codeS _ (InitializeR e) = showString (initialize e)
codeS _ (InitializeK e) = showString (initialize e)
codeS _ (InitializeS e) = showString (initialize e)
codeS _ (FreeR e) = showString (free e)
codeS _ (FreeK e) = showString (free e)

latexS :: Int -> Statement -> ShowS
latexS _ (AssignR x y) = latexPrec 0 (var x y) . showString " = " . latexPrec 0 y
latexS _ (AssignK x y) = latexPrec 0 (var x y) . showString " = " . latexPrec 0 y
latexS _ (AssignS x y) = latexPrec 0 (var x y) . showString " = " . latexPrec 0 y
latexS _ (InitializeR e) = showString (initialize e)
latexS _ (InitializeK e) = showString (initialize e)
latexS _ (InitializeS e) = showString (initialize e)
latexS _ (FreeR e) = showString (free e)
latexS _ (FreeK e) = showString (free e)

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
                    ns = map (FreeR . r_var) fe
                    check n = (hasexpression (r_var n) xn) || (hasexpression (k_var n) xn)
          freeHelper v (x@(AssignK _ xn):xs) | or $ map check v = ns ++ [x] ++ freeHelper (v \\ fe) xs
              where fe = filter check v
                    ns = map (FreeK . k_var) fe
                    check n = (hasexpression (r_var n) xn) || (hasexpression (k_var n) xn)
          freeHelper v (x:xs) = [x] ++ freeHelper v xs
          freeHelper _ [] = []

vecInMem :: [Statement] -> [String]
vecInMem s = filter isFreeVec $ map ini s
     where ini (InitializeR (Expression (R x))) = x
           ini (InitializeK (Expression (K x))) = x
           ini _ = ""
           fre (FreeR (Expression (R x))) = x
           fre (FreeK (Expression (K x))) = x
           fre _ = ""
           isFreeVec x = not $ elem x $ map fre s

\end{code}


\begin{code}
{-
data Statement a where
     (:=) :: (Type a, Code a) => String -> Expression a -> Statement (Expression a)
     (:?=) :: (Type a, Code a) => String -> Expression a -> Statement (Expression a)
     (:+=) :: (Type a, Code a) => String -> Expression a-> Statement (Expression a)
     Initialize :: (Type a) => Expression a -> Statement (Expression a)
     Return :: a -> Statement a
     (:>>) :: Statement b -> Statement a -> Statement a
infixl 1 :>>
infix 4 :=, :?=, :+=
instance Monad Statement where
  Return _ >> x = x
  (a :>> Return _) >> b = a >> b
  a >> b = a :>> b
  (>>=) = thenS
  return = Return
instance Show (Statement a) where
    showsPrec = showsS
showsS :: Int -> Statement a -> ShowS
showsS _ (x := y) = showString x . showString " := " . showsPrec 0 y
showsS _ (x :?= y) = showString (x ++ "_" ++ hash (show y)) . showString " :?= " . showsPrec 0 y
showsS _ (x :+= y) = showString x . showString " :+= " . showsPrec 0 y
showsS _ (Initialize _) = showString "initialize ???"
showsS _ (Return _) = showString "return ???" 
showsS p (x :>> y) = showsS p x . showString "\n" . showsS p y

instance Code (Statement a) where
  codePrec = codeS
  latexPrec = latexS
codeS :: Int -> Statement a -> ShowS
codeS _ (x := y) = showString (prefix "" y) . 
                    codePrec 0 (var x y) . showString " = " . codePrec 0 y . showString ";" .
                    showString (postfix y)
codeS _ (x :?= y) = showString (prefix x y) .
                     codePrec 0 (var x y) . showString " := " . codePrec 0 y . showString ";" .
                     showString (postfix y)
codeS _ (x :+= y) = showString (prefix "" y) .
                     codePrec 0 (var x y) . showString " += " . codePrec 0 y . showString ";" .
                     showString (postfix y)
codeS _ (Initialize e) = showString (initialize e)
codeS _ (Return _) = showString ""
codeS p (x :>> y) = codeS p x . showString "\n" . codeS p y

latexS :: Int -> Statement a -> ShowS
latexS _ (x := y) = latexPrec 0 (var x y) . showString " = " . latexPrec 0 y
latexS _ (x :?= y) = latexPrec 0 (var x y) . showString " := " . latexPrec 0 y
latexS _ (x :+= y) =  latexPrec 0 (var x y) . showString " += " . latexPrec 0 y
latexS _ (Initialize e) = showString (initialize e)
latexS _ (Return _) = showString ""
latexS p (x :>> y) = latexS p x . showString "\n" . latexS p y


thenS :: Statement a -> (a -> Statement b) -> Statement b
thenS f g = f :>> g (evalS f)

evalS :: Statement a -> a
evalS (s := e) = var s e
evalS (s :?= e) = var (s ++ "_" ++ hash (show e)) e
evalS (s :+= e) = var s e
evalS (Initialize a) = a
evalS (Return a) = a
evalS (_ :>> b) = evalS b
{-
combineS :: Statement a -> Statement b -> (a -> b -> a) -> Statement a
combineS (s1 := e1) (_ := e2) op = (s1 := (op e1 e2))
combineS (s1 :?= e1) (_ := e2) op = (s1 :?= (op e1 e2))
combineS s1 (s2 :?= e) op = combineS s1 (s2 := e) op
combineS s1 (s2 :>> s3) op = s2 :>> (combineS s1 s3 op)
combineS s1 _ _ = s1
-}
-}
\end{code}

\begin{code}

functionCode :: String -> String -> [(String, String)] -> String -> String
functionCode "" "" [] "" = ""
functionCode "" "" (x:xs) "" = if xs == [] 
                               then fst x ++ " " ++ snd x
                               else fst x ++ " " ++ snd x ++ ", " ++ functionCode "" "" xs ""
functionCode n t a b = t ++ " " ++ n ++ "(" ++ functionCode "" "" a "" ++ ") const {\n" ++ b ++ "\n}\n"



classCode :: Expression RealSpace -> Maybe (Expression RealSpace) -> String -> String
classCode e arg n = "class " ++ n ++ " : public FunctionalInterface {\npublic:\n" ++ n ++ codeA arg ++ "  {\n\thave_integral = true;\n}\n" ++
                functionCode "I_have_analytic_grad" "bool" [] "\treturn false;" ++
                functionCode "integral" "double" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n\tdouble output=0;\n" ++ codeStatements codeIntegrate ++ "\t// " ++ show (countFFT codeIntegrate) ++ " Fourier transform used.\n" ++ "\treturn output*gd.dvolume;\n" ) ++
                functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n" ++ codeStatements codeVTransform  ++ "\t// " ++ show (countFFT codeVTransform) ++ " Fourier transform used.\n" ++ "\treturn output;\n")  ++
                functionCode "transform" "double" [("double", "kT"), ("double", "x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(x==x); // to avoid an unused parameter error\n\t" ++ codeStatements codeDTransform ++ "\n\treturn output;\n") ++
                functionCode "append_to_name" "bool" [("", "std::string")] "\treturn false;" ++
                functionCode "derive" "double" [("double", "kT"), ("double", "x")] ("\tassert(kT==kT);\n\tassert(x==x);\n" ++ codeStatements codeDerive ++ "\n\treturn output;\n") ++
                functionCode "d_by_dT" "double" [("double", ""), ("double", "")] "\tassert(0); // fail\n\treturn 0;\n" ++
                functionCode "derive_homogeneous" "Expression" [("const Expression &", "")] "\tassert(0); // fail\n\treturn Expression(0);\n" ++
                functionCode "grad" "Functional" [("const Functional", "&ingrad"), ("const Functional", "&x"), ("bool", "")] "\tassert(&ingrad==&ingrad);\n\tassert(&x==&x);\n\treturn ingrad;" ++
                functionCode "grad_T" "Functional" [("const Functional", "&ingradT")] "\tassert(&ingradT==&ingradT);\n\treturn ingradT;" ++
                functionCode "grad" "void" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x"), ("const VectorXd", "&ingrad"), ("VectorXd", "*outgrad"), ("VectorXd", "*outpgrad")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n\tassert(outpgrad==outpgrad);\n" ++ codeStatements codeGrad ++ "\t// " ++ show (countFFT codeGrad) ++ " Fourier transform used.\n") ++
                functionCode "printme" "Expression" [("const Expression", "&x")] ("\treturn funexpr(\"" ++ n ++ "()\")(x);") ++
                functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")] "\tFunctionalInterface::print_summary(prefix, energy, name);" ++
                "private:\n"++ codeArgInit arg  ++"}; // End of " ++ n ++ " class"
    where
      codeIntegrate = freeVectors (st ++ [AssignS "output" e'])
          where (st, e') = simp2 (integrate e)
      codeVTransform = freeVectors (InitializeR e: st ++ [AssignR "output" e'])
          where (st, e') = simp2 e
      codeDTransform = freeVectors [InitializeS (makeHomogeneous e), AssignS "output" (makeHomogeneous e)]
      codeDerive = freeVectors [InitializeS (makeHomogeneous (derive (R "x") 1 e)), AssignS "output" (makeHomogeneous (derive (R "x") 1 e))]
      codeGrad = freeVectors (st ++ [AssignR "(*outgrad)" (r_var "(*outgrad)" + e')])
          where (st, e') = simp2 (derive (R "x") (r_var "ingrad") e )
      codeA (Just rep) = "(double " ++ code (makeHomogeneous rep) ++ "_arg) : " ++ code (makeHomogeneous rep) ++ "(" ++ code (makeHomogeneous rep) ++ "_arg)"
      codeA Nothing = "()"
      codeArgInit (Just rep) = code (InitializeS (makeHomogeneous rep))
      codeArgInit Nothing = ""

generateHeader :: Expression RealSpace -> Maybe (Expression RealSpace) -> String -> String
generateHeader e arg n = "// -*- mode: C++; -*-\n\n#include \"MinimalFunctionals.h\"\n#include \"utilities.h\"\n#include \"handymath.h\"\n\n" ++ 
                     classCode e arg (n ++ "_type") ++
                     "\n\nFunctional " ++ n ++"(" ++ codeA arg ++ ") {\n\treturn Functional(new " ++ n ++ "_type(" ++ codeA' arg ++ "), \"" ++ n ++ "\");\n}\n"
    where codeA (Just rep) = "double " ++ code (makeHomogeneous rep)
          codeA Nothing = ""
          codeA' (Just rep) = code (makeHomogeneous rep)
          codeA' Nothing = ""



\end{code}
