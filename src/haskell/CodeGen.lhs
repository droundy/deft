\begin{code}
{-# LANGUAGE GADTs, PatternGuards #-}
module CodeGen ( RealSpace(..), r_var,
                 KSpace(..), k_var, kx, ky, kz, k, ksqr,
                 Scalar, s_var,
                 fft, ifft, integrate, grad, derive,
                 Expression,                 
                 Statement( (:=), (:?=) ),
                 Type, 
                 makeHomogeneous, isConstant, -- for debugging only!!!
                 generateHeader, code, latex, simp, setZero, codeStatement )
       where

import Debug.Trace

import qualified Data.Map as Map
import Hash ( hash )
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
  latexPrec _ (R v) = showString v
  latexPrec _ (IFFT ke) = showString "ifft(" . latexPrec 0 ke . showString ")"
instance Type RealSpace where
  isRealSpace _ = Same
  derivativeHelper v ddr r | Same <- isRealSpace (Expression v), v == r = ddr
  derivativeHelper v ddr (IFFT ke) = derive v (fft ddr) ke
  derivativeHelper _ _ _ = 0
  zeroHelper v x | Same <- isRealSpace (Expression v), v == x = 0
  zeroHelper v (IFFT ke) = ifft (setZero v ke)
  zeroHelper _ x = Expression x
  simpHelper (R x) = (return (), r_var x)
  simpHelper (IFFT ksp@(Expression (K _))) = (do st
                                                 Initialize (r_var temp)
                                                 temp := ifft ksp'
                                                 return ()
                                            , r_var temp)
      where (st, ksp') = simp ksp
            Expression (R temp) = evalS ("tempIFFT" :?= ifft ksp')
  simpHelper (IFFT ksp) = (do st
                              Initialize (k_var tempk)
                              tempk := ksp'
                              Initialize (r_var temp)
                              temp := ifft ksp''
                              return ()
                         , r_var temp)
      where (st, ksp') = simp ksp
            ksp''@(Expression (K tempk)) = evalS ("tempR" :?= ksp')
            Expression (R temp) = evalS ("tempIFFT" :?= ifft ksp'')
  var v (Scalar _) = s_var v
  var v _ = r_var v
  prefix "" (Expression (IFFT _)) = ""
  prefix "" _ = "for (int i=0; i<gd.NxNyNz; i++) {    \n"
  prefix v _ = "Grid " ++ v ++ "(gd);\nfor (int i=0; i<gd.NxNyNz; i++) {    \n"
  postfix (Expression (IFFT _)) = ""
  postfix _ = "\n}\n"
  codeStatementHelper a op (Expression (IFFT (Expression (K v)))) = a ++ op ++ "ifft(gd, " ++ v ++ ");\n"
  codeStatementHelper _ _ (Expression (IFFT _)) = error "It is a bug to generate code for a non-var input to ifft"
  codeStatementHelper a op e = prefix "" e ++ a ++ "[i]" ++ op ++ code e ++ ";\n" ++ postfix e
  initialize (Expression (R x)) = "VectorXd " ++ x ++ "(gd.NxNyNz);\n"
  initialize _ = "VectorXd output(gd.NxNyNz);\n"
  toScalar (R v) = s_var v
  toScalar (IFFT ke) = makeHomogeneous $ setZero (S "dr") $ setZero Kx $ setZero Ky $ setZero Kz ke


instance Code KSpace where
  codePrec _ (K v) = showString (v ++ "[i]")
  codePrec _ Kx = showString "k_i[0]"
  codePrec _ Ky = showString "k_i[1]"
  codePrec _ Kz = showString "k_i[2]"
  codePrec _ Delta = showString "delta(k?)"
  codePrec _ (FFT r) = showString "fft(gd, " . codePrec 0 (makeHomogeneous r) . showString ")"
  latexPrec _ (K v) = showString v
  latexPrec _ Kx = showString "kx"
  latexPrec _ Ky = showString "ky"
  latexPrec _ Kz = showString "kz"
  latexPrec _ Delta = showString "delta(k)"
  latexPrec _ (FFT r) = showString "fft(" . latexPrec 0 r . showString ")"
instance Type KSpace where
  isKSpace _ = Same
  derivativeHelper v ddk kk | Same <- isKSpace (Expression v), kk == v = ddk
  derivativeHelper v ddk (FFT r) = derive v (ifft ddk) r
  derivativeHelper _ _ _ = 0
  zeroHelper v x | Same <- isKSpace (Expression v), x == v = 0
  zeroHelper v (FFT r) = fft (setZero v r)
  zeroHelper _ x = Expression x
  simpHelper (FFT rsp@(Expression (R _))) = (do st
                                                Initialize (k_var temp)
                                                temp := fft rsp'
                                                return ()
                                            , k_var temp)
      where (st, rsp') = simp rsp
            Expression (K temp) = evalS ("tempFFT" :?= fft rsp')
  simpHelper (FFT rsp) = (do st
                             Initialize (r_var tempr)
                             tempr := rsp'
                             Initialize (k_var temp)
                             temp := fft rsp''
                             return ()
                         , k_var temp)
      where (st, rsp') = simp rsp
            rsp''@(Expression (R tempr)) = evalS ("tempR" :?= rsp')
            Expression (K temp) = evalS ("tempFFT" :?= fft rsp'')
  simpHelper ksp = (return (), Expression ksp)
  var v (Scalar _) = s_var v
  var v _ = k_var v
  prefix "" (Expression (FFT _)) = ""
  prefix "" _ = "for (int i=0; i<gd.NxNyNzOver2; i++) {\n\tconst int z = i % gd.NzOver2;\n\tconst int n = (i-z)/gd.NzOver2;\n\tconst int y = n % gd.Ny;\n\tconst int x = (n-y)/gd.Ny;\n\tconst RelativeReciprocal rvec((x>gd.Nx/2) ? x - gd.Nx : x, (y>gd.Ny/2) ? y - gd.Ny : y, z);\n\tconst Reciprocal k_i = gd.Lat.toReciprocal(rvec);\n\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);\n\t"
  prefix v _ = "ReciprocalGrid " ++ v ++ "(gd);\nfor (int i=0; i<gd.NxNyNz; i++) {\n\t"
  postfix (Expression (FFT _)) = ""
  postfix _ = "}\n"
  codeStatementHelper a op (Expression (FFT (Expression (R v)))) = a ++ op ++ "fft(gd, " ++ v ++ ");\n"
  codeStatementHelper _ _ (Expression (FFT _)) = error "It is a bug to generate code for a non-var input to fft"
  codeStatementHelper a op e = prefix "" e ++ "if (i == 0) {\n\t" ++ a ++ "[0]" ++ op ++ code (setZero Kz (setZero Ky (setZero Kx e))) ++ ";\n\t}\n\telse {\n\t" ++ a ++ "[i]"  ++ op ++ code e ++ ";\n\t}" ++ postfix e
  initialize (Expression (K x)) = "VectorXcd " ++ x ++ "(gd.NxNyNz);\n"
  initialize _ = "VectorXcd output(gd.NxNyNz);\n"
  toScalar (K v) = s_var v
  toScalar Delta = 1
  toScalar Kx = 0
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
mapExpression f (Product p) = trace "mapExpression Product" product $ map ff $ product2list p
  where ff (e,n) = (mapExpression f e) ** toExpression n
mapExpression f (Sum s) = trace "mapExpression Sum" pairs2sum $ map ff $ sum2list_pair s
  where ff (x,y) = (x, mapExpression f y)
mapExpression f (Expression x) = f x

setZero :: (Type a, Type b) => b -> Expression a -> Expression a
setZero v (Scalar e) = Scalar (setZero v e)
setZero v (Cos e) = cos (setZero v e)
setZero v (Sin e) = sin (setZero v e)
setZero v (Exp e) = exp (setZero v e)
setZero v (Log e) = log (setZero v e)
setZero v (Abs e) = abs (setZero v e)
setZero v (Signum e) = signum (setZero v e)
setZero v (Product p) | product2denominator p == [] = product $ map ff $ product2list p
  where ff (e,n) = (setZero v e) ** toExpression n
setZero v (Product p) = 
  if trace ("setZero Product " ++ latex (Product p)) zd /= 0
  then zn / zd
  else if trace "L'Hopital's rule" zn /= 0
       then error "L'Hopital's rule failure"
       else               case isKSpace (Expression v) of
                            Same -> 
                              case isKSpace n of
                                Same -> trace ("dtop is " ++ latex dtop ++ " dbot is " ++ latex dbot) 
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
                                        Same -> setZero v (derive v 1 n / derive v 1 d)
                                        _ -> error "oopsies"
                                    _ -> error "oops" -- setZero v (derive v 1 a / derive v 1 b)
  where d = product $ product2denominator p
        n = product $ product2numerator p
        zn = setZero v n
        zd = setZero v d
setZero v (Sum s) = pairs2sum $ map sz $ sum2list_pair s
  where sz (f,x) = (f, setZero v x)
setZero v (Expression x) = zeroHelper v x

instance Code Scalar where
  codePrec _ (S v) = showString v
  codePrec _ (Integrate r) = showString "(" . codePrec 0 r . showString ") * gd.dvolume"
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
  simpHelper sc = (return (), Expression sc)
  var v _ = s_var v
  codeStatementHelper a " = " (Expression (Integrate e)) = a ++ " = 0;\nfor (int i=0; i<gd.NxNyNz; i++) {    \n" ++ a ++ "+=" ++ code e ++ ";\n}\n"
  codeStatementHelper a " += " (Expression (Integrate e)) = "for (int i=0; i<gd.NxNyNz; i++) {    \n" ++ a ++ " += " ++ code e ++ ";\n}\n"
  codeStatementHelper a op e = a ++ op ++ code e ++ ";\n"
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
fft r = Expression (FFT r)

ifft :: Expression KSpace -> Expression RealSpace
ifft ke = Expression (IFFT ke)
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

sum2list :: Type a => Map.Map (Expression a) Double -> [Expression a]
sum2list s = map toe $ sum2list_pair s
  where toe (f,x) = Sum $ Map.singleton x f

sum2list_pair :: Type a => Map.Map (Expression a) Double -> [(Double, Expression a)]
sum2list_pair s = map rev $ Map.assocs s
  where rev (a,b) = (b,a)

pairs2sum :: Type a => [(Double, Expression a)] -> Expression a
pairs2sum s = helper $ filter ((/= 0) . snd) $ filter ((/= 0) . fst) s
  where rev (a,b) = (b,a)
        helper [] = 0
        helper [(1,e)] = e
        helper [(x,Sum y)] | [(x2,y')] <- sum2list_pair y = helper [(x*x2, y')]
        helper es = Sum $ Map.fromList $ map rev es

product2list :: Type a => Map.Map (Expression a) Double -> [(Expression a, Double)]
product2list s = Map.assocs s

prod :: Type a => Map.Map (Expression a) Double -> Expression a
prod p | Map.size p == 1 = case product2list p of [(e,1)] -> e
                                                  _ -> Product p
prod p = Product p

list2product :: Type a => [(Expression a, Double)] -> Expression a
list2product [(e,1)] = e
list2product xs = prod $ Map.fromList xs

product2numerator :: Type a => Map.Map (Expression a) Double -> [Expression a]
product2numerator s = map f $ product2numerator_pair s
  where f (a,b) = a ** (toExpression b)

product2numerator_pair :: Type a => Map.Map (Expression a) Double -> [(Expression a, Double)]
product2numerator_pair s = filter ((>=0) . snd) $ product2list s

product2denominator :: Type a => Map.Map (Expression a) Double -> [Expression a]
product2denominator s = map n $ filter ((<0) . snd) $ product2list s
  where n (a,b) = a **(toExpression $ -b)

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
codeE pree (Product p) = showParen (pree > 7) $
                         case den of
                         [] -> codesimple num
                         [a] -> codesimple num . showString "/" . codeE 8 a
                         _ -> codesimple num . showString " / ( " . codeE 0 (product den) . showString " )"
  where codesimple [] = showString "1"
        codesimple [(a,n)] = codee a n
        codesimple [(a,n),(b,m)] = codee a n . showString "*" . codee b m
        codesimple ((a,n):es) = codee a n . showString "*" . codesimple es
        num = product2numerator_pair p
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
codeE p (Sum s) = showParen (p > 6) (showString me)
  where me = foldl addup "" $ sum2list_pair s
        addup "" (1,e) = codeE 6 e ""
        addup "" (f,e) = if e == 1
                         then show f
                         else show f ++ "*" ++ codeE 6 e ""
        addup rest (1,e) = codeE 6 e (showString " + " $ rest)
        addup rest (f,e) = show f ++ "*" ++ codeE 6 e (showString " + " $ rest)

latexE :: (Type a, Code a) => Int -> Expression a -> ShowS
latexE p (Scalar x) = latexPrec p x
latexE p (Expression x) = latexPrec p x
latexE _ (Cos x) = showString "cos(" . latexE 0 x . showString ")"
latexE _ (Sin x) = showString "sin(" . latexE 0 x . showString ")"
latexE _ (Exp x) = showString "exp(" . latexE 0 x . showString ")"
latexE _ (Log x) = showString "log(" . latexE 0 x . showString ")"
latexE _ (Abs x) = showString "fabs(" . latexE 0 x . showString ")"
latexE _ (Signum _) = undefined
latexE p (Product x) | Map.size x == 1 && length (product2denominator x) == 0 =
  case product2list x of
    [(_,0)] -> showString "1" -- this shouldn't happen...
    [(_, n)] | n < 0 -> error "shouldn't have negative power here"
    [(e, 1)] -> latexE p e
    [(e, 0.5)] -> showString "sqrt(" . latexE 0 e . showString ")"
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
latexE pree (Product p) = showParen (pree > 7) $
                          case den of
                          [] -> ltexsimple num
                          a -> showString "\\frac{" . ltexsimple num . showString "}{" . 
                                 latexE 7 (product a) . showString "}"
  where ltexsimple [] = showString "1"
        ltexsimple [a,b] = latexE 7 a . showString "*" . latexE 7 b
        ltexsimple (a:es) = latexE 7 a . showString "*" . ltexsimple es
        num = product2numerator p
        den = product2denominator p
latexE p (Sum s) = showParen (p > 6) (showString me)
  where me = foldl addup "" $ sum2list_pair s
        addup "" (1,e) = latexE 6 e ""
        addup "" (f,e) = if e == 1
                         then show f
                         else show f ++ " " ++ latexE 6 e ""
        addup rest (f,e) = show f ++ " " ++ latexE 6 e (showString " + " $ rest)


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
isConstant (Sum s) = case sum2list_pair s of
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
  simpHelper :: a -> (Statement (), Expression a)
  var :: String -> Expression a -> Expression a
  codeStatementHelper :: String -> String -> Expression a -> String
  prefix :: String -> Expression a -> String
  postfix :: Expression a -> String
  postfix _ = ""
  initialize :: Expression a -> String
  initialize _ = ""
  toScalar :: a -> Expression Scalar

codeStatement :: Statement a -> String
codeStatement (a := e) = codeStatementHelper a " = " e
codeStatement (a :?= e) = codeStatementHelper (a ++ "_" ++ hash (show e)) " = " e
codeStatement (a :+= e) = codeStatementHelper a " += " e
codeStatement (Initialize e) = initialize e
codeStatement (Return _) = ""
codeStatement (a :>> b) = codeStatement a ++ codeStatement b


makeHomogeneous :: Type a => Expression a -> Expression Scalar
makeHomogeneous = mapExpression toScalar

instance Type a => Num (Expression a) where
  x + y | Just 0 == isConstant x = y
        | Just 0 == isConstant y = x
  Sum a + Sum b = sumup a $ sum2list_pair b
      where sumup x [] = Sum x
            sumup x ((f,y):ys) = case Map.lookup y x of
                                 Nothing -> sumup (Map.insert y f x) ys
                                 Just f' -> if f + f' == 0
                                            then sumup (Map.delete y x) ys
                                            else sumup (Map.insert y (f+f') x) ys
  Sum a + b = case Map.lookup b a of
                Just fac -> if fac + 1 == 0
                            then 
                              case Map.size a of 
                                1 -> 0
                                2 -> head $ sum2list deleted
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
  Sum x * y | Just c <- isConstant y = pairs2sum $ map (f c) $ sum2list_pair x
                      where f c (a,b) = (c*a, b)
  y * Sum x | Just c <- isConstant y = pairs2sum $ map (f c) $ sum2list_pair x
                      where f c (a,b) = (c*a, b)
  Sum x * y | [(f,a)] <- sum2list_pair x = pairs2sum [(f, a*y)]
  y * Sum x | [(f,a)] <- sum2list_pair x = pairs2sum [(f, a*y)]
  Product a * Product b = puttogether a (product2list b)
      where puttogether x [] = prod x
            puttogether x ((y,n):ys) =
                  case Map.lookup y x of
                    Just n' -> if n + n' == 0
                               then puttogether (Map.delete y x) ys
                               else puttogether (Map.insert y (n+n') x) ys
                    Nothing -> puttogether (Map.insert y n x) ys
  Product a * b = case Map.lookup b a of
                    Just n -> if n + 1 == 0
                             then prod deleted
                             else prod $ Map.insert b (n+1) a
                    Nothing -> prod $ Map.insert b 1 a
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
  (Product x) ** c | Just n <- isConstant c = list2product $ map (p n) $ product2list x
                       where p n (e,n2) = (e,n2*n)
  x ** c | Just n <- isConstant c = list2product [(x,n)]
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

simp :: Type a => Expression a -> (Statement (), Expression a)
simp (Expression e) = simpHelper e
simp (Sum s) = (sequence_ $ map fst simped, pairs2sum $ map snd simped)
  where es = sum2list_pair s
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
simp (Product p) = (sequence_ $ map fst simped, product $ map snd simped)
  where es = product2list p
        simped = map simpme es
        simpme (x,n) = (st, x' ** toExpression n)
          where (st,x') = simp x
simp _ = error "simp incomplete"

derive :: (Type a, Type b) => b -> Expression a -> Expression a -> Expression b
derive v dda (Sum s) = sum $ map dbythis $ sum2list_pair s
  where dbythis (f,x) = toExpression f * derive v dda x
derive v dda (Product p) = sum $ map dbythis $ product2list p
  where dbythis (x,n) = derive v (Product p*toExpression n*dda/x) x
derive v _ (Scalar x) = derive v 1 x -- FIXME
derive v dda (Cos e) = derive v (-dda*sin e) e
derive v dda (Sin e) = derive v (dda*cos e) e
derive v dda (Exp e) = derive v (dda*exp e) e
derive v dda (Log e) = derive v (dda/e) e
derive _ _ (Abs _) = error "I didn't think we'd need abs"
derive _ _ (Signum _) = error "I didn't think we'd need signum"
derive v dda (Expression e) = derivativeHelper v dda e
\end{code}

The statement type allows us to string together a sequence of
definitions and assignments.

\begin{code}
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
                functionCode "integral" "double" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\ndouble output=0;\n" ++ trace "hi there" codeStatement codeIntegrate ++ "return output*gd.dvolume;\n" ) ++
                functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n" ++ codeStatement codeVTransform ++ ";\nreturn output;\n")  ++
                functionCode "transform" "double" [("double", "kT"), ("double", "x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(x==x); // to avoid an unused parameter error\n" ++ code codeDTransform ++ ";\nreturn output;\n") ++
                functionCode "append_to_name" "bool" [("", "std::string")] "\treturn false;" ++
                functionCode "derive" "double" [("double", "kT"), ("double", "x")] ("\tassert(kT==kT);\n\tassert(x==x);\n" ++ code codeDerive ++ ";\nreturn output;\n") ++
                functionCode "d_by_dT" "double" [("double", ""), ("double", "")] "\tassert(0); // fail\n\treturn 0;\n" ++
                functionCode "derive_homogeneous" "Expression" [("const Expression &", "")] "\tassert(0); // fail\n\treturn Expression(0);\n" ++
                functionCode "grad" "Functional" [("const Functional", "&ingrad"), ("const Functional", "&x"), ("bool", "")] "\tassert(&ingrad==&ingrad);\n\tassert(&x==&x);\n\treturn ingrad;" ++
                functionCode "grad_T" "Functional" [("const Functional", "&ingradT")] "\tassert(&ingradT==&ingradT);\n\treturn ingradT;" ++
                functionCode "grad" "void" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x"), ("const VectorXd", "&ingrad"), ("VectorXd", "*outgrad"), ("VectorXd", "*outpgrad")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n\tassert(outpgrad==outpgrad);\n" ++ codeStatement codeGrad) ++
                functionCode "printme" "Expression" [("const Expression", "&x")] ("\treturn funexpr(\"" ++ n ++ "()\")(x);") ++
                functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")] "\tFunctionalInterface::print_summary(prefix, energy, name);" ++
                "private:\n"++ codeArgInit arg  ++"}; // End of " ++ n ++ " class"
    where
      codeIntegrate = do st
                         "output" :+= trace "outputting codeIntegrate stuff" e'
          where (st, e') = trace "codeIntegrate" simp (integrate e)
      codeVTransform = do Initialize e
                          st
                          "output" := e'
          where (st, e') = trace "codeVTransform" simp e
      codeDTransform = do Initialize (trace "codeDTransform" makeHomogeneous e)
                          "output" := (trace "codeDTransform" makeHomogeneous e)
      codeDerive = do Initialize (trace "codeDerive" makeHomogeneous (derive (R "x") 1 e))
                      "output" := trace "codeDerive" makeHomogeneous (derive (R "x") 1 e)
      codeGrad = do st
                    "(*outgrad)" :+= e'
          where (st, e') = trace "codeGrad" simp (derive (R "x") (r_var "ingrad") e )
      codeA (Just rep) = "(double " ++ trace "codeA" code (makeHomogeneous rep) ++ "_arg) : " ++ trace "codeA" code (makeHomogeneous rep) ++ "(" ++ code (makeHomogeneous rep) ++ "_arg)"
      codeA Nothing = "()"
      codeArgInit (Just rep) = code (Initialize (makeHomogeneous rep))
      codeArgInit Nothing = ""
      


generateHeader :: Expression RealSpace -> Maybe (Expression RealSpace) -> String -> String
generateHeader e arg n = "// -*- mode: C++; -*-\n\n#pragma once\n\n#include \"MinimalFunctionals.h\"\n#include \"utilities.h\"\n#include \"handymath.h\"\n\n" ++ 
                     classCode e arg (n ++ "_type") ++
                     "\n\ninline Functional " ++ n ++"(" ++ codeA arg ++ ") {\n\treturn Functional(new " ++ n ++ "_type(" ++ codeA' arg ++ "), \"" ++ n ++ "\");\n}\n"
    where codeA (Just rep) = "double " ++ code (makeHomogeneous rep)
          codeA Nothing = ""
          codeA' (Just rep) = code (makeHomogeneous rep)
          codeA' Nothing = ""
\end{code}
