\begin{code}
{-# LANGUAGE GADTs, PatternGuards #-}
module CodeGen ( RealSpace(..), r_var,
                 KSpace(..), k_var, kx, ky, kz, k, ksqr,
                 Scalar, s_var,
                 fft, ifft, integrate, grad, derive,
                 Expression,                 
                 Statement( (:=), (:?=) ),
                 Type,
                 generateHeader, code, latex, simp, setZero, codeStatement )
       where

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
data Scalar = Constant Double | 
              S String |
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
mapExpression f (a ::** n) = (mapExpression f a) ** mapExpression f n
mapExpression f (a ::+ b) = mapExpression f a + mapExpression f b
mapExpression f (a ::* b) = mapExpression f a * mapExpression f b
mapExpression f (a ::/ b) = mapExpression f a / mapExpression f b
mapExpression f (Expression x) = f x

setZero :: (Type a, Type b) => b -> Expression a -> Expression a
setZero v (Scalar e) = Scalar (setZero v e)
setZero v (Cos e) = cos (setZero v e)
setZero v (Sin e) = sin (setZero v e)
setZero v (Exp e) = exp (setZero v e)
setZero v (Log e) = log (setZero v e)
setZero v (Abs e) = abs (setZero v e)
setZero v (Signum e) = signum (setZero v e)
setZero v (a ::** n) = (setZero v a) ** setZero v n
setZero v (a ::+ b) = setZero v a + setZero v b
setZero v (a ::* b) = setZero v a * setZero v b
setZero v (a ::/ b) =
    if zb == 0
    then case isKSpace (Expression v) of
           Same -> case isKSpace a  of
                     Same -> setZero v (derive v 1 a / derive v 1 b)
                     _ -> error "oopsies"
           _ -> case isScalar (Expression v) of
                  Same -> case isScalar a  of
                            Same -> setZero v (derive v 1 a / derive v 1 b)
                            _ -> error "oopsies"
                  _ -> case isRealSpace (Expression v) of
                         Same -> case isRealSpace a  of
                                   Same -> setZero v (derive v 1 a / derive v 1 b)
                                   _ -> error "oopsies"
                         _ -> error "oops" -- setZero v (derive v 1 a / derive v 1 b)
                      else za / zb
    where za = setZero v a
          zb = setZero v b
setZero v (Expression x) = zeroHelper v x

kzeroValue :: Type a => KSpace -> Expression a
kzeroValue (K v) = s_var (v ++ "[0]")
kzeroValue Kx = 0
kzeroValue Ky = 0
kzeroValue Kz = 0
kzeroValue Delta = error "Can't find zero value of delta function!"
kzeroValue (FFT _) = error "Hard to find zero value of an fft"

instance Code Scalar where
  codePrec _ (S v) = showString v
  codePrec p (Constant d) = showsPrec p d
  codePrec _ (Integrate r) = showString "(" . codePrec 0 r . showString ") * gd.dvolume"
  latexPrec _ (S v) = showString v
  latexPrec p (Constant d) = showsPrec p d
  latexPrec _ (Integrate r) = showString "integrate(" . latexPrec 0 r . showString ")"
instance Type Scalar where
  toExpression = Expression . Constant . fromRational . toRational
  s_var = Expression . S
  isConstant (Expression (Constant x)) = Just x
  isConstant (Scalar x) = isConstant x
  isConstant _ = Nothing
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
--integrate (x ::+ y) = Expression (Integrate x) + Expression (Integrate y)
integrate x = Expression (Integrate x)

fft :: Expression RealSpace -> Expression KSpace
fft (Scalar e) = Scalar e * Expression Delta
fft r = Expression (FFT r)

ifft :: Expression KSpace -> Expression RealSpace
ifft ke = case factor (Expression Delta) ke of
            Just k' -> mapExpression kzeroValue k'
            Nothing -> Expression (IFFT ke)
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
                    Expression a ::** Expression a |
                    Expression a ::+ Expression a |
                    Expression a ::* Expression a |
                    Expression a ::/ Expression a
              deriving (Eq, Ord, Show)


infixr 8 ::**
infixl 7 ::*, ::/
infixl 6 ::+

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
codeE p (x ::** y) =
  case isConstant y of
    Just 0 -> showString "1" -- this shouldn't happen...
    Just 1 -> codeE p x
    Just 0.5 -> showString "sqrt(" . codePrec 0 x . showString ")"
    Just n
      | n > 0 && n == fromInteger (floor n) && 
        fromInteger (floor (n/2)) /= fromInteger(floor n)/(2 :: Double) ->
          codeE p (x ::* x ::** toExpression (n-1))
      | n > 0 && n == fromInteger (floor n) ->
          codeE p (x ::** toExpression (n/2) ::* x ::** toExpression (n/2))
      | n > 0 && n - fromInteger (floor n) == 0.5 ->
          codeE p (x ::** toExpression (floor n :: Integer) ::* x ::** 0.5)
      | n < 0 -> codeE p (1 ::/ x ::** toExpression (-n))
    _ -> showString "pow(" . codeE 0 x . showString ", " . codeE 0 y . showString ")"
codeE p (x ::+ y) = showParen (p > 6) (codeE 6 x . showString " + " . codeE 7 y)
codeE p (x ::* y) = showParen (p > 7) (codeE 7 x . showString " * " . codeE 8 y)
codeE p (x ::/ y) = showParen (p > 7) (codeE 7 x . showString " / " . codeE 8 y)


latexE :: (Type a, Code a) => Int -> Expression a -> ShowS
latexE p (Scalar x) = latexPrec p x
latexE p (Expression x) = latexPrec p x
latexE _ (Cos x) = showString "cos(" . latexE 0 x . showString ")"
latexE _ (Sin x) = showString "sin(" . latexE 0 x . showString ")"
latexE _ (Exp x) = showString "exp(" . latexE 0 x . showString ")"
latexE _ (Log x) = showString "log(" . latexE 0 x . showString ")"
latexE _ (Abs x) = showString "fabs(" . latexE 0 x . showString ")"
latexE _ (Signum _) = undefined
latexE _ (x ::** y) = latexE 8 x . showString "^{" . latexE 0 y . showString "}"
latexE p (x ::+ y) = showParen (p > 6) (latexE 6 x . showString " + " . latexE 7 y)
latexE p (x ::* y) = showParen (p > 7) (latexE 7 x . showString " * " . latexE 8 y)
latexE p (x ::/ y) = showParen (p > 7) (latexE 7 x . showString " / " . latexE 8 y)


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

class (Eq a, Show a, Code a) => Type a where 
  toExpression :: Real x => x -> Expression a
  toExpression = Scalar . Expression . Constant . fromRational . toRational
  isScalar :: Expression a -> Same a Scalar
  isScalar _ = Different
  isRealSpace :: Expression a -> Same a RealSpace
  isRealSpace _ = Different
  isKSpace :: Expression a -> Same a KSpace
  isKSpace _ = Different
  s_var :: String -> Expression a
  s_var = Scalar . s_var
  isConstant :: Expression a -> Maybe Double
  isConstant (Scalar x) = isConstant x
  isConstant _ = Nothing
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
  (+) = \x y -> case (x, y) of
                  _ | y == 0 -> x
                    | x == 0 -> y
                  (Scalar s1 ::* a, Scalar s2 ::* b) | s1 == s2 -> Scalar s1 * (a + b)
                  (Scalar a, Scalar b) -> Scalar (a+b)
                  _ | Just xx <- isConstant x, Just yy <- isConstant y -> toExpression (xx + yy)
                    | otherwise ->
                        case factor x y of
                         Just y' | x /= 1 -> (1+y')*x
                         _ ->
                           case factor y x of
                             Just x' | y /= 1 -> (1+x')*y
                             _ -> x ::+ y
  (-) = \x y -> x + (-y)
  -- the fromRational on the following line is needed to avoid a
  -- tricky infinite loop where -1 intepreted as (negate $ fromRational 1)
  negate = \x -> (fromRational $ -1) * x
  (*) = \x y -> case (x, y) of
                  _ | x == 1 -> y
                    | y == 1 -> x
                    | x == 0 || y == 0 -> 0
                    | x == y -> x ** 2
                    | Just a <- isConstant x, Just b <- isConstant y -> toExpression (a*b)
                  (Scalar a ::* c, Scalar b) -> Scalar (a*b) * c
                  (Scalar a, Scalar s ::* b) -> Scalar (a*s) * b
                  (Scalar s1 ::* a, Scalar s2 ::* b) -> Scalar (s1*s2) * a * b
                  (_, a ::* b) -> x * a * b
                  (_, a ::/ b) -> case factor x b of
                                  Just b' -> a / b'
                                  Nothing -> 
                                    case factor b x of
                                      Just x' -> x' * a
                                      Nothing -> (x*a)/b
                  (a ::/ b, _) -> case factor y b of
                                  Just b' -> a / b'
                                  Nothing -> 
                                    case factor b y of
                                      Just y' -> a * y'
                                      Nothing -> (a*y)/b
                  (Scalar a, Scalar b) -> Scalar (a*b)
                  (_, Scalar _) -> y * x
                  _ -> x ::* y
  fromInteger = \x -> toExpression (fromInteger x :: Rational)
  abs = undefined
  signum = undefined

-- | Factor is used to do some crude simplification in the * and /
-- operators.  I should note that we *don't* fully simplify ratios,
-- but I expect that this may be "good enough," whatever that means.

factor :: Type a => Expression a -> Expression a -> Maybe (Expression a)
factor x (y ::* z) | x == y = Just z
                   | x == z = Just y
                   | otherwise = case factor x y of
                                 Just y' -> Just (y' * z)
                                 _ -> case factor x z of
                                      Just z' -> Just (y * z')
                                      _ -> Nothing
factor x (y ::** n) | x == y = Just (y ** (n-1))
factor x y | x == y = Just 1
factor _ _ = Nothing

instance Type a => Fractional (Expression a) where
  (/) = \x y -> case (x, y) of
                  _ | y == 1 -> x
                    | x == 0 -> 0
                    | Just a <- isConstant x, Just b <- isConstant y -> toExpression (a/b)
                    | Just b <- isConstant y -> x * toExpression (1.0 / b)
                  (Scalar a, (Scalar s) ::* b) -> Scalar (a/s) / b
                  ((Scalar r) ::* a, (Scalar s) ::* b) -> Scalar (r/s) * a / b
                  (Scalar a, Scalar b) -> Scalar (a/b)
                  (a ::/ b, c ::/ d) ->
                    case factor b c of
                      Just c' -> a * c' / d
                      Nothing -> 
                        case factor d a of
                          Just a' -> a' * c / b
                          Nothing -> 
                            case factor c b of
                              Just b' -> a / d / b'
                              Nothing ->
                                case factor a d of
                                  Just d' -> c / b / d'
                                  Nothing -> (a * c) / (b * d)
                  (_, a ::/ b) -> 
                    case factor x a of
                      Just a' -> b / a'
                      Nothing ->
                        case factor a x of
                          Just x' -> x'*b
                          Nothing -> (x*b)/a
                  (a ::/ b, _) -> 
                    case factor a y of
                      Just y' -> 1/(b*y')
                      Nothing ->
                        case factor y a of
                          Just a' -> a' / b
                          Nothing -> a/(b*y)
                  _ -> x ::/ y
  fromRational = toExpression

instance Type a => Floating (Expression a) where
  pi = Scalar $ Expression $ Constant pi
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
  (a ::* b) ** y | a == b = a ** (2*y)
  (x ::** b) ** c = x ::** (b*c)
  x ** y | isConstant x == Just 0 = 0
         | isConstant x == Just 1 = 1
         | isConstant y == Just 0 = 1
         | isConstant y == Just 1 = x
         | Just xv <- isConstant x, Just yv <- isConstant y = toExpression (xv**yv)
  x ** y = x ::** y
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
simp (e1 ::+ e2) = (st1 >> st2, f1 + f2)
    where (st1, f1) = simp e1
          (st2, f2) = simp e2
simp (e1 ::* e2) = (st1 >> st2, f1 * f2)
    where (st1, f1) = simp e1
          (st2, f2) = simp e2
simp (e1 ::/ e2) = (st1 >> st2, f1 / f2)
    where (st1, f1) = simp e1
          (st2, f2) = simp e2
simp (e1 ::** e2) = (st1 >> st2, f1 ** f2)
    where (st1, f1) = simp e1
          (st2, f2) = simp e2
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
simp x = error ("I haven't finished simp" ++ show x)

derive :: (Type a, Type b) => b -> Expression a -> Expression a -> Expression b
derive v dda (x ::+ y) = derive v dda x + derive v dda y
derive v dda (x ::* y) = derive v (dda*x) y + derive v (dda*y) x
derive v dda (x ::/ y) = derive v (dda / y) x + derive v (-dda*x/y**2) y
derive v _ (Scalar x) = derive v 1 x -- FIXME
derive v dda (Cos e) = derive v (-dda*sin e) e
derive v dda (Sin e) = derive v (dda*cos e) e
derive v dda (Exp e) = derive v (dda*exp e) e
derive v dda (Log e) = derive v (dda/e) e
derive _ _ (Abs _) = error "I didn't think we'd need abs"
derive _ _ (Signum _) = error "I didn't think we'd need signum"
derive v dda (e ::** n) = derive v (dda*n*e**(n-1)) e
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
                functionCode "integral" "double" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\ndouble output=0;\n" ++ codeStatement codeIntegrate ++ "return output*gd.dvolume;\n" ) ++
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
                         "output" :+= e'
          where (st, e') = simp (integrate e)
      codeVTransform = do Initialize e
                          st
                          "output" := e'
          where (st, e') = simp e
      codeDTransform = do Initialize (makeHomogeneous e)
                          "output" := (makeHomogeneous e)
      codeDerive = do Initialize (makeHomogeneous (derive (R "x") 1 e))
                      "output" := makeHomogeneous (derive (R "x") 1 e)
      codeGrad = do st
                    "(*outgrad)" :+= e'
          where (st, e') = simp (derive (R "x") (r_var "ingrad") e )
      codeA (Just rep) = "(double " ++ code (makeHomogeneous rep) ++ "_arg) : " ++ code (makeHomogeneous rep) ++ "(" ++ code (makeHomogeneous rep) ++ "_arg)"
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
