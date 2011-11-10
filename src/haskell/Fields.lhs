\begin{code}
module Fields ( RealSpaceField, r_var, 
                ReciprocalSpaceField ) 
       where

import Prelude hiding ((**))
\end{code}

The \verb!RealSpaceField! data type describes a field in real space.

\begin{code}
data RealSpaceField = R_scalar Double |
                      R_variable String |
                      (:+) RealSpaceField RealSpaceField |
                      (:*) RealSpaceField RealSpaceField |
                      (:/) RealSpaceField RealSpaceField |
                      RF Function RealSpaceField |
                      IFFT ReciprocalSpaceField
                    deriving ( Eq, Ord, Show )

infixl 7 :*, :/
infixl 6 :+

r_var :: String -> RealSpaceField
r_var = R_variable

fft :: RealSpaceField -> ReciprocalSpaceField
fft (R_scalar a :* b) = K_scalar a * (fft b)
fft a = FFT a

ifft :: ReciprocalSpaceField -> RealSpaceField
ifft (K_scalar a :*: b) = R_scalar a * (ifft b)
ifft a = IFFT a

instance Num RealSpaceField where
  (+) = \x y -> case (x, y) of
                  (R_scalar 0, _) -> y
                  (_, R_scalar 0) -> x
                  (R_scalar (-1) :* a, R_scalar (-1) :* b) -> negate (a + b)
                  (R_scalar a, R_scalar b) -> R_scalar (a+b)
                  _ -> x :+ y
  (-) = \x y -> x + (-y)
  negate = \x -> case x of
    R_scalar s -> R_scalar (-s)
    (R_scalar s :* a) -> R_scalar (-s) * a
    _ -> -1 * x
  (*) = \x y -> case (x, y) of
                  (R_scalar 1, _) -> y
                  (_, R_scalar 1) -> x
                  (R_scalar 0, _) -> R_scalar 0
                  (_, R_scalar 0) -> R_scalar 0
                  (R_scalar a :* c, R_scalar b) -> R_scalar (a*b) * c
                  (R_scalar a, R_scalar s :* b) -> R_scalar (a*s) * b
                  (R_scalar s1 :* a, R_scalar s2 :* b) -> R_scalar (s1*s2) * a * b
                  (_, a :/ b) -> (x*a)/b
                  (a :/ b, _) -> (a*y)/b
                  (R_scalar a, R_scalar b) -> R_scalar (a*b)
                  (_, R_scalar _) -> y * x
                  _ -> x :* y
  fromInteger = R_scalar . fromInteger
  abs = undefined
  signum = undefined

instance Fractional RealSpaceField where
  (/) = \x y -> case (x, y) of
                  (_, R_scalar 1) -> x
                  (R_scalar 0, _) -> R_scalar 0
                  (R_scalar a, R_scalar s :* b) -> R_scalar (a/s) / b
                  (R_scalar r :* a, R_scalar s :* b) -> R_scalar (r/s) * a / b
                  (R_scalar a, R_scalar b) -> R_scalar (a/b)
                  (a :/ b, c :/ d) -> (a * c) / (b * d)
                  (_, a :/ b) -> (x*b)/a
                  (a :/ b, _) -> a/(b*y)
                  _ -> x :/ y
  fromRational = R_scalar . fromRational

instance Floating RealSpaceField where
  pi = R_scalar pi
  exp = \x -> case x of
    R_scalar 0 -> R_scalar 1
    _ -> RF Exp x
  log = \x -> case x of
    R_scalar 1 -> R_scalar 0
    _ -> RF Log x
  sinh = \x -> case x of
    R_scalar 0 -> R_scalar 0
    _ -> 0.5*(exp x - exp (-x))
  cosh = \x -> case x of
    R_scalar 0 -> R_scalar 1
    _ -> 0.5*(exp x + exp (-x))
  sin = \x -> case x of
    R_scalar 0 -> R_scalar 0
    _ -> RF Sin x
  cos = undefined
  asin = undefined
  acos = undefined
  atan = undefined
  asinh = undefined
  acosh = undefined
  atanh = undefined
instance Pow RealSpaceField where
  (**) = \x n -> 
    case (x, n) of
      (_, 0) -> R_scalar 1
      (_, 1) -> x
      (R_scalar s, _) -> R_scalar (s ^^ n)
      _ -> RF (Pow n) x
\end{code}

The \verb!ReciprocalSpaceField! data type describes a field in real space.

\begin{code}
data ReciprocalSpaceField = K_scalar Double |
                            K_variable String |
                            (:+:) ReciprocalSpaceField ReciprocalSpaceField |
                            (:*:) ReciprocalSpaceField ReciprocalSpaceField |
                            (:/:) ReciprocalSpaceField ReciprocalSpaceField |
                            KF Function ReciprocalSpaceField |
                            FFT RealSpaceField
                          deriving ( Eq, Ord, Show )

infixl 7 :*:, :/:
infixl 6 :+:

instance Num ReciprocalSpaceField where
  (+) = \x y -> case (x, y) of
                  (K_scalar 0, _) -> y
                  (_, K_scalar 0) -> x
                  (K_scalar (-1) :*: a, K_scalar (-1) :*: b) -> negate (a + b)
                  (K_scalar a, K_scalar b) -> K_scalar (a+b)
                  _ -> x :+: y
  (-) = \x y -> x + (-y)
  negate = \x -> case x of
    K_scalar s -> K_scalar (-s)
    (K_scalar s :*: a) -> K_scalar (-s) * a
    _ -> -1 * x
  (*) = \x y -> case (x, y) of
                  (K_scalar 1, _) -> y
                  (_, K_scalar 1) -> x
                  (K_scalar 0, _) -> K_scalar 0
                  (_, K_scalar 0) -> K_scalar 0
                  (K_scalar a :*: c, K_scalar b) -> K_scalar (a*b) * c
                  (K_scalar a, K_scalar s :*: b) -> K_scalar (a*s) * b
                  (K_scalar s1 :*: a, K_scalar s2 :*: b) -> K_scalar (s1*s2) * a * b
                  (_, a :/: b) -> (x*a)/b
                  (a :/: b, _) -> (a*y)/b
                  (K_scalar a, K_scalar b) -> K_scalar (a*b)
                  (_, K_scalar _) -> y * x
                  _ -> x :*: y
  fromInteger = K_scalar . fromInteger
  abs = undefined
  signum = undefined
  
instance Fractional ReciprocalSpaceField where
  (/) = \x y -> case (x, y) of
                  (_, K_scalar 1) -> x
                  (K_scalar 0, _) -> K_scalar 0
                  (K_scalar a, K_scalar s :*: b) -> K_scalar (a/s) / b
                  (K_scalar r :*: a, K_scalar s :*: b) -> K_scalar (r/s) * a / b
                  (K_scalar a, K_scalar b) -> K_scalar (a/b)
                  (a :/: b, c :/: d) -> (a * c) / (b * d)
                  (_, a :/: b) -> (x*b)/a
                  (a :/: b, _) -> a/(b*y)
                  _ -> x :/: y
  fromRational = K_scalar . fromRational

instance Floating ReciprocalSpaceField where
  pi = K_scalar pi
  exp = \x -> case x of
    K_scalar 0 -> K_scalar 1
    _ -> KF Exp x
  log = \x -> case x of
    K_scalar 1 -> K_scalar 0
    _ -> KF Log x
  sinh = \x -> case x of
    K_scalar 0 -> K_scalar 0
    _ -> 0.5*(exp x - exp (-x))
  cosh = \x -> case x of
    K_scalar 0 -> K_scalar 1
    _ -> 0.5*(exp x + exp (-x))
  sin = undefined
  cos = undefined
  asin = undefined
  acos = undefined
  atan = undefined
  asinh = undefined
  acosh = undefined
  atanh = undefined
instance Pow ReciprocalSpaceField where
  (**) = \x n -> 
    case (x, n) of
      (_, 0) -> K_scalar 1
      (_, 1) -> x
      (K_scalar s, _) -> K_scalar (s ^^ n)
      _ -> KF (Pow n) x
\end{code}

The \verb!Function! data type describes a simple function of one or
more variable.

\begin{code}
data Function = Constant Double |
                Variable String |
                (:.) Function Function | -- compose two functions
                Cos |
                Sin |
                Exp |
                Log |
                Abs |
                Signum |
                Pow Int |
                (::+) Function Function |
                (::*) Function Function |
                (::/) Function Function
              deriving (Show, Eq, Ord)

infixr 9 :.
infixr 8 **
infixl 7 ::*, ::/
infixl 6 ::+

class Pow p where
  (**) :: p -> Int -> p

instance Pow Function where
  (**) = \x n -> 
    case (x, n) of
      (Constant 0, _) -> Constant 0
      (_, 0) -> Constant 1
      (_, 1) -> x
      _ -> Pow n :. x
instance Num Function where
  (+) = \x y -> case (x, y) of
                  (Constant 0, _) -> y
                  (_, Constant 0) -> x
                  (Constant (-1) ::* a, Constant (-1) ::* b) -> negate (a + b)
                  (Constant a, Constant b) -> Constant (a+b)
                  _ -> x ::+ y
  (-) = \x y -> x + (-y)
  negate = \x -> case x of
    Constant s -> Constant (-s)
    (Constant s ::* a) -> Constant (-s) * a
    _ -> -1 * x
  (*) = \x y -> case (x, y) of
                  (Constant 1, _) -> y
                  (_, Constant 1) -> x
                  (Constant 0, _) -> Constant 0
                  (_, Constant 0) -> Constant 0
                  (Constant a ::* c, Constant b) -> Constant (a*b) * c
                  (Constant a, Constant s ::* b) -> Constant (a*s) * b
                  (Constant s1 ::* a, Constant s2 ::* b) -> Constant (s1*s2) * a * b
                  (_, a ::/ b) -> (x*a)/b
                  (a ::/ b, _) -> (a*y)/b
                  (Constant a, Constant b) -> Constant (a*b)
                  (_, Constant _) -> y * x
                  _ -> x ::* y
  fromInteger = Constant . fromInteger
  abs = undefined
  signum = undefined
  
instance Fractional Function where
  (/) = \x y -> case (x, y) of
                  (_, Constant 1) -> x
                  (Constant 0, _) -> Constant 0
                  (Constant a, Constant s ::* b) -> Constant (a/s) / b
                  (Constant r ::* a, Constant s ::* b) -> Constant (r/s) * a / b
                  (Constant a, Constant b) -> Constant (a/b)
                  (a ::/ b, c ::/ d) -> (a * c) / (b * d)
                  (_, a ::/ b) -> (x*b)/a
                  (a ::/ b, _) -> a/(b*y)
                  _ -> x ::/ y
  fromRational = Constant . fromRational

instance Floating Function where
  pi = Constant pi
  exp = \x -> case x of
    Constant 0 -> Constant 1
    _ -> Exp :. x
  log = \x -> case x of
    Constant 1 -> Constant 0
    _ -> Log :. x
  sinh = \x -> case x of
    Constant 0 -> Constant 0
    _ -> 0.5*(exp x - exp (-x))
  cosh = \x -> case x of
    Constant 0 -> Constant 1
    _ -> 0.5*(exp x + exp (-x))
  sin = undefined
  cos = undefined
  asin = undefined
  acos = undefined
  atan = undefined
  asinh = undefined
  acosh = undefined
  atanh = undefined
\end{code}
