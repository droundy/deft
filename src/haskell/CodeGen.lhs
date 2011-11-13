\begin{code}
{-# LANGUAGE GADTs #-}
module CodeGen ( RealSpace, r_var,
                 KSpace, k_var,
                 Scalar, s_var,
                 (**),
                 fft, ifft, integrate, grad, derive,
                 Expression, 
                 Statement( (:=), (:?=) ) ) 
       where

import Hash ( hash )
import Prelude hiding ((**))
\end{code}

The \verb!RealSpaceField! data type describes a field in real space.

\begin{code}
data RealSpace = R String |
                 IFFT (Expression KSpace)
               deriving ( Eq, Ord )
data KSpace = K String |
              Delta | -- handy for FFT of homogeneous systems
              FFT (Expression RealSpace)
            deriving ( Eq, Ord )
data Scalar = Constant Double | 
              S String |
              Integrate (Expression RealSpace)
            deriving ( Eq, Ord )

instance Show RealSpace where
  showsPrec _ (R v) = showString (v ++ "[i]")
  showsPrec _ (IFFT k) = showString "ifft(" . showsPrec 0 k . showString ")"
instance Type RealSpace where
  derivativeHelper = deriveR
  var v (Scalar _) = s_var v
  var v _ = r_var v
  prefix "" _ = "for (int i=0; i<gd.NxNyNz; i++) {    \n"
  prefix v _ = "Grid " ++ v ++ "(gd);\nfor (int i=0; i<gd.NxNyNz; i++) {    \n"
  postfix _ = "\n}\n"
deriveR :: String -> Expression RealSpace -> RealSpace -> Expression RealSpace
deriveR v dda (R v') | v == v' = dda
                     | otherwise = 0
deriveR v ddr (IFFT k) = derive v (fft ddr) k -- CHECK THIS!

instance Show KSpace where
  showsPrec _ (K v) = showString (v ++ "[i]")
  showsPrec _ Delta = showString "delta(k?)"
  showsPrec _ (FFT r) = showString "fft(" . showsPrec 0 r . showString ")"
instance Type KSpace where
  derivativeHelper = deriveK
  var v (Scalar _) = s_var v
  var v _ = k_var v
  prefix "" _ = "for (int i=0; i<gd.NxNyNz; i++) {    \n"
  prefix v _ = "ReciprocalGrid " ++ v ++ "(gd);\nfor (int i=0; i<gd.NxNyNz; i++) {    \n"
  postfix _ = "\n}\n"
deriveK :: String -> Expression KSpace -> KSpace -> Expression RealSpace
deriveK _ _ (K _) = 0
deriveK _ _ Delta = 0
deriveK v ddk (FFT r) = derive v (ifft ddk) r -- CHECK THIS!

mapExpression :: (Type a, Type b) => (a -> Expression b) -> Expression a -> Expression b
mapExpression _ (Scalar e) = Scalar e
mapExpression f (Cos e) = cos (mapExpression f e)
mapExpression f (Sin e) = sin (mapExpression f e)
mapExpression f (Exp e) = exp (mapExpression f e)
mapExpression f (Log e) = log (mapExpression f e)
mapExpression f (Abs e) = abs (mapExpression f e)
mapExpression f (Signum e) = signum (mapExpression f e)
mapExpression f (a ::** n) = (mapExpression f a) ** n
mapExpression f (a ::+ b) = mapExpression f a + mapExpression f b
mapExpression f (a ::* b) = mapExpression f a * mapExpression f b
mapExpression f (a ::/ b) = mapExpression f a / mapExpression f b
mapExpression f (Expression x) = f x

kzeroValue :: KSpace -> Expression RealSpace
kzeroValue (K v) = s_var (v ++ "[0]")
kzeroValue Delta = error "Can't find zero value of delta function!"
kzeroValue (FFT _) = error "Hard to find zero value of an fft"

instance Show Scalar where
  showsPrec _ (S v) = showString v
  showsPrec p (Constant d) = showsPrec p d
  showsPrec _ (Integrate r) = showString "integrate(" . showsPrec 0 r . showString ")"
instance Type Scalar where
  toExpression = Expression . Constant . fromRational
  s_var = Expression . S
  isConstant (Expression (Constant x)) = Just x
  isConstant (Scalar x) = isConstant x
  isConstant _ = Nothing
  derivativeHelper = deriveS
  var v _ = s_var v
  prefix "" _ = ""
  prefix x e = "double " ++ show (var x e) ++ ";\n"
deriveS :: String -> Expression Scalar -> Scalar -> Expression RealSpace
deriveS _ _ (S _) = 0
deriveS _ _ (Constant _) = 0
deriveS v dds (Integrate e) = derive v (Scalar dds*s_var "dV") e

r_var :: String -> Expression RealSpace
r_var v = Expression (R v)

k_var :: String -> Expression KSpace
k_var v = Expression (K v)

integrate :: Expression RealSpace -> Expression Scalar
integrate = Expression . Integrate

fft :: Expression RealSpace -> Expression KSpace
fft (Scalar e) = Scalar e * Expression Delta
fft r = Expression (FFT r)

ifft :: Expression KSpace -> Expression RealSpace
ifft k = case factor (Expression Delta) k of
           Just k' -> mapExpression kzeroValue k'
           Nothing -> Expression (IFFT k)
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
                    Expression a ::** Int |
                    Expression a ::+ Expression a |
                    Expression a ::* Expression a |
                    Expression a ::/ Expression a
              deriving (Eq, Ord)


infixr 8 ::**, **
infixl 7 ::*, ::/
infixl 6 ::+

instance (Eq a, Show a) => Show (Expression a) where
  showsPrec = showsE
showsE :: Show a => Int -> Expression a -> ShowS
showsE p (Scalar x) = showsPrec p x
showsE p (Expression x) = showsPrec p x
showsE _ (Cos x) = showString "cos(" . showsE 0 x . showString ")"
showsE _ (Sin x) = showString "sin(" . showsE 0 x . showString ")"
showsE _ (Exp x) = showString "exp(" . showsE 0 x . showString ")"
showsE _ (Log x) = showString "log(" . showsE 0 x . showString ")"
showsE _ (Abs x) = showString "fabs(" . showsE 0 x . showString ")"
showsE _ (Signum _) = undefined
showsE _ (_ ::** 0) = showString "0" -- this shouldn't happen...
showsE p (x ::** 1) = showsE p x       -- this also shouldn't happen
showsE p (x ::** n) | n `mod` 1 == 1 = showsE p (x ::* (x ::** (n-1)))
showsE p (x ::** n) = showsE p (x2 ::* x2)
  where x2 = x ::** (n `div` 2)
showsE p (x ::+ y) = showParen (p > 6) (showsE 6 x . showString " + " . showsE 7 y)
showsE p (x ::* y) = showParen (p > 7) (showsE 7 x . showString "*" . showsE 8 y)
showsE p (x ::/ y) = showParen (p > 7) (showsE 7 x . showString "/" . showsE 8 y)

class Pow p where
  (**) :: p -> Int -> p

class (Eq a, Show a) => Type a where 
  toExpression :: Rational -> Expression a
  toExpression = Scalar . Expression . Constant . fromRational
  s_var :: String -> Expression a
  s_var = Scalar . s_var
  isConstant :: Expression a -> Maybe Double
  isConstant (Scalar x) = isConstant x
  isConstant _ = Nothing
  derivativeHelper :: String -> Expression a -> a -> Expression RealSpace
  var :: String -> Expression a -> Expression a
  prefix :: String -> Expression a -> String
  postfix :: Expression a -> String
  postfix _ = ""

instance Type a => Pow (Expression a) where
  (**) = \x n -> 
    case (x, n) of
      (0, _) -> 0
      (_, 0) -> 1
      (_, 1) -> x
      _ -> x ::** n
instance Type a => Num (Expression a) where
  (+) = \x y -> case (x, y) of
                  _ | y == 0 -> x
                    | x == 0 -> y
                  (Scalar s1 ::* a, Scalar s2 ::* b) | s1 == s2 -> Scalar s1 * (a + b)
                  (Scalar a, Scalar b) -> Scalar (a+b)
                  _ -> 
                    case isConstant x of
                      Just xx ->
                        case isConstant y of
                          Just yy -> toExpression $ toRational (xx + yy)
                          _ -> x ::+ y
                      _ ->
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
  fromInteger = toExpression . fromInteger
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
factor x y | x == y = Just 1
factor _ _ = Nothing

instance Type a => Fractional (Expression a) where
  (/) = \x y -> case (x, y) of
                  _ | y == 1 -> x
                    | x == 0 -> 0
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
  asin = undefined
  acos = undefined
  atan = undefined
  asinh = undefined
  acosh = undefined
  atanh = undefined

-- | grad takes the gradient of a scalar-valued expression with
-- respect to a particular realspace variable.
grad :: String -> Expression Scalar -> Expression RealSpace
grad v e = derive v 1 e
-- grad v (Expression (Integrate (Expression (R v')))) 
--   | v == v' = s_var "dV"
--   | otherwise = 0
-- grad v (Expression (Integrate (x ::+ y))) =
--   grad v (integrate x) + grad v (integrate y)
-- grad v (Expression (Integrate (x ::** n))) =
--   (fromIntegral n)* x**(n-1) * grad v (integrate x)
-- grad v (Expression (Integrate (x ::* y))) =
--   x * grad v (integrate y) + y * grad v (integrate x)
-- grad v (Expression (Integrate (x ::/ y))) =
--   grad v (integrate x) / y - (x / y**2) * grad v (integrate y)
-- grad _ (Expression (Integrate (Expression (IFFT _)))) =
--   r_var "insert ifft grad here"
-- grad v (Expression (Integrate (Scalar e))) =
--   s_var "volume" * grad v e
-- grad v (x ::+ y) = grad v x + grad v y
-- grad v (x ::* y) = grad v x * Scalar y + Scalar x * grad v y
-- grad v (x ::/ y) = grad v x / Scalar y - Scalar (x / y**2) * grad v y
-- grad v (x ::** n) = Scalar ((fromIntegral n)*x**(n-2))*grad v x
-- grad v (Scalar e) = grad v e
-- grad _ (Expression (S _)) = 0
-- grad _ (Expression (Constant _)) = 0
-- grad _ x = error $ "need to implement grad " ++ show x

derive :: Type a => String -> Expression a -> Expression a -> Expression RealSpace
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
derive v dda (e ::** n) = derive v (dda*(fromIntegral n)*e**(n-1)) e
derive v dda (Expression e) = derivativeHelper v dda e
\end{code}

The statement type allows us to string together a sequence of
definitions and assignments.

\begin{code}
data Statement a where
     (:=) :: Type a => String -> Expression a -> Statement (Expression a)
     (:?=) :: Type a => String -> Expression a -> Statement (Expression a)
     Return :: a -> Statement a
     (:>>) :: Statement b -> Statement a -> Statement a
instance Monad Statement where
  (>>=) = thenS
  return = Return
instance Show (Statement a) where
  showsPrec = showsS
showsS :: Int -> Statement a -> ShowS
showsS _ (x := y) = showString (prefix "" y) . 
                    showsPrec 0 (var x y) . showString " = " . showsPrec 0 y . showString ";" .
                    showString (postfix y)
showsS _ (x :?= y) = showString (prefix x y) .
                     showsPrec 0 (var x y) . showString " := " . showsPrec 0 y . showString ";" .
                     showString (postfix y)
showsS _ (Return _) = id
showsS p (x :>> y) = showsS p x . showString "\n" . showsS p y

thenS :: Statement a -> (a -> Statement b) -> Statement b
thenS f g = f :>> g (evalS f)

evalS :: Statement a -> a
evalS (s := e) = var s e
evalS (s :?= e) = var (s ++ "-" ++ hash (show e)) e
evalS (Return a) = a
evalS (_ :>> b) = evalS b
\end{code}
