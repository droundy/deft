\begin{code}
{-# LANGUAGE GADTs, PatternGuards #-}
module CodeGen ( RealSpace, r_var,
                 KSpace, k_var, kx, ky, kz, k, ksqr,
                 Scalar, s_var,
                 fft, ifft, integrate, grad, derive,
                 Expression,                 
                 Statement( (:=), (:?=) ),
                 generateHeader, code, codeDouble, codeVector ) 
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
  codePrec _ (IFFT ke) = showString "ifft(" . showsPrec 0 ke . showString ")"
  codeDouble (R v) = v
  codeDouble (IFFT ke) = "ifft(" ++ codeDouble ke ++ ")"
instance Type RealSpace where
  derivativeHelper = deriveR
  var v (Scalar _) = s_var v
  var v _ = r_var v
  prefix "" _ = "for (int i=0; i<gd.NxNyNz; i++) {    \n"
  prefix v _ = "Grid " ++ v ++ "(gd);\nfor (int i=0; i<gd.NxNyNz; i++) {    \n"
  postfix _ = "\n}\n"
  toScalar (R v) = s_var v
  toScalar (IFFT ke) = makeHomogeneous ke
deriveR :: String -> Expression RealSpace -> RealSpace -> Expression RealSpace
deriveR v dda (R v') | v == v' = dda
                     | otherwise = 0
deriveR v ddr (IFFT ke) = derive v (fft ddr) ke -- CHECK THIS!

instance Code KSpace where
  codePrec _ (K v) = showString (v ++ "[i]")
  codePrec _ Kx = showString "khere[0]"
  codePrec _ Ky = showString "khere[1]"
  codePrec _ Kz = showString "khere[2]"
  codePrec _ Delta = showString "delta(k?)"
  codePrec _ (FFT r) = showString "fft(" . showsPrec 0 r . showString ")"
instance Type KSpace where
  derivativeHelper = deriveK
  var v (Scalar _) = s_var v
  var v _ = k_var v
  prefix "" _ = "for (int i=0; i<gd.NxNyNz; i++) {    \n"
  prefix v _ = "ReciprocalGrid " ++ v ++ "(gd);\nfor (int i=0; i<gd.NxNyNz; i++) {    \n"
  postfix _ = "\n}\n"
  toScalar (K v) = s_var v
  toScalar Delta = 1
  toScalar Kx = 0
  toScalar Ky = 0
  toScalar Kz = 0
  toScalar (FFT e) = makeHomogeneous e
deriveK :: String -> Expression KSpace -> KSpace -> Expression RealSpace
deriveK _ _ (K _) = 0
deriveK _ _ Kx = 0
deriveK _ _ Ky = 0
deriveK _ _ Kz = 0
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
mapExpression f (a ::** n) = (mapExpression f a) ** mapExpression f n
mapExpression f (a ::+ b) = mapExpression f a + mapExpression f b
mapExpression f (a ::* b) = mapExpression f a * mapExpression f b
mapExpression f (a ::/ b) = mapExpression f a / mapExpression f b
mapExpression f (Expression x) = f x

kzeroValue :: KSpace -> Expression RealSpace
kzeroValue (K v) = s_var (v ++ "[0]")
kzeroValue Kx = 0
kzeroValue Ky = 0
kzeroValue Kz = 0
kzeroValue Delta = error "Can't find zero value of delta function!"
kzeroValue (FFT _) = error "Hard to find zero value of an fft"

instance Code Scalar where
  codePrec _ (S v) = showString v
  codePrec p (Constant d) = showsPrec p d
  codePrec _ (Integrate r) = showString "output += (" . codePrec 0 r . showString ") * gd.dvolume"
instance Type Scalar where
  toExpression = Expression . Constant . fromRational . toRational
  s_var = Expression . S
  isConstant (Expression (Constant x)) = Just x
  isConstant (Scalar x) = isConstant x
  isConstant _ = Nothing
  derivativeHelper = deriveS
  var v _ = s_var v
  prefix "" _ = ""
  prefix x e = "double " ++ code (var x e) ++ ";\n"
  toScalar (Integrate r) = makeHomogeneous r
  toScalar v = Expression v
  
deriveS :: String -> Expression Scalar -> Scalar -> Expression RealSpace
deriveS _ _ (S _) = 0
deriveS _ _ (Constant _) = 0
deriveS v dds (Integrate e) = derive v (Scalar dds*s_var "dV") e

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
                   
  codeDouble (x ::+ y) = "(" ++ codeDouble x ++ ") + (" ++ codeDouble y ++ ")"
  codeDouble (x ::* y) = "(" ++ codeDouble x ++ ") * (" ++ codeDouble y ++ ")"
  codeDouble (x ::/ y) = "(" ++ codeDouble x ++ ") / (" ++ codeDouble y ++ ")"
  codeDouble (x ::** y) = 
    case isConstant y of
      Just 1 -> codeDouble x
      Just 0.5 -> "sqrt(" ++ codeDouble x ++ ")"
      Just n
        | n > 0 && n == fromInteger (floor n) && 
          fromInteger (floor (n/2)) /= fromInteger(floor n)/(2 :: Double) ->
            codeDouble (x ::* x ::** toExpression (n-1))
        | n > 0 && n == fromInteger (floor n) ->
            codeDouble (x ::** toExpression (n/2) ::* x ::** toExpression (n/2))
        | n > 0 && n - fromInteger (floor n) == 0.5 ->
            codeDouble (x ::** toExpression (floor n :: Integer) ::* x ::** 0.5)
        | n < 0 -> codeDouble (1 ::/ x ::** toExpression (-n))
      _ -> "pow(" ++ codeDouble x ++ ", " ++ codeDouble y ++ ")"
  codeDouble (Expression x) = codeDouble x
  codeDouble (Log x) = "log(" ++ codeDouble x ++ ")"
  codeDouble x = code x

  codeVector (Expression x) = "for (int i=0; i<gd.NxNyNz; i++) {\n\t" ++ code x ++ ";\n}\n"
  codeVector (Log x) = "for (int i=0; i<gd.NxNyNz; i++) {\n\toutput[i] = " ++ code (Log x) ++ ";\n}\n"
  codeVector (x ::+ y) = "for (int i=0; i<gd.NxNyNz; i++) {\n\toutput[i] = " ++ code (x+y) ++ ";\n}\n"
  codeVector (x ::* y) = "for (int i=0; i<gd.NxNyNz; i++) {\n\toutput[i] = " ++ code (x*y) ++ ";\n}\n"
  codeVector (x ::/ y) = "for (int i=0; i<gd.NxNyNz; i++) {\n\toutput[i] = " ++ code (x/y) ++ ";\n}\n"
  codeVector (x ::** y) = "for (int i=0; i<gd.NxNyNz; i++) {\n\toutput[i] = " ++ code (x**y) ++ ";\n}\n"
  codeVector x = code x

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

class Code a  where
    codePrec  :: Int -> a -> ShowS
    codePrec _ x s = code x ++ s
    code      :: a -> String 
    code x = codePrec 0 x ""
    codeDouble :: a -> String
    codeDouble x = code x
    codeVector :: a -> String
    codeVector x = code x

class (Eq a, Show a) => Type a where 
  toExpression :: Real x => x -> Expression a
  toExpression = Scalar . Expression . Constant . fromRational . toRational
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
  toScalar :: a -> Expression Scalar

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
  (**) = \x y ->
    case isConstant x of
      Just 0 -> 0
      Just 1 -> 1
      _ -> case isConstant y of
             Just 0 -> 1
             Just 1 -> x
             Just 2 -> x ::* x
             _ -> case x of
                    e ::** n -> e ** (n * y)
                    _ -> x ::** y
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
derive v dda (e ::** n) = derive v (dda*n*e**(n-1)) e
derive v dda (Expression e) = derivativeHelper v dda e
\end{code}

The statement type allows us to string together a sequence of
definitions and assignments.

\begin{code}
data Statement a where
     (:=) :: (Type a, Code a) => String -> Expression a -> Statement (Expression a)
     (:?=) :: (Type a, Code a) => String -> Expression a -> Statement (Expression a)
     Return :: a -> Statement a
     (:>>) :: Statement b -> Statement a -> Statement a
infixl 1 :>>
infix 4 :=, :?=
instance Monad Statement where
  (>>=) = thenS
  return = Return
instance Show (Statement a) where
    showsPrec = showsS
showsS :: Int -> Statement a -> ShowS
showsS _ (x := y) = showString x . showString " := " . showsPrec 0 y
showsS _ (x :?= y) = showString x . showString " :?= " . showsPrec 0 y
showsS _ (Return _) = showString "return ???" -- . showsPrec 0 x
showsS p (x :>> y) = showsS p x . showString "\n" . showsS p y

instance Code (Statement a) where
	codePrec = codeS
codeS :: Int -> Statement a -> ShowS
codeS _ (x := y) = showString (prefix "" y) . 
                    codePrec 0 (var x y) . showString " = " . codePrec 0 y . showString ";" .
                    showString (postfix y)
codeS _ (x :?= y) = showString (prefix x y) .
                     codePrec 0 (var x y) . showString " := " . codePrec 0 y . showString ";" .
                     showString (postfix y)
codeS _ (Return _) = id
codeS p (x :>> y) = codeS p x . showString "\n" . codeS p y

thenS :: Statement a -> (a -> Statement b) -> Statement b
thenS f g = f :>> g (evalS f)

evalS :: Statement a -> a
evalS (s := e) = var s e
evalS (s :?= e) = var (s ++ "-" ++ hash (show e)) e
evalS (Return a) = a
evalS (_ :>> b) = evalS b
\end{code}

\begin{code}
functionCode :: String -> String -> [(String, String)] -> String -> String
functionCode "" "" [] "" = ""
functionCode "" "" (x:xs) "" = if xs == [] 
                               then fst x ++ " " ++ snd x
                               else fst x ++ " " ++ snd x ++ ", " ++ functionCode "" "" xs ""
functionCode n t a b = t ++ " " ++ n ++ "(" ++ functionCode "" "" a "" ++ ") const {\n" ++ b ++ "\n}\n"



classCode :: Expression RealSpace -> String -> String
classCode e n = "class " ++ n ++ " : public FunctionalInterface {\npublic:\n" ++ n ++ "()  {\n\thave_integral = true;\n}\n" ++
                functionCode "I_have_analytic_grad" "bool" [] "\treturn false;" ++
                functionCode "integral" "double" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n\tdouble output = 0;\n\t" ++ codeVector (integrate e) ++ "\n\treturn output;\n")  ++
                functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n\tVectorXd output(gd.NxNyNz);\n\t" ++ codeVector e ++ "\n\treturn output;\n")  ++
                functionCode "transform" "double" [("double", "kT"), ("double", "x")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(x==x); // to avoid an unused parameter error\n\treturn " ++ codeDouble e ++ ";\n") ++
                functionCode "append_to_name" "bool" [("", "std::string")] "\treturn false;" ++
                functionCode "derive" "double" [("double", "kT"), ("double", "x")] ("\tassert(kT==kT);\n\tassert(x==x);\n\treturn " ++ codeDouble (derive "x" 1 e) ++ ";\n") ++
                functionCode "d_by_dT" "double" [("double", ""), ("double", "")] "\tassert(0); // fail\n\treturn 0;\n" ++
                functionCode "derive_homogeneous" "Expression" [("const Expression &", "")] "\tassert(0); // fail\n\treturn Expression(0);\n" ++
                functionCode "grad" "Functional" [("const Functional", "&ingrad"), ("const Functional", "&x"), ("bool", "")] "\tassert(&ingrad==&ingrad);\n\tassert(&x==&x);\n\treturn ingrad;" ++
                functionCode "grad_T" "Functional" [("const Functional", "&ingradT")] "\tassert(&ingradT==&ingradT);\n\treturn ingradT;" ++
                functionCode "grad" "void" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x"), ("const VectorXd", "&ingrad"), ("VectorXd", "*outgrad"), ("VectorXd", "*outpgrad")] ("\tassert(kT==kT); // to avoid an unused parameter error\n\tassert(&gd); // to avoid an unused parameter error\n\tassert(&x); // to avoid an unused parameter error\n\tassert(outpgrad==outpgrad);\n\tfor (int i=0; i<gd.NxNyNz; i++) {\n\t\t(*outgrad)[i] += (" ++ code (derive "x" 1 e) ++ ") * ingrad[i];\n\t}" ) ++
                functionCode "printme" "Expression" [("const Expression", "&x")] ("\treturn funexpr(\"" ++ n ++ "()\")(x);") ++
                functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")] "\tFunctionalInterface::print_summary(prefix, energy, name);" ++
                "private:\n}; // End of " ++ n ++ " class"

generateHeader :: Expression RealSpace -> String -> String
generateHeader e n = "// -*- mode: C++; -*-\n\n#pragma once\n\n#include \"MinimalFunctionals.h\"\n#include \"utilities.h\"\n#include \"handymath.h\"\n\n" ++ 
                     classCode e (n ++ "_type") ++
                     "\n\ninline Functional " ++ n ++"() {\n\treturn Functional(new " ++ n ++ "_type(), \"" ++ n ++ "\");\n}\n"

\end{code}
