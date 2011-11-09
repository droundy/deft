\begin{code}
module Fields ( RealSpaceField, r_var, 
                ReciprocalSpaceField ) 
       where
\end{code}

The \verb!RealSpaceField! data type describes a field in real space.

\begin{code}
data RealSpaceField = R_scalar Double |
                      R_variable String |
                      (:+) RealSpaceField RealSpaceField |
                      (:-) RealSpaceField RealSpaceField |
                      (:*) RealSpaceField RealSpaceField |
                      (:/) RealSpaceField RealSpaceField |
                      (:**) RealSpaceField Int |
                      R_neg RealSpaceField |
                      R_abs RealSpaceField |
                      R_exp RealSpaceField |
                      R_log RealSpaceField
                    deriving ( Eq, Ord, Show )

infixr 8 :**
infixl 7 :*, :/
infixl 6 :+, :-

r_var :: String -> RealSpaceField
r_var = R_variable

r_neg :: RealSpaceField -> RealSpaceField
r_neg (R_neg a) = a
r_neg (a :- b) = b :- a
r_neg (R_scalar a) = R_scalar (-a)
r_neg x = R_neg x

instance Num RealSpaceField where
  (+) = \x y -> case (x, y) of
                  (R_scalar 0, _) -> y
                  (_, R_scalar 0) -> x
                  (R_neg a, R_neg b) -> r_neg (a + b)
                  (_, R_neg b) -> x - b
                  (R_neg a, _) -> a - y
                  (R_scalar a, R_scalar b) -> R_scalar (a+b)
                  _ -> x :+ y
  (-) = \x y -> case (x, y) of
                  (R_scalar 0, _) -> r_neg y
                  (_, R_scalar 0) -> x
                  (R_neg a, R_neg b) -> b - a
                  (_, R_neg b) -> x + b
                  (R_neg a, _) -> r_neg (a+y)
                  (R_scalar a, R_scalar b) -> R_scalar (a-b)
                  _ -> x :- y
  (*) = \x y -> case (x, y) of
                  (R_scalar 1, _) -> y
                  (_, R_scalar 1) -> x
                  (R_scalar 0, _) -> R_scalar 0
                  (_, R_scalar 0) -> R_scalar 0
                  (R_neg a, R_scalar b) -> R_scalar (-b) * a
                  (R_scalar a, R_neg b) -> R_scalar (-a) * b
                  (R_neg a, R_neg b) -> a * b
                  (_, R_neg b) -> r_neg (x*b)
                  (R_neg a, _) -> r_neg (a*y)
                  (R_scalar a, R_scalar b) -> R_scalar (a*b)
                  _ -> x :* y
  fromInteger = R_scalar . fromInteger
  abs = R_abs
  signum = undefined
\end{code}

The \verb!ReciprocalSpaceField! data type describes a field in real space.

\begin{code}
data ReciprocalSpaceField = K_scalar Double |
                      K_variable String |
                      (:+:) ReciprocalSpaceField ReciprocalSpaceField |
                      (:-:) ReciprocalSpaceField ReciprocalSpaceField |
                      (:*:) ReciprocalSpaceField ReciprocalSpaceField |
                      (:/:) ReciprocalSpaceField ReciprocalSpaceField |
                      (:**:) ReciprocalSpaceField Int |
                      K_neg ReciprocalSpaceField |
                      K_exp ReciprocalSpaceField |
                      K_log ReciprocalSpaceField
                    deriving ( Eq, Ord, Show )

infixr 8 :**:
infixl 7 :*:, :/:
infixl 6 :+:, :-:

\end{code}
