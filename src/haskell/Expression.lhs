\begin{code}
{-# LANGUAGE PatternGuards, Rank2Types #-}

module Expression (Exprn(..),
                   break_real_from_imag,
                   RealSpace(..), r_var, dV, dVscalar, dr,
                   rx, ry, rz, rmag, rvec,
                   lat1, lat2, lat3, rlat1, rlat2, rlat3, volume,
                   KSpace(..), k_var, imaginary, kvec, kx, ky, kz, k, ksqr, setkzero,
                   Scalar(..), scalar,
                   Vector, t_var, cross, dot, (/.), (*.), (.*), xhat, yhat, zhat,
                   erf, heaviside,
                   nameVector, vfft, vifft,
                   fft, ifft, integrate, grad, derive, scalarderive, deriveVector, realspaceGradient,
                   Expression(..), joinFFTs, (===), var, protect, vprotect,
                   Type(..), Code(..), IsTemp(..),
                   makeHomogeneous, isConstant,
                   setZero, cleanvars, factorize, factorOut,
                   initializeE, freeE, newinitializeE, newfreeE,
                   nameE, newdeclareE,
                   sum2pairs, pairs2sum, codeStatementE, newcodeStatementE,
                   product2pairs, pairs2product, product2denominator,
                   hasActualFFT, hasFFT, hasexpression, hasExprn,
                   searchExpression, searchExpressionDepthFirst,
                   findRepeatedSubExpression, findNamedScalars, findOrderedInputs, findInputs,
                   findTransforms, transform, Symmetry(..),
                   MkBetter(..), Monoid(..), mconcat,
                   countexpression, substitute, countAfterRemoval,
                   substituteE, countAfterRemovalE,
                   mapExprn,
                   mapExpressionShortcut, -- just to avoid unused warning
                   countVars, varSet)
    where

import Debug.Trace

import qualified Data.Map as Map
import qualified Data.Set as Set
import LatexDouble ( latexDouble )

data Symmetry = Spherical { dk, kmax :: Double,
                            rresolution, rmax :: Expression Scalar } |
                VectorS { dk, kmax :: Double,
                          rresolution, rmax :: Expression Scalar }
                deriving ( Eq, Ord, Show )

data RealSpace = IFFT (Expression KSpace)
               | Rx | Ry | Rz
               deriving ( Eq, Ord, Show )
data KSpace = Delta | -- handy for FFT of homogeneous systems
              Complex (Expression Scalar) (Expression Scalar) |
              Kx | Ky | Kz |
              SetKZeroValue (Expression KSpace) (Expression KSpace) |
              SphericalFourierTransform Symmetry (Expression Scalar) |
              FFT (Expression RealSpace)
            deriving ( Eq, Ord, Show )
data Scalar = Summate (Expression RealSpace)
            deriving ( Eq, Ord, Show )

data Vector a = Vector (Expression a) (Expression a) (Expression a)
     deriving ( Eq, Ord, Show )

kinversion :: Expression KSpace -> Expression KSpace
kinversion (Var _ _ _ _ (Just e)) = kinversion e
kinversion e@(Var _ _ _ _ Nothing) = e
kinversion (Scalar e) = Scalar e
kinversion (Heaviside e) = heaviside (kinversion e)
kinversion (Cos e) = cos (kinversion e)
kinversion (Sin e) = sin (kinversion e)
kinversion (Erf e) = erf (kinversion e)
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

xhat, yhat, zhat :: Type a => Vector a
xhat = Vector 1 0 0
yhat = Vector 0 1 0
zhat = Vector 0 0 1

vector :: Type a => Expression a -> Expression a -> Expression a -> Vector a
vector a b c = Vector a b c

cross :: Type a => Vector a -> Vector a -> Vector a
cross (Vector x y z) (Vector x' y' z') =
  vector (y*z' - z*y') (z*x' - x*z') (x*y' - y*x')

dot :: Type a => Vector a -> Vector a -> Expression a
dot (Vector a b c) (Vector x y z) = a*x + b*y + c*z

(/.) :: Type a => Vector a -> Expression a -> Vector a
Vector x y z /. s = vector (x/s) (y/s) (z/s)

(*.) :: Type a => Vector a -> Expression a -> Vector a
Vector x y z *. s = vector (x*s) (y*s) (z*s)

(.*) :: Type a => Expression a -> Vector a -> Vector a
s .* Vector x y z = vector (x*s) (y*s) (z*s)
infixl 7 `cross`, `dot`, /., .*, *.

instance Type a => Code (Vector a) where
  codePrec _ (Vector a b c) = showString ("<" ++ code a ++ ", " ++ code b ++ ", " ++ code c ++">")
  latexPrec p (Vector a b c) = showParen (p > 6) (showString $
                                                  codePrec 6 a "" ++ " \\mathbf{\\hat{x}} + " ++
                                                  codePrec 6 b "" ++ " \\mathbf{\\hat{y}} + " ++
                                                  codePrec 6 c "" ++ " \\mathbf{\\hat{z}}")
instance Type a => Num (Vector a) where
  Vector x y z + Vector x' y' z' = vector (x+x') (y+y') (z+z')
  Vector x y z - Vector x' y' z' = vector (x-x') (y-y') (z-z')
  negate (Vector x y z) = Vector (-x) (-y) (-z)
  fromInteger 0 = Vector 0 0 0
  fromInteger _ = error "cannot convert non-zero integer to Vector!"
  abs _ = error "cannot use abs on Vector (it has wrong type in Haskell!)"
  signum x = x /. sqrt (x `dot` x)
  (*) = error "cannot use * on Vectors (use dot or cross instead!)"

instance Code RealSpace where
  codePrec _ (IFFT (Var _ _ ksp _ Nothing)) = showString ("ifft(gd, " ++ksp++ ")")
  codePrec _ (IFFT ke) = showString "ifft(gd, " . codePrec 0 ke . showString ")"
  codePrec _ Rx = showString "rx"
  codePrec _ Ry = showString "ry"
  codePrec _ Rz = showString "rz"
  newcodePrec _ (IFFT (Var _ _ ksp _ Nothing)) = showString ("ifft(" ++ksp++ ")")
  newcodePrec _ (IFFT ke) = showString "ifft(" . codePrec 0 ke . showString ")"
  newcodePrec _ Rx = showString "rx"
  newcodePrec _ Ry = showString "ry"
  newcodePrec _ Rz = showString "rz"
  latexPrec _ (IFFT ke) = showString "\\text{ifft}\\left(" . latexPrec 0 ke . showString "\\right)"
  latexPrec _ Rx = showString "\\textbf{r_x}"
  latexPrec _ Ry = showString "\\textbf{r_y}"
  latexPrec _ Rz = showString "\\textbf{r_z}"
instance Type RealSpace where
  amRealSpace _ = True
  mkExprn = ER
  derivativeHelper v ddr (IFFT ke) = derive v (fft ddr) (kinversion ke)
  derivativeHelper _ _ Rx = 0
  derivativeHelper _ _ Ry = 0
  derivativeHelper _ _ Rz = 0
  scalarderivativeHelper v (IFFT ke) = ifft (scalarderive v ke)
  scalarderivativeHelper _ Rx = 0
  scalarderivativeHelper _ Ry = 0
  scalarderivativeHelper _ Rz = 0
  zeroHelper v (IFFT ke) = ifft (setZero v ke)
  zeroHelper _ Rx = rx
  zeroHelper _ Ry = ry
  zeroHelper _ Rz = rz
  codeStatementHelper (Var _ _ a _ _) op (Expression (IFFT (Var _ _ v _ Nothing))) =
    a ++ op ++ "ifft(gd, " ++ v ++ ");\n"
  codeStatementHelper _ _ (Expression (IFFT e)) =
    error ("It is a bug to generate code for a non-var input to ifft\n"++ latex e)
  codeStatementHelper a op (Var _ _ _ _ (Just e)) = codeStatementHelper a op e
  codeStatementHelper a op e =
    unlines ["for (int i=0; i<gd.NxNyNz; i++) {",
             codes (1 :: Int) e,
             "\t}"]
      where codes n x = case findRepeatedSubExpression x of
              MB (Just (_,x')) -> "\t\tconst double t"++ show n ++ " = " ++ code x' ++ ";\n" ++
                                  codes (n+1) (substitute x' (s_var ("t"++show n)) x)
              MB Nothing -> "\t\t" ++ code a ++ op ++ code x ++ ";"
  --
  newcodeStatementHelper (Var _ _ a _ _) op (Expression (IFFT (Var _ _ v _ Nothing))) =
    a ++ op ++ "ifft(Nx,Ny,Nz,dV," ++ v ++ ");\n"
  newcodeStatementHelper _ _ (Expression (IFFT e)) =
    error ("It is a bug to generate newcode for a non-var input to ifft\n"++ latex e)
  newcodeStatementHelper a op (Var _ _ _ _ (Just e)) = newcodeStatementHelper a op e
  newcodeStatementHelper a op e =
    unlines ["for (int i=0; i<Nx*Ny*Nz; i++) {",
             newcodes (1 :: Int) e,
             "\t}"]
      where newcodes n x = case findRepeatedSubExpression x of
              MB (Just (_,x')) -> "\t\tconst double t"++ show n ++ " = " ++ newcode x' ++ ";\n" ++
                                  newcodes (n+1) (substitute x' (s_var ("t"++show n)) x)
              MB Nothing -> "\t\t" ++ newcode a ++ op ++ newcode (cleanvars x) ++ ";"
  initialize (Var IsTemp _ x _ Nothing) = "VectorXd " ++ x ++ "(gd.NxNyNz);"
  initialize _ = error "VectorXd output(gd.NxNyNz);"
  free (Var IsTemp _ x _ Nothing) = x ++ ".resize(0); // Realspace"
  free e = error $ trace "free error" ("free error " ++ show e)
  newdeclare _ = "Vector"
  newinitialize (Var _ _ x _ Nothing) = "Vector " ++ x ++ "(Nx*Ny*Nz); // RS"
  newinitialize x = error ("oops newinitializeE: " ++ show x)
  newfree (Var IsTemp _ x _ Nothing) = x ++ ".free(); // Realspace"
  newfree e = error $ trace "free error" ("free error " ++ show e)
  toScalar (IFFT ke) = makeHomogeneous ke
  toScalar Rx = 0
  toScalar Ry = 0
  toScalar Rz = 0
  mapExpressionHelper' f (IFFT ke) = ifft (f ke)
  mapExpressionHelper' _ Rx = rx
  mapExpressionHelper' _ Ry = ry
  mapExpressionHelper' _ Rz = rz
  subAndCountHelper x y (IFFT ke) = case subAndCount x y ke of (ke', n) -> (ifft ke', n)
  subAndCountHelper _ _ Rx = (rx, 0)
  subAndCountHelper _ _ Ry = (ry, 0)
  subAndCountHelper _ _ Rz = (rz, 0)
  searchHelper f (IFFT e) = f e
  searchHelper _ Rx = mempty
  searchHelper _ Ry = mempty
  searchHelper _ Rz = mempty
  joinFFThelper (Sum s0 _) = joinup Map.empty $ sum2pairs s0
        where joinup m [] = sum $ map toe $ Map.toList m
                where toe (rs, Right ks) = rs * ifft (joinFFTs $ pairs2sum ks)
                      toe (_, Left (_,e)) = e
              joinup m ((f,e):es) =
                case isifft e of
                  Nothing -> toExpression f * e + joinup m es
                  Just (ks,rs) -> joinup (Map.insert rs ks' m) es
                    where ks' = case Map.lookup rs m of
                                Nothing -> Left ([(f, ks)],
                                                 toExpression f * e)
                                Just (Left (ks0,_)) -> Right $ (f, ks) : ks0
                                Just (Right ks0) -> Right $ (f, ks) : ks0
  joinFFThelper e = e
  safeCoerce a _ = case mkExprn a of
                    ER a' -> Just a'
                    _ -> Nothing

-- In the following, we assume that when we have k == 0, then any
-- imaginary terms will vanish.  This simplifies the situation for odd
-- functions (whose fourier transforms are pure imaginary), since our
-- setZero sometimes has trouble employing L'Hopital's rule.
setKequalToZero :: Expression KSpace -> Expression KSpace
setKequalToZero e = setZero (EK kz) $ expand kz $ setZero (EK ky) $ setZero (EK kx) $
                    setZero (EK imaginary) $
                    -- trace ("setKequalToZero\n    "++code e)
                    e

instance Code KSpace where
  codePrec _ Kx = showString "k_i[0]"
  codePrec _ Ky = showString "k_i[1]"
  codePrec _ Kz = showString "k_i[2]"
  codePrec _ Delta = showString "delta(k?)"
  codePrec p (SetKZeroValue _ val) = codePrec p val
  codePrec _ (FFT r) = showString "fft(gd, " . codePrec 0 (makeHomogeneous r) . showString ")"
  codePrec _ (SphericalFourierTransform s r) =
    showString ("transform(" ++ show s ++ ",") . codePrec 0 (makeHomogeneous r) . showString ")"
  codePrec _ (Complex 0 1) = showString "complex(0,1)"
  codePrec _ (Complex a b) = showString "complex(" . codePrec 0 a . showString ", " .
                                                     codePrec 0 b . showString ")"
  latexPrec _ Kx = showString "k_{x}"
  latexPrec _ Ky = showString "k_{y}"
  latexPrec _ Kz = showString "k_{z}"
  latexPrec _ Delta = showString "\\delta(k)"
  latexPrec p (SetKZeroValue _ val) = latexPrec p val
  latexPrec _ (FFT r) = showString "\\text{fft}\\left(" . latexPrec 0 r . showString "\\right)"
  latexPrec _ (SphericalFourierTransform s r) =
    showString ("\\mathcal{F}(" ++ show s ++ ",") . latexPrec 0 r . showString ")"
  latexPrec p (Complex a b) = latexPrec p (a + s_var "i" * b)
instance Type KSpace where
  amKSpace _ = True
  mkExprn = EK
  derivativeHelper v ddk (FFT r) = derive v (ifft ddk) r
  derivativeHelper v ddk (SphericalFourierTransform _ e) =
    if derive v 1 e == 0 || ddk == 0
    then 0
    else error "do not have an implementation for derivativeHelper of SphericalFourierTransform"
  derivativeHelper v ddk (SetKZeroValue _ e) = derive v (setkzero 0 ddk) e -- FIXME: how best to handle k=0 derivative?
  derivativeHelper _ _ Kx = 0
  derivativeHelper _ _ Ky = 0
  derivativeHelper _ _ Kz = 0
  derivativeHelper _ _ Delta = 0
  derivativeHelper v ddk (Complex a b) = derive v (real_part ddk) a -
                                         derive v (imag_part ddk) b
  scalarderivativeHelper v (FFT r) = fft (scalarderive v r)
  scalarderivativeHelper v (SphericalFourierTransform s e) = transform s (scalarderive v e)
  scalarderivativeHelper v (SetKZeroValue z e) = setkzero (scalarderive v z) (scalarderive v e)
  scalarderivativeHelper v (Complex a b) = complex (scalarderive v a) (scalarderive v b)
  scalarderivativeHelper _ Kx = 0
  scalarderivativeHelper _ Ky = 0
  scalarderivativeHelper _ Kz = 0
  scalarderivativeHelper _ Delta = 0
  zeroHelper v (FFT r) = fft (setZero v r)
  zeroHelper v (SphericalFourierTransform s e) = transform s (setZero v e)
  zeroHelper _ Kx = Expression Kx
  zeroHelper _ Ky = Expression Ky
  zeroHelper _ Kz = Expression Kz
  zeroHelper _ Delta = Expression Delta
  zeroHelper v (Complex a b) = complex (setZero v a) (setZero v b)
  zeroHelper v e@(SetKZeroValue val _) | v == EK kz = val -- slightly hokey... assuming that if we set kz = 0 then we set kx and ky = 0
                                       | otherwise = Expression e
  codeStatementHelper (Var _ _ a _ _) op (Expression (FFT (Var _ _ v _ Nothing))) =
    a ++ op ++ "fft(gd, " ++ v ++ ");\n"
  codeStatementHelper _ _ (Expression (FFT _)) =
    error "It is a bug to generate code for a non-var input to fft"
  codeStatementHelper a op (Var _ _ _ _ (Just e)) = codeStatementHelper a op e
  codeStatementHelper (Var _ _ a _ _) op e =
    unlines [setzero,
             "\tfor (int i=1; i<gd.NxNyNzOver2; i++) {",
             "\t\tconst int z = i % gd.NzOver2;",
             "\t\tconst int n = (i-z)/gd.NzOver2;",
             "\t\tconst int y = n % gd.Ny;",
             "\t\tconst int xa = (n-y)/gd.Ny;",
             "\t\tconst RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);",
             "\t\tconst Reciprocal k_i = gd.Lat.toReciprocal(rvec);",
             codes (1 :: Int) e,
             "\t}"]
      where codes n x = case findRepeatedSubExpression x of
              MB (Just (_,x')) -> "\t\tconst complex t"++ show n ++ " = " ++ code x' ++ ";\n" ++
                                  codes (n+1) (substitute x' (s_var ("t"++show n)) x)
              MB Nothing -> "\t\t" ++ a ++ "[i]" ++ op ++ code x ++ ";"
            setzero = case code $ setKequalToZero e of
                      "0.0" -> a ++ "[0]" ++ op ++ "0;"
                      k0code -> unlines ["\t{",
                                         "\t\tconst int i = 0;",
                                         "\t\tconst Reciprocal k_i = Reciprocal(0,0,0);",
                                         "\t\t" ++ a ++ "[0]" ++ op ++ k0code ++ ";",
                                         "\t}"]

  codeStatementHelper _ _ _ = error "Illegal input to codeStatementHelper for kspace"


  newcodeStatementHelper (Var _ _ a _ _) op (Expression (FFT (Var _ _ v _ Nothing))) =
    a ++ op ++ "fft(Nx,Ny,Nz,dV," ++ v ++ ");\n"
  newcodeStatementHelper _ _ (Expression (FFT _)) =
    error "It is a bug to generate newcode for a non-var input to fft"
  newcodeStatementHelper a op (Var _ _ _ _ (Just e)) = newcodeStatementHelper a op e
  newcodeStatementHelper (Var _ _ a _ _) op e =
    unlines [setzero++
             "\tfor (int i=1; i<Nx*Ny*(int(Nz)/2+1); i++) {",
             "\t\tconst int _z = i % (int(Nz)/2+1);",
             "\t\tconst int _n = (i-_z)/(int(Nz)/2+1);",
             "\t\tconst int _y = _n % int(Ny);",
             "\t\tconst int _x = (_n-_y)/int(Ny);",
             "\t\tconst Vector k_i = Vector(" ++ code (xhat `dot` k_i) ++ ", " ++
                                                 code (yhat `dot` k_i) ++ ", " ++
                                                 code (zhat `dot` k_i) ++ ");",
             newcodes (1 :: Int) e,
             "\t}"]
      where k_i = cleanvec $ s_var "_x" / s_var "Nx" .* rlat1 +
                             s_var "_y" / s_var "Ny" .* rlat2 +
                             s_var "_z" / s_var "Nz" .* rlat3
            cleanvec (Vector ea eb ec) = vector (cleanvars ea) (cleanvars eb) (cleanvars ec)
            newcodes n x = case findRepeatedSubExpression x of
              MB (Just (_,x')) ->
                  case break_real_from_imag x' of
                    Expression (Complex r 0) ->
                            "\t\tconst double t"++ show n ++ " = " ++ newcode r ++ ";\n" ++
                            newcodes (n+1) (substitute x' (s_var ("t"++show n)) x)
                    _ -> error "don't hanld e complex well"
              MB Nothing -> "\t\t" ++ a ++ "[i].real()" ++ op ++ newcode (real_part x) ++ ";\n" ++
                            "\t\t" ++ a ++ "[i].imag()" ++ op ++ newcode (imag_part x) ++ ";"
            setzero = case newcode $ setKequalToZero e of
                      "0.0" -> a ++ "[0]" ++ op ++ "0;\n"
                      k0newcode -> unlines ["{",
                                         "\t\tconst int i = 0;",
                                         "\t\t" ++ a ++ "[0]" ++ op ++ k0newcode ++ ";",
                                         "\t\t" ++ a ++ "[1]" ++ op ++ "0;",
                                         "\t}"]

  newcodeStatementHelper _ _ _ = error "Illegal input to newcodeStatementHelper for kspace"
  initialize (Var IsTemp _ x _ Nothing) = "VectorXcd " ++ x ++ "(gd.NxNyNzOver2);"
  initialize _ = error "VectorXcd output(gd.NxNyNzOver2);"
  free (Var IsTemp _ x _ Nothing) = x ++ ".resize(0); // KSpace"
  free _ = error "free error"
  newdeclare _ = "Vector"
  newinitialize (Var _ _ x _ Nothing) = "ComplexVector " ++ x ++ "(Nx*Ny*(int(Nz)/2+1)); // KS"
  newinitialize _ = error "oops newinitialize"
  newfree (Var IsTemp _ x _ Nothing) = x ++ ".free(); // KSpace"
  newfree _ = error "free error"
  toScalar Delta = 1
  toScalar Kx = s_var "_kx"
  toScalar Ky = 0
  toScalar Kz = 0
  toScalar (SetKZeroValue val _) = makeHomogeneous val
  toScalar (FFT e) = makeHomogeneous e
  toScalar (SphericalFourierTransform _ _) = error "need to do spherical transform for toScalar, really need to just evaluate the FT once, and then make this resultingarray[0]"
  toScalar (Complex a _) = a
  mapExpressionHelper' f (FFT e) = fft (f e)
  mapExpressionHelper' f (SphericalFourierTransform s e) = transform s (f e)
  mapExpressionHelper' f (Complex a b) = complex (f a) (f b)
  mapExpressionHelper' f (SetKZeroValue z v) = setkzero (f z) (f v)
  mapExpressionHelper' _ kk = Expression kk
  subAndCountHelper x y (FFT e) = case subAndCount x y e of (e', n) -> (fft e', n)
  subAndCountHelper x y (SphericalFourierTransform s e) =
             case subAndCount x y e of (e', n) -> (transform s e', n)
  subAndCountHelper x y (SetKZeroValue z e) = (setkzero z' e', n1+n2)
        where (z',n1) = subAndCount x y z
              (e',n2) = subAndCount x y e
  subAndCountHelper x y (Complex a b) = (complex a' b', na+nb)
        where (a', na) = subAndCount x y a
              (b', nb) = subAndCount x y b
  subAndCountHelper _ _ Kx = (kx, 0)
  subAndCountHelper _ _ Ky = (ky, 0)
  subAndCountHelper _ _ Kz = (kz, 0)
  subAndCountHelper _ _ Delta = error "unhandled case: Delta"
  searchHelper f (FFT e) = f e
  searchHelper f (SphericalFourierTransform _ e) = f e
  searchHelper f (SetKZeroValue _ e) = f e
  searchHelper f (Complex a b) = mappend (f a) (f b)
  searchHelper _ Kx = mempty
  searchHelper _ Ky = mempty
  searchHelper _ Kz = mempty
  searchHelper _ Delta = mempty
  joinFFThelper (Sum s0 _) = joinup Map.empty $ sum2pairs s0
        where joinup m [] = sum $ map toe $ Map.toList m
                where toe (rs, Right ks) = rs * fft (joinFFTs $ pairs2sum ks)
                      toe (_, Left (_,e)) = e
              joinup m ((f,e):es) =
                case isfft e of
                  Nothing -> toExpression f * e + joinup m es
                  Just (ks,rs) -> joinup (Map.insert rs ks' m) es
                    where ks' = case Map.lookup rs m of
                                Nothing -> Left ([(f,ks)],
                                                 toExpression f * e)
                                Just (Left (ks0,_)) -> Right $ (f, ks) : ks0
                                Just (Right ks0) -> Right $ (f, ks) : ks0
  joinFFThelper e = e
  safeCoerce a _ = case mkExprn a of
                    EK a' -> Just a'
                    _ -> Nothing

mapExprn :: (forall a. Type a => Expression a -> c) -> Exprn -> c
mapExprn f (ES e) = f e
mapExprn f (ER e) = f e
mapExprn f (EK e) = f e

mapExpression :: (Type a, Type b) => (a -> Expression b) -> Expression a -> Expression b
mapExpression f (Var t a b c (Just e)) = Var t a b c $ Just $ mapExpression f e
mapExpression _ (Var tt c v t Nothing) = Var tt c v t Nothing
mapExpression _ (Scalar e) = Scalar e
mapExpression f (Heaviside e) = heaviside (mapExpression f e)
mapExpression f (Cos e) = cos (mapExpression f e)
mapExpression f (Sin e) = sin (mapExpression f e)
mapExpression f (Erf e) = erf (mapExpression f e)
mapExpression f (Exp e) = exp (mapExpression f e)
mapExpression f (Log e) = log (mapExpression f e)
mapExpression f (Abs e) = abs (mapExpression f e)
mapExpression f (Signum e) = signum (mapExpression f e)
mapExpression f (Product p _) = product $ map ff $ product2pairs p
  where ff (e,n) = (mapExpression f e) ** toExpression n
mapExpression f (Sum s _) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, mapExpression f y)
mapExpression f (Expression x) = f x

mapExpression' :: Type a => (forall b. Type b => Expression b -> Expression b) -> Expression a -> Expression a
mapExpression' f (Var IsTemp a b c (Just e)) = f $ Var IsTemp a b c $ Just $ mapExpression' f e
mapExpression' f e@(Var CannotBeFreed _ _ _ (Just _)) = f e
mapExpression' f (Var tt c v t Nothing) = f $ Var tt c v t Nothing
mapExpression' f (Scalar e) = f $ Scalar (mapExpression' f e)
mapExpression' f (Heaviside e) = f $ heaviside (mapExpression' f e)
mapExpression' f (Cos e) = f $ cos (mapExpression' f e)
mapExpression' f (Sin e) = f $ sin (mapExpression' f e)
mapExpression' f (Erf e) = f $ erf (mapExpression' f e)
mapExpression' f (Exp e) = f $ exp (mapExpression' f e)
mapExpression' f (Log e) = f $ log (mapExpression' f e)
mapExpression' f (Abs e) = f $ abs (mapExpression' f e)
mapExpression' f (Signum e) = f $ signum (mapExpression' f e)
mapExpression' f (Product p _) = f $ pairs2product $ map ff $ product2pairs p
  where ff (e,n) = (mapExpression' f e, n)
mapExpression' f (Sum s _) = f $ pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, mapExpression' f y)
mapExpression' f (Expression x) = f $ mapExpressionHelper' (mapExpression' f) x

mapExpressionShortcut :: Type a => (forall b. Type b => Expression b -> Maybe (Expression b))
                         -> Expression a -> Expression a
mapExpressionShortcut f e | Just e' <- f e = e'
mapExpressionShortcut f (Var t a b c (Just e)) = Var t a b c $ Just $ mapExpressionShortcut f e
mapExpressionShortcut _ (Var tt c v t Nothing) = Var tt c v t Nothing
mapExpressionShortcut f (Scalar e) = Scalar (mapExpressionShortcut f e)
mapExpressionShortcut f (Heaviside e) = heaviside (mapExpressionShortcut f e)
mapExpressionShortcut f (Cos e) = cos (mapExpressionShortcut f e)
mapExpressionShortcut f (Sin e) = sin (mapExpressionShortcut f e)
mapExpressionShortcut f (Erf e) = erf (mapExpressionShortcut f e)
mapExpressionShortcut f (Exp e) = exp (mapExpressionShortcut f e)
mapExpressionShortcut f (Log e) = log (mapExpressionShortcut f e)
mapExpressionShortcut f (Abs e) = abs (mapExpressionShortcut f e)
mapExpressionShortcut f (Signum e) = signum (mapExpressionShortcut f e)
mapExpressionShortcut f (Product p _) = pairs2product $ map ff $ product2pairs p
  where ff (e,n) = (mapExpressionShortcut f e, n)
mapExpressionShortcut f (Sum s _) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, mapExpressionShortcut f y)
mapExpressionShortcut f (Expression x) = mapExpressionHelper' (mapExpressionShortcut f) x

-- The first argument, i, is a set of variables.  We search for a
-- subexpression that satisfies the function f (which is the second
-- argument), which itself contains the set of variables i.
searchExpression :: Type a => Set.Set String -> (forall b. Type b => Expression b -> Maybe Exprn)
                    -> Expression a -> Maybe Exprn
searchExpression _ f e | Just c <- f e = Just c
searchExpression i _ e | not $ Set.isSubsetOf i (varSet e) = Nothing
searchExpression _ _ (Var _ _ _ _ Nothing) = Nothing
searchExpression i f v@(Var IsTemp _ _ _ (Just e)) =
  case searchExpression i f e of
    Nothing -> Nothing
    Just e' | mkExprn e == e' -> Just $ mkExprn v
            | otherwise -> Just e'
searchExpression _ _ (Var CannotBeFreed _ _ _ (Just _)) = Nothing
searchExpression i f (Scalar e) = searchExpression i f e
searchExpression i f (Heaviside e) = searchExpression i f e
searchExpression i f (Cos e) = searchExpression i f e
searchExpression i f (Sin e) = searchExpression i f e
searchExpression i f (Erf e) = searchExpression i f e
searchExpression i f (Exp e) = searchExpression i f e
searchExpression i f (Log e) = searchExpression i f e
searchExpression i f (Abs e) = searchExpression i f e
searchExpression i f (Signum e) = searchExpression i f e
searchExpression i f (Product p _) = mconcat $ map (searchExpression i f . fst) $ product2pairs p
searchExpression i f (Sum s _) = mconcat $ map (searchExpression i f . snd) $ sum2pairs s
searchExpression i f (Expression x) = searchHelper (searchExpression i f) x

searchExpressionDepthFirst :: Type a => Set.Set String
                              -> (forall b. Type b => Expression b -> Maybe Exprn)
                              -> Expression a -> Maybe Exprn
searchExpressionDepthFirst i _ e | not $ Set.isSubsetOf i (varSet e) = Nothing
searchExpressionDepthFirst _ f e@(Var _ _ _ _ Nothing) = f e
searchExpressionDepthFirst i f x@(Var IsTemp _ _ _ (Just e)) =
  case searchExpressionDepthFirst i f e of
    Nothing -> f x
    Just e' | mkExprn e == e' -> Just $ mkExprn x
            | otherwise -> Just e'
searchExpressionDepthFirst _ f x@(Var CannotBeFreed _ _ _ (Just _)) = f x
searchExpressionDepthFirst i f x@(Scalar e) = searchExpressionDepthFirst i f e `mor` f x
searchExpressionDepthFirst i f x@(Heaviside e) = searchExpressionDepthFirst i f e `mor` f x
searchExpressionDepthFirst i f x@(Cos e) = searchExpressionDepthFirst i f e `mor` f x
searchExpressionDepthFirst i f x@(Sin e) = searchExpressionDepthFirst i f e `mor` f x
searchExpressionDepthFirst i f x@(Erf e) = searchExpressionDepthFirst i f e `mor` f x
searchExpressionDepthFirst i f x@(Exp e) = searchExpressionDepthFirst i f e `mor` f x
searchExpressionDepthFirst i f x@(Log e) = searchExpressionDepthFirst i f e `mor` f x
searchExpressionDepthFirst i f x@(Abs e) = searchExpressionDepthFirst i f e `mor` f x
searchExpressionDepthFirst i f x@(Signum e) = searchExpressionDepthFirst i f e `mor` f x
searchExpressionDepthFirst i f x@(Product p _) =
  se (map (searchExpressionDepthFirst i f . fst) $ product2pairs p) `mor` f x
  where se [] = Nothing
        se (Just c:_) = Just c
        se (_:cs) = se cs
searchExpressionDepthFirst i f x@(Sum s _) =
  se (map (searchExpressionDepthFirst i f . snd) $ sum2pairs s)  `mor` f x
  where se [] = Nothing
        se (Just c:_) = Just c
        se (_:cs) = se cs
searchExpressionDepthFirst i f x@(Expression e) = searchHelper (searchExpressionDepthFirst i f) e `mor` f x

mor :: Maybe a -> Maybe a -> Maybe a
mor (Just x) _ = Just x
mor _ y = y

cleanvars :: Type a => Expression a -> Expression a
cleanvars = mapExpression' helper
    where helper :: Type a => Expression a -> Expression a
          helper (Var IsTemp b c d (Just e)) | ES _ <- mkExprn e = Var IsTemp b c d (Just e)
                                             | otherwise = e
          helper e = e

isEven :: (Type a) => Exprn -> Expression a -> Double
isEven v e | v == mkExprn e = -1
isEven v (Var _ _ _ _ (Just e)) = isEven v e
isEven _ (Var _ _ _ _ Nothing) = 1
isEven v (Scalar e) = isEven v e
isEven _ (Cos _) = 1
isEven v (Sin e) = isEven v e
isEven v (Erf e) = isEven v e
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
isEven v e = case mkExprn e of
             EK (Expression (FFT r)) -> isEven v r
             EK _ -> 1
             ER (Expression (IFFT ks)) -> isEven v ks
             ES (Expression (Summate x)) -> isEven v x
             _ -> 1 -- Expression _.  Technically, it might be good to recurse into this

-- expand does a few terms in a taylor expansion, and is intended only
-- to be a tool for setting a variable (k, in particular) to zero.
expand :: Type a => Expression a -> Expression a -> Expression a
expand v (Var a b c d (Just e)) = Var a b c d (Just $ expand v e)
expand _ e@(Var _ _ _ _ Nothing) = e
expand _ (Scalar e) = Scalar e
expand v (Heaviside e) = heaviside e'
     where e' = expand v e
expand v (Cos e) = if setZero (mkExprn v) e == 0
                   then 1 - e'**2/2 + e'**4/4/3/2
                   else cos e'
     where e' = expand v e
expand v (Sin e) = if setZero (mkExprn v) e == 0
                   then e' - e'**3/3/2
                   else sin e'
     where e' = expand v e
expand v (Erf e) = if setZero (mkExprn v) e == 0
                   then (2/sqrt pi)*(e' - e'**3/3)
                   else erf e'
     where e' = expand v e
expand v (Exp e) = if setZero (mkExprn v) e == 0
                   then 1 + e' + e'**2/2 + e'**3/3/2
                   else exp e'
     where e' = expand v e
expand v (Log e) = if setZero (mkExprn v) e == 1
                   then - (1 - e') - (1-e')**2/2 - (1-e')**3/3
                   else log e'
     where e' = expand v e
expand v (Abs e) = abs e'
  where e' = expand v e
expand _ (Signum _) = error "ugh signum"
expand v (Sum s _) = pairs2sum $ map ex $ sum2pairs s
  where ex (f,e) = (f, expand v e)
expand v (Product p _) = distribute $ pairs2product $ map ex $ product2pairs p
  where ex (e,n) = (expand v e, n)
expand _ (Expression e) = Expression e

setZero :: Type a => Exprn -> Expression a -> Expression a
setZero v e | v == mkExprn e = 0
            | isEven v e == -1 = 0
setZero v (Var t a b c (Just e)) = case isConstant e' of Nothing -> Var t a b c (Just e')
                                                         Just _ -> e'
  where e' = setZero v e
setZero _ e@(Var _ _ _ _ Nothing) = e
setZero v (Scalar e) = case isConstant $ setZero v e of
                         Just c -> toExpression c
                         Nothing -> Scalar (setZero v e)
setZero v (Heaviside e) = heaviside (setZero v e)
setZero v (Cos e) = cos (setZero v e)
setZero v (Sin e) = sin (setZero v e)
setZero v (Erf e) = erf (setZero v e)
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
         then error ("L'Hopital's rule failure:\n" ++
                     latex n ++ "\n /\n  " ++ latex d ++ "\n\n\n"
                     ++ latex (Product p i) ++ "\n\n\n" ++ latex zn)
         else setZero v (scalarderive v n / scalarderive v d)
  where d = product2denominator p
        n = product $ product2numerator p
        zn = setZero v n
        zd = setZero v d
setZero _ (Sum s i) | Sum s i == 0 = 0
setZero v (Sum s i) =
    case factorize (Sum s i) of
      x | x == Sum s i -> pairs2sum $ map sz $ sum2pairs s
      x -> setZero v x
  where sz (f,x) = (f, setZero v x)
setZero v (Expression x) = zeroHelper v x

instance Code Scalar where
  codePrec _ (Summate r) = showString "integrate(" . codePrec 0 r . showString ")"
  latexPrec _ (Summate r) = showString "\\int " . latexPrec 0 r
instance Type Scalar where
  s_var ("complex(0,1)") = Var CannotBeFreed "complex(0,1)" "complex(0,1)" "i" Nothing
  s_var v@['d',_] = Var CannotBeFreed v v v Nothing -- for differentials
  s_var v = Var CannotBeFreed v v (cleanTex v) Nothing
  s_tex vv tex = Var CannotBeFreed vv vv tex Nothing
  amScalar _ = True
  mkExprn = ES
  derivativeHelper v dds (Summate e) = derive v (scalar dds) e
  scalarderivativeHelper v (Summate e) = summate (scalarderive v e)
  zeroHelper v (Summate e) = summate (setZero v e)
  codeStatementHelper a " = " (Var _ _ _ _ (Just e)) = codeStatementHelper a " = " e
  codeStatementHelper a " = " (Expression (Summate e)) =
    code a ++ " = 0;\n\tfor (int i=0; i<gd.NxNyNz; i++) {\n\t\t" ++
    code a ++ " += " ++ code e ++
    ";\n\t}\n"
  codeStatementHelper _ op (Expression (Summate _)) = error ("Haven't implemented "++op++" for integrate...")
  codeStatementHelper a op e = code a ++ op ++ code e ++ ";"

  newcodeStatementHelper a " = " (Var _ _ _ _ (Just e)) = newcodeStatementHelper a " = " e
  newcodeStatementHelper a " = " (Expression (Summate e)) =
    newcode a ++ " = 0;\n\tfor (int i=0; i<Nx*Ny*Nz; i++) {\n\t\t" ++
    newcode a ++ " += " ++ newcode e ++
    ";\n\t}\n"
  newcodeStatementHelper _ op (Expression (Summate _)) = error ("Haven't implemented "++op++" for integrate...")
  newcodeStatementHelper a op e = newcode a ++ op ++ newcode e ++ ";"
  initialize (Var _ _ x _ Nothing) = "double " ++ x ++ " = 0;\n"
  initialize v = error ("bug in initialize Scalar: "++show v)
  newdeclare _ = "double"
  toScalar (Summate r) = makeHomogeneous (r/dV)
  fromScalar = id
  mapExpressionHelper' f (Summate e) = summate (f e)
  subAndCountHelper x y (Summate e) = case subAndCount x y e of (e', n) -> (summate e', n)
  searchHelper f (Summate e) = f e
  safeCoerce a _ = case mkExprn a of
                    ES a' -> Just a'
                    _ -> Nothing

scalar :: Type a => Expression Scalar -> Expression a
scalar (Scalar e) = scalar e
scalar e | Just c <- isConstant e = toExpression c
scalar e = fromScalar e

cleanTex :: String -> String
cleanTex [a] = [a]
cleanTex ('_':v) = cleanTex v
cleanTex (a:v) = a : '_' : '{' : v ++ "}"
cleanTex "" = error "cleanTex on empty string"

r_var :: String -> Expression RealSpace
r_var v | take 5 v == "rtemp" = Var IsTemp (v++"[i]") ('r':v) v Nothing
r_var v | take 4 v == "temp" = Var IsTemp ("r"++v++"[i]") ('r':v) v Nothing
r_var v = Var CannotBeFreed (v++"[i]") v (cleanTex v) Nothing

k_var :: String -> Expression KSpace
k_var v | take 5 v == "ktemp" = Var IsTemp (v++"[i]") ('r':v) v Nothing
k_var v | take 4 v == "temp" = Var IsTemp ("k"++v++"[i]") ('r':v) v Nothing
k_var v = Var CannotBeFreed (v++"[i]") v (cleanTex v) Nothing

t_var :: String -> Vector Scalar
t_var v = Vector (s_var $ v++"x") (s_var $ v++"y") (s_var $ v++"z")

imaginary :: Expression KSpace
imaginary = complex 0 1

infix 4 ===, `nameVector`

(===) :: Type a => String -> Expression a -> Expression a
--_ === e = e
v@(a:r@(_:_)) === e = Var IsTemp c v ltx (Just e)
  where ltx = a : "_{"++r++"}"
        c = if amScalar e then v else v ++ "[i]"
v === e = Var IsTemp c v v (Just e)
  where c = if amScalar e then v else v ++ "[i]"

var :: Type a => String -> String -> Expression a -> Expression a
var v ltx e = Var IsTemp c v ltx (Just e)
  where c = if amScalar e then v else v ++ "[i]"

protect :: Type a => String -> String -> Expression a -> Expression a
protect v ltx e = Var CannotBeFreed c v ltx (Just e)
  where c = if amScalar e then v else v ++ "[i]"

vprotect :: Type a => String -> (String -> String) -> Vector a -> Vector a
vprotect v ltx (Vector x y z) = Vector (protect (v++"x") (ltx "x") x) 
                                       (protect (v++"y") (ltx "y") y) 
                                       (protect (v++"z") (ltx "z") z)

rmag :: Expression RealSpace
rmag = sqrt (rx**2 + ry**2 + rz**2)

rvec :: Vector RealSpace
rvec = Vector rx ry rz

rx :: Expression RealSpace
rx = Expression Rx

ry :: Expression RealSpace
ry = Expression Ry

rz :: Expression RealSpace
rz = Expression Rz

kx :: Expression KSpace
kx = Expression Kx

ky :: Expression KSpace
ky = Expression Ky

kz :: Expression KSpace
kz = Expression Kz

kvec :: Vector KSpace
kvec = Vector kx ky kz

k :: Expression KSpace
k = sqrt ksqr

ksqr :: Expression KSpace
ksqr = kx**2 + ky**2 + kz**2

setkzero :: Expression KSpace -> Expression KSpace -> Expression KSpace
setkzero zeroval otherval = Expression $ SetKZeroValue zeroval otherval

integrate :: Type a => Expression RealSpace -> Expression a
integrate f = summate (f*dV)

summate :: Type a => Expression RealSpace -> Expression a
summate (Sum s _) = pairs2sum $ map i $ sum2pairs s
  where i (f, e) = (f, summate e)
summate x = fromScalar $ Expression $ Summate x

transform :: Symmetry -> Expression Scalar -> Expression KSpace
transform s e = if e == 0
                then 0
                else Expression (SphericalFourierTransform s e)

fft :: Expression RealSpace -> Expression KSpace
fft (Scalar e) = Scalar e * Expression Delta
fft r | r == 0 = 0
      | otherwise = Expression (FFT r)

ifft :: Expression KSpace -> Expression RealSpace
ifft ke | ke == 0 = 0
        | otherwise = Expression (IFFT ke)

vfft :: Vector RealSpace -> Vector KSpace
vfft (Vector x y z) = Vector (fft x) (fft y) (fft z)

vifft :: Vector KSpace -> Vector RealSpace
vifft (Vector x y z) = Vector (ifft x) (ifft y) (ifft z)

complex :: Expression Scalar -> Expression Scalar -> Expression KSpace
complex a b | b == 0 = scalar a
            | otherwise = Expression (Complex a b)

break_real_from_imag :: Expression KSpace -> Expression KSpace
break_real_from_imag = brfi
  where brfi (Product p _) = handle 1 0 $ product2pairs p
          where handle re im [] = Expression (Complex re im)
                handle re im ((x,n):es)
                   | imx == 0 = handle (re* rex ** toExpression n) (im * rex ** toExpression n) es
                   | n == 1 = handle (re*rex - im*imx) (re*imx + im*rex) es
                   | otherwise = error ("haven't finished break_real_from_imag " ++ code x ++ " ** " ++ show n
                                       ++ "\nsee also " ++ code rex ++ " and imaginary " ++ code imx)
                     where Expression (Complex rex imx) = brfi x
        brfi (Sum s _) = handle [] [] $ sum2pairs s
           where handle re im ((f,x):xs) = handle ((f,r):re) ((f,i):im) xs
                   where Expression (Complex r i) = brfi x
                 handle re im [] = Expression (Complex (pairs2sum re) (pairs2sum im))
        brfi (Erf e) = case brfi e of
                         Expression (Complex r 0) -> Expression (Complex (erf r) 0)
                         _ -> error "cerazinesss in erf"
        brfi (Exp e) = case brfi e of
                         Expression (Complex r 0) -> Expression (Complex (exp r) 0)
                         Expression (Complex r i) -> Expression (Complex (exp r * cos i) (exp r * sin i))
                         _ -> error "ceraziness"
        brfi (Sin e) = case brfi e of
                         Expression (Complex r 0) -> Expression (Complex (sin r) 0)
                         _ -> error "ceraziness in sin"
        brfi (Cos e) = case brfi e of
                         Expression (Complex r 0) -> Expression (Complex (cos r) 0)
                         _ -> error "ceraziness in cos"
        brfi (Var t _ b tex Nothing) = Expression $ Complex (Var t (b++"[i].real()") b tex Nothing)
                                                            (Var t (b++"[i].imag()") b tex Nothing)
        brfi (Var a b c d (Just e))
          | i == 0 && r /= 0 = Expression $ Complex (Var a b c d (Just r)) 0
          | otherwise = Expression $ Complex r i
          where Expression (Complex r i) = brfi e
        brfi (Expression Kx) = Expression $ Complex (s_var "k_i[0]") 0
        brfi (Expression Ky) = Expression $ Complex (s_var "k_i[1]") 0
        brfi (Expression Kz) = Expression $ Complex (s_var "k_i[2]") 0
        brfi (Scalar e) = Expression $ Complex e 0
        brfi (Expression (Complex a b)) = Expression (Complex a b)
        brfi e = error ("brfi doesn't handle " ++ show e)

real_part :: Expression KSpace -> Expression Scalar
real_part x = case break_real_from_imag x of
                Expression (Complex r _) -> r
                _ -> error "oopsies"

imag_part :: Expression KSpace -> Expression Scalar
imag_part x = case break_real_from_imag x of
                Expression (Complex _ i) -> i
                _ -> error "oopsies"
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

-- An expression statement holds a mathematical expression, which is
-- be any one of several different types: RealSpace (as in n(\vec{r})),
-- KSpace (as in w(\vec{k})), or Scalar (as in an ordinary scalar value
-- like k_BT.
data Expression a = Scalar (Expression Scalar) |
                    Var IsTemp String String String (Maybe (Expression a)) | -- A variable with a possible value
                    Expression a |
                    Heaviside (Expression a) |
                    Cos (Expression a) |
                    Sin (Expression a) |
                    Exp (Expression a) |
                    Log (Expression a) |
                    Erf (Expression a) |
                    Abs (Expression a) |
                    Signum (Expression a) |
                    Product (Map.Map (Expression a) Double) (Set.Set String) |
                    Sum (Map.Map (Expression a) Double) (Set.Set String)
              deriving (Eq, Ord, Show)

-- An "E" type is *any* type of expression.  This is useful when we
-- want to specify subexpressions, for instance, when we might not
-- want to care in advance which sort of subexpression it is.  Also
-- useful for comparing two expressions that might not be the same
-- type.
data Exprn = EK (Expression KSpace) | ER (Expression RealSpace) | ES (Expression Scalar)
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
        fl a ((f,x):xs) = fl (Map.alter upd x a) xs
          where upd Nothing = Just f
                upd (Just f') = if f + f' == 0
                                then Nothing
                                else Just (f+f')

map2sum :: Type a => Map.Map (Expression a) Double -> Expression a
map2sum s | Map.size s == 1 =
  case sum2pairs s of [(1,e)] -> e
                      _ -> Sum s (Set.unions $ map varSet $ Map.keys s)
map2sum s = Sum s (Set.unions $ map varSet $ Map.keys s)

product2pairs :: Type a => Map.Map (Expression a) Double -> [(Expression a, Double)]
product2pairs s = Map.assocs s

-- map2sum converts a Map to a product.  It handles the case of a
-- singleton map properly.  It also appropriately combines any
-- constant values in the product, and eliminates any terms with a
-- zero power (which are thus one).  Unline pairs2product, map2product
-- *doesn't* need to check for a given expression occurring more than
-- once, which makes it a bit simpler.
map2product :: Type a => Map.Map (Expression a) Double -> Expression a
map2product p | Map.size p == 1 =
  case product2pairs p of [(e,1)] -> e
                          _ -> Product p (Set.unions $ map varSet $ Map.keys p)
map2product p = helper 1 (Map.empty) $ product2pairs p
  where helper 1 a [] = Product a (vl a)
        helper f a [] = Sum (Map.singleton (Product a i) f) i
          where i = vl a
        helper f a ((Sum x _,n):xs) | [(f',x')] <- sum2pairs x = helper (f*f'**n) a ((x',n):xs)
        helper f a ((x,n):xs)
          | n == 0 = helper f a xs
          | x == 1 = helper f a xs
          | otherwise = helper f (Map.insert x n a) xs
        vl = Set.unions . map varSet . Map.keys

-- pairs2product combines a list of expressions and their powers into
-- a product expression.  It calls map2product internally, so it
-- automatically handles anything map2product can handle.  In
-- addition, it handles situations where the same expression occurs
-- several times, such as x**2/y/x.
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

instance Code Exprn where
  codePrec p = mapExprn (codePrec p)
  latexPrec p = mapExprn (latexPrec p)
instance (Type a, Code a) => Code (Expression a) where
  codePrec _ (Var _ c _ _ Nothing) = showString c
  codePrec p (Var _ _ _ _ (Just e)) = codePrec p e
  codePrec p (Scalar x) = codePrec p x
  codePrec p (Expression x) = codePrec p x
  codePrec _ (Heaviside x) = showString "heaviside(" . codePrec 0 x . showString ")"
  codePrec _ (Cos x) = showString "cos(" . codePrec 0 x . showString ")"
  codePrec _ (Sin x) = showString "sin(" . codePrec 0 x . showString ")"
  codePrec _ (Exp x) = showString "exp(" . codePrec 0 x . showString ")"
  codePrec _ (Erf x) = showString "erf(" . codePrec 0 x . showString ")"
  codePrec _ (Log x) = showString "log(" . codePrec 0 x . showString ")"
  codePrec _ (Abs x) = showString "fabs(" . codePrec 0 x . showString ")"
  codePrec _ (Signum _) = undefined
  codePrec _ (Product p i) | Product p i == 1 = showString "1.0"
  codePrec pree (Product p _) = showParen (pree > 7) $
                           if den == 1
                           then codesimple num
                           else codesimple num . showString "/" . codePrec 8 den
    where codesimple [] = showString "1.0"
          codesimple [(a,n)] = codee a n
          codesimple [(a,n),(b,m)] = codee a n . showString "*" . codee b m
          codesimple ((a,n):es) = codee a n . showString "*" . codesimple es
          num = product2numerator_pairs p
          den = product2denominator p
          codee _ 0 = showString "1.0" -- this shouldn't happen...
          codee _ n | n < 0 = error "shouldn't have negative power here"
          codee x 1 = codePrec 7 x
          codee x 0.5 = showString "sqrt(" . codePrec 0 x . showString ")"
          codee x nn
            | fromInteger n2 == 2*nn && odd n2 = codee x 0.5 . showString "*" . codee x (nn-0.5)
            | fromInteger n == nn && odd n = codee x 1 . showString "*" . codee x (nn-1)
            | fromInteger n == nn =
              showParen (nn/2>1) (codee x (nn / 2)) . showString "*" . showParen (nn/2>1) (codee x (nn / 2))
            where n2 = floor (2*nn)
                  n = floor nn
          codee x n = showString "pow(" . codePrec 0 x . showString (", " ++ show n ++ ")")
  codePrec _ (Sum s i) | Sum s i == 0 = showString "0.0"
  codePrec p (Sum s _) = showParen (p > 6) (showString me)
    where me = foldl addup "" $ sum2pairs s
          addup "" (1,e) = codePrec 6 e ""
          addup "" (f,e) = if e == 1
                           then show f
                           else show f ++ "*" ++ codePrec 7 e ""
          addup rest (1,e) = codePrec 6 e (showString " + " $ rest)
          addup rest (f,e) = show f ++ "*" ++ codePrec 7 e (showString " + " $ rest)
  latexPrec p (Var _ _ "" "" (Just e)) = latexPrec p e
  latexPrec _ (Var _ _ c "" _) = showString c
  latexPrec _ (Var _ _ _ t _) = showString t
  latexPrec _ x | Just xx <- isConstant x = showString (latexDouble xx)
  latexPrec p (Scalar x) = latexPrec p x
  latexPrec p (Expression x) = latexPrec p x
  latexPrec _ (Heaviside x) = showString "\\Theta(" . latexPrec 0 x . showString ")"
  latexPrec _ (Cos x) = showString "\\cos(" . latexPrec 0 x . showString ")"
  latexPrec _ (Sin x) = showString "\\sin(" . latexPrec 0 x . showString ")"
  latexPrec _ (Exp x) = showString "\\exp\\left(" . latexPrec 0 x . showString "\\right)"
  latexPrec _ (Erf x) = showString "\\textrm{erf}\\left(" . latexPrec 0 x . showString "\\right)"
  latexPrec _ (Log x) = showString "\\log(" . latexPrec 0 x . showString ")"
  latexPrec _ (Abs x) = showString "\\left|" . latexPrec 0 x . showString "\\right|"
  latexPrec _ (Signum _) = undefined
  latexPrec p (Product x _) | Map.size x == 1 && product2denominator x == 1 =
    case product2pairs x of
      [(_,0)] -> showString "1" -- this shouldn't happen...
      [(_, n)] | n < 0 -> error "shouldn't have negative power here"
      [(e, 1)] ->   latexPrec p e
      [(e, 0.5)] -> showString "\\sqrt{" . latexPrec 0 e . showString "}"
      [(e, n)] | floor n == (ceiling n :: Int) && n < 10 -> latexPrec 8 e . showString ("^" ++ latexDouble n)
               | otherwise -> latexPrec 8 e . showString ("^{" ++ latexDouble n ++ "}")
      _ -> error "This really cannot happen."
  latexPrec pree (Product p _) | product2denominator p == 1 = latexParen (pree > 7) $ ltexsimple $ product2numerator p
    where ltexsimple [] = showString "1"
          ltexsimple [a] = latexPrec 7 a
          ltexsimple [a,b] = latexPrec 7 a . showString " " . latexPrec 7 b
          ltexsimple (a:es) = latexPrec 7 a . showString " " . ltexsimple es
  latexPrec pree (Product p _) = latexParen (pree > 7) $
              showString "\\frac{" . latexPrec 0 num . showString "}{" .
                                     latexPrec 0 den . showString "}"
    where num = product $ product2numerator p
          den = product2denominator p
  latexPrec p (Sum s _) = latexParen (p > 6) (showString me)
    where me = foldl addup "" $ sum2pairs s
          addup "" (1,e) = latexPrec 6 e ""
          addup "" (f,e) | f < 0 = "-" ++ addup "" (-f, e)
          addup "" (f,e) = if e == 1
                           then latexDouble f
                           else latexDouble f ++ " " ++ latexPrec 6 e ""
          addup ('-':rest) (1,e) = latexPrec 6 e (" - " ++ rest)
          addup ('-':rest) (-1,e) = "-" ++ latexPrec 6 e (" - " ++ rest)
          addup ('-':rest) (f,e) = latexDouble f ++ " " ++ latexPrec 7 e (showString " - " $ rest)
          addup rest (-1,e) = "-" ++ latexPrec 6 e (" + " ++ rest)
          addup rest (1,e) = latexPrec 6 e (" + " ++ rest)
          addup rest (f,e) = latexDouble f ++ " " ++ latexPrec 7 e (showString " + " $ rest)

latexParen :: Bool -> ShowS -> ShowS
latexParen False x = x
latexParen True x = showString "\\left(" . x . showString "\\right)"

class Code a  where
    codePrec  :: Int -> a -> ShowS
    codePrec _ x s = code x ++ s
    code      :: a -> String
    code x = codePrec 0 x ""
    newcodePrec  :: Int -> a -> ShowS
    newcodePrec = codePrec
    newcode      :: a -> String
    newcode x = newcodePrec 0 x ""
    latexPrec :: Int -> a -> ShowS
    latexPrec _ x s = latex x ++ s
    latex     :: a -> String
    latex x = latexPrec 0 x ""

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

class Monoid m where
  mappend :: m -> m -> m
  mempty :: m
instance Monoid (Maybe a) where
  mappend (Just x) _ = Just x
  mappend Nothing y = y
  mempty = Nothing
instance Ord a => Monoid (Set.Set a) where
  mappend x y = Set.union x y
  mempty = Set.empty
mconcat :: Monoid m => [m] -> m
mconcat (x:xs) = mappend x (mconcat xs)
mconcat [] = mempty

class (Ord a, Show a, Code a) => Type a where 
  amScalar :: Expression a -> Bool
  amScalar _ = False
  amRealSpace :: Expression a -> Bool
  amRealSpace _ = False
  amKSpace :: Expression a -> Bool
  amKSpace _ = False
  mkExprn :: Expression a -> Exprn
  s_var :: String -> Expression a
  s_var = Scalar . s_var
  s_tex :: String -> String -> Expression a
  s_tex v t = Scalar (s_tex v t)
  derivativeHelper :: Type b => Expression b -> Expression a -> a -> Expression b
  scalarderivativeHelper :: Exprn -> a -> Expression a
  zeroHelper :: Exprn -> a -> Expression a
  codeStatementHelper :: Expression a -> String -> Expression a -> String
  newcodeStatementHelper :: Expression a -> String -> Expression a -> String
  initialize :: Expression a -> String
  newinitialize :: Expression a -> String
  newinitialize e@(Var _ _ x _ _) = newdeclare e ++ " " ++ x ++ ";"
  newinitialize _ = error "bad newinitialize"
  newdeclare :: Expression a -> String
  free :: Expression a -> String
  free x = error ("free nothing " ++ show x)
  newfree :: Expression a -> String
  newfree x = error ("free nothing " ++ show x)
  toScalar :: a -> Expression Scalar
  fromScalar :: Expression Scalar -> Expression a
  fromScalar = Scalar
  mapExpressionHelper' :: (forall b. Type b => Expression b -> Expression b) -> a -> Expression a
  joinFFThelper :: Expression a -> Expression a
  joinFFThelper = id
  safeCoerce :: Type b => Expression b -> Expression a -> Maybe (Expression a)
  subAndCountHelper :: Type b => Expression b -> Expression b -> a -> (Expression a, Int)
  searchHelper :: Monoid c => (forall b. Type b => Expression b -> c) -> a -> c

initializeE :: Exprn -> String
initializeE (ES e) = initialize e
initializeE (EK e) = initialize e
initializeE (ER e) = initialize e
freeE :: Exprn -> String
freeE (ES e) = free e
freeE (ER e) = free e
freeE (EK e) = free e

newinitializeE :: Exprn -> String
newinitializeE (ES e) = newinitialize e
newinitializeE (EK e) = newinitialize e
newinitializeE (ER e) = newinitialize e

newdeclareE :: Exprn -> String
newdeclareE (ES e) = newdeclare e
newdeclareE (EK e) = newdeclare e
newdeclareE (ER e) = newdeclare e

nameE :: Exprn -> String
nameE (ES (Var _ _ v _ Nothing)) = v
nameE (EK (Var _ _ v _ Nothing)) = v
nameE (ER (Var _ _ v _ Nothing)) = v
nameE e = show e

newfreeE :: Exprn -> String
newfreeE (ES e) = newfree e
newfreeE (ER e) = newfree e
newfreeE (EK e) = newfree e

codeStatementE :: Exprn -> String -> Exprn -> String
codeStatementE (ES a) op (ES b) = codeStatementHelper a op b
codeStatementE (EK a) op (EK b) = codeStatementHelper a op b
codeStatementE (ER a) op (ER b) = codeStatementHelper a op b
codeStatementE _ _ _ = error "bug revealed by codeStatementE"

newcodeStatementE :: Exprn -> String -> Exprn -> String
newcodeStatementE (ES a) op (ES b) = newcodeStatementHelper a op b
newcodeStatementE (EK a) op (EK b) = newcodeStatementHelper a op b
newcodeStatementE (ER a) op (ER b) = newcodeStatementHelper a op b
newcodeStatementE _ _ _ = error "bug revealed by newcodeStatementE"

makeHomogeneous :: Type a => Expression a -> Expression Scalar
makeHomogeneous ee =
  scalarScalar $ setZero (ES (s_var "_kx")) $ expand (s_var "_kx") $ mapExpression toScalar ee
  where scalarScalar :: Expression Scalar -> Expression Scalar
        scalarScalar (Var t a b c (Just e)) = Var t a b c (Just $ scalarScalar e)
        scalarScalar (Var _ _ c l Nothing) = Var CannotBeFreed c c l Nothing
        scalarScalar (Scalar s) = s
        scalarScalar (Sum x _) = pairs2sum $ map f $ sum2pairs x
          where f (a,b) = (a, scalarScalar b)
        scalarScalar (Product x _) = pairs2product $ map f $ product2pairs x
          where f (a,b) = (scalarScalar a, b)
        scalarScalar (Expression e) = Expression e -- FIXME
        scalarScalar (Heaviside x) = heaviside (scalarScalar x)
        scalarScalar (Cos x) = cos (scalarScalar x)
        scalarScalar (Sin x) = sin (scalarScalar x)
        scalarScalar (Exp x) = exp (scalarScalar x)
        scalarScalar (Erf x) = erf (scalarScalar x)
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

erf :: Type a => Expression a -> Expression a
erf x = case x of 0 -> 0
                  _ -> Erf x

heaviside :: Type a => Expression a -> Expression a
heaviside x = case isConstant x of
              Just n -> if n >= 0 then 1 else 0
              Nothing -> Heaviside x

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

countVars :: [Exprn] -> Int
countVars s = Set.size $ mconcat $ map (mapExprn varSet) s

varSet :: Type a => Expression a -> Set.Set String
varSet e@(Expression _) = case mkExprn e of
                          EK (Expression (FFT e')) -> varSet e'
                          EK (Expression (SetKZeroValue _ e')) -> varSet e'
                          ER (Expression (IFFT e')) -> varSet e'
                          ES (Expression (Summate e')) -> varSet e'
                          _ -> Set.empty
varSet (Var _ _ _ _ (Just e)) = varSet e
varSet (Var IsTemp _ c _ Nothing) = Set.singleton c
varSet (Var CannotBeFreed _ _ _ Nothing) = Set.empty
varSet (Heaviside e) = varSet e
varSet (Sin e) = varSet e
varSet (Cos e) = varSet e
varSet (Log e) = varSet e
varSet (Exp e) = varSet e
varSet (Erf e) = varSet e
varSet (Abs e) = varSet e
varSet (Signum e) = varSet e
varSet (Sum _ i) = i
varSet (Product _ i) = i
varSet (Scalar _) = Set.empty

hasFFT :: Type a => Expression a -> Bool
hasFFT e@(Expression _) = case mkExprn e of
  EK (Expression (FFT _)) -> True
  EK (Expression (SetKZeroValue z e')) -> hasFFT e' || hasFFT z
  ER (Expression (IFFT _)) -> True
  ES (Expression (Summate _)) -> True -- a bit weird... rename this function?
  _ -> False
hasFFT (Var _ _ _ _ (Just e)) = hasFFT e
hasFFT (Sum s _) = or $ map (hasFFT . snd) (sum2pairs s)
hasFFT (Product p _) = or $ map (hasFFT . fst) (product2pairs p)
hasFFT (Heaviside e) = hasFFT e
hasFFT (Sin e) = hasFFT e
hasFFT (Cos e) = hasFFT e
hasFFT (Log e) = hasFFT e
hasFFT (Exp e) = hasFFT e
hasFFT (Erf e) = hasFFT e
hasFFT (Abs e) = hasFFT e
hasFFT (Signum e) = hasFFT e
hasFFT (Var _ _ _ _ Nothing) = False
hasFFT (Scalar e) = hasFFT e

hasActualFFT :: Type a => Expression a -> Bool
hasActualFFT e@(Expression _) = case mkExprn e of
  EK (Expression (FFT _)) -> True
  EK (Expression (SetKZeroValue z e')) -> hasActualFFT e' || hasActualFFT z
  ER (Expression (IFFT _)) -> True
  ES (Expression (Summate e')) -> hasActualFFT e'
  _ -> False
hasActualFFT (Var _ _ _ _ (Just e)) = hasActualFFT e
hasActualFFT (Sum s _) = or $ map (hasActualFFT . snd) (sum2pairs s)
hasActualFFT (Product p _) = or $ map (hasActualFFT . fst) (product2pairs p)
hasActualFFT (Heaviside e) = hasActualFFT e
hasActualFFT (Sin e) = hasActualFFT e
hasActualFFT (Cos e) = hasActualFFT e
hasActualFFT (Log e) = hasActualFFT e
hasActualFFT (Exp e) = hasActualFFT e
hasActualFFT (Erf e) = hasActualFFT e
hasActualFFT (Abs e) = hasActualFFT e
hasActualFFT (Signum e) = hasActualFFT e
hasActualFFT (Var _ _ _ _ Nothing) = False
hasActualFFT (Scalar e) = hasActualFFT e

isfft :: Expression KSpace -> Maybe (Expression RealSpace, Expression KSpace)
isfft (Expression (FFT e)) = Just (e,1)
isfft (Product p _) = tofft 1 1 Nothing $ product2pairs p
  where tofft _ _ Nothing [] = Nothing
        tofft sc ks (Just e) [] = Just (sc*e, ks)
        tofft sc ks Nothing ((Expression (FFT e),1):xs) = tofft sc ks (Just e) xs
        tofft sc ks myfft ((Scalar s,n):xs) = tofft (sc * Scalar s ** toExpression n) ks myfft xs
        tofft sc ks myfft ((kk,n):xs) = tofft sc (ks * kk ** toExpression n) myfft xs
isfft (Var _ _ _ _ (Just e)) = isfft e
isfft _ = Nothing

isifft :: Expression RealSpace -> Maybe (Expression KSpace, Expression RealSpace)
isifft (Expression (IFFT e)) = Just (e,1)
isifft (Product p _) = tofft 1 1 Nothing $ product2pairs p
  where tofft _ _ Nothing [] = Nothing
        tofft sc rs (Just e) [] = Just (sc*e, rs)
        tofft sc rs Nothing ((Expression (IFFT e),1):xs) = tofft sc rs (Just e) xs
        tofft sc rs myfft ((Scalar s,n):xs) = tofft (sc * Scalar s ** toExpression n) rs myfft xs
        tofft sc rs myfft ((r,n):xs) = tofft sc (rs * r ** toExpression n) myfft xs
isifft (Var _ _ _ _ (Just e)) = isifft e
isifft _ = Nothing

joinFFTs :: Type a => Expression a -> Expression a
joinFFTs = mapExpression' joinFFThelper

factorOut :: Type a => Expression a -> Expression a -> Maybe Double
factorOut e e' | e == e' = Just 1
factorOut e (Sum s i) | Set.isSubsetOf (varSet e) i = factorOutSumPairs e (sum2pairs s)
factorOut e (Product p _) = Map.lookup e p
factorOut _ _ = Nothing

factorOutSumPairs :: Type a => Expression a -> [(Double, Expression a)] -> Maybe Double
factorOutSumPairs _ [] = Nothing
factorOutSumPairs e pairs = findpower 0 pairs
  where findpower smallestpower ((_,x):xs) = case factorOut e x of
                                             Just n | n*smallestpower < 0 -> Nothing
                                                    | smallestpower == 0 -> findpower n xs
                                                    | n > 0 -> findpower (min smallestpower n) xs
                                                    | n < 0 -> findpower (max smallestpower n) xs
                                             _ -> Nothing
        findpower smallestpower [] = Just smallestpower

factorizeSumPairs :: Type a => [(Double, Expression a)] -> Expression a
factorizeSumPairs [] = 0
factorizeSumPairs p@([_]) = pairs2sum p
factorizeSumPairs xs@((_,x):_) = fsp [] (Set.toList $ allFactors x) xs
  where fsp prefac (f:fs) pairs = case factorOutSumPairs f pairs of
                                     Nothing -> fsp prefac fs pairs
                                     Just n -> fsp ((f,n):prefac) fs (map divbyme pairs)
                                       where divbyme (a,e) = (a, e / fton)
                                             fton = f ** toExpression n
        fsp prefac [] pairs = pairs2product prefac * pairs2sum pairs

allFactors :: Type a => Expression a -> Set.Set (Expression a)
allFactors (Product p _) = Map.keysSet p
allFactors (Sum s _) | [(_,e)] <- sum2pairs s = allFactors e
allFactors e = Set.singleton e

factorize :: Type a => Expression a -> Expression a
--factorize = mapExpressionShortcut factorizeHelper

-- The following factorizes more thoroughly, which leads to quicker
-- code generation, but ends up needing more memory when running
-- (although it also runs faster).

factorize = mapExpression' helper
     where helper e | Just e' <- factorizeHelper e = e'
                    | otherwise = e

factorizeHelper :: Type a => Expression a -> Maybe (Expression a)
factorizeHelper (Sum s _) = Just $ fac (Set.toList $ Set.unions $ map (allFactors . snd) $ sum2pairs s) $
                            map toe $ sum2pairs s
  where toe (f,e) = toExpression f * e
        fac _ [] = 0
        fac [] pairs = sum pairs
        fac (f:fs) pairs = collect [] [] [] pairs
          where collect pos none neg (x:xs) = case factorOut f x of
                                                Just n | n <= -1 -> collect pos none (x:neg) xs
                                                       | n >= 1 -> collect (x:pos) none neg xs
                                                _ -> collect pos (x:none) neg xs
                collect pos none neg []
                  | length pos <= 1 && length neg <= 1 =
                    fac fs (pos ++ none ++ neg)
                  | length pos <= 1 = fac fs (pos ++ none) +
                                      (fac (f:fs) (map (*f) neg))/f
                  | length neg <= 1 = fac fs (neg ++ none) +
                                      f*(fac (f:fs) (map (/f) pos))
                  | otherwise = f * fac (f:fs) (map (/f) pos) +
                                fac fs none +
                                (fac (f:fs) (map (*f) neg))/f
factorizeHelper _ = Nothing


-- distribute is the complement of factorize, but it currently is only
-- used as a helper in creating a Taylor expansion in "expand".
distribute :: Type a => Expression a -> Expression a
distribute = mapExpression' helper
  where helper (Product p _) = sum $ dist $ product2pairs p
          where dist :: Type a => [(Expression a, Double)] -> [Expression a]
                dist [] = [1]
                dist ((Sum s i,n):es) | n > 1 = dist $ (Sum s i, 1):(Sum s i, n-1):es
                dist ((Sum s _,1):es) = concatMap (\e -> map (e *) es') $ sum2es s
                  where es' = dist es
                dist ((e,n):es) = map ((e ** toExpression n) *) $ dist es
        helper e = e
        sum2es :: Type a => Map.Map (Expression a) Double -> [Expression a]
        sum2es = map (\(f,x) -> toExpression f * x) . sum2pairs

hasExpressionInFFT :: (Type a, Type b) => Expression b -> Expression a -> Bool
hasExpressionInFFT v e | not (hasexpression v e) = False
hasExpressionInFFT v (Expression e) = case mkExprn (Expression e) of
                                      EK (Expression (FFT e')) -> hasexpression v e'
                                      EK (Expression (SetKZeroValue _ e')) -> hasExpressionInFFT v e'
                                      ER (Expression (IFFT e')) -> hasexpression v e'
                                      ES (Expression (Summate e')) -> hasExpressionInFFT v e'
                                      EK (Expression Kx) -> False
                                      EK (Expression Ky) -> False
                                      EK (Expression Kz) -> False
                                      _ -> error "inexhaustive pattern in hasExpressionInFFT"
hasExpressionInFFT v (Var _ _ _ _ (Just e)) = hasExpressionInFFT v e
hasExpressionInFFT v (Sum s _) = or $ map (hasExpressionInFFT v . snd) (sum2pairs s)
hasExpressionInFFT v (Product p _) = or $ map (hasExpressionInFFT v . fst) (product2pairs p)
hasExpressionInFFT v (Heaviside e) = hasExpressionInFFT v e
hasExpressionInFFT v (Sin e) = hasExpressionInFFT v e
hasExpressionInFFT v (Cos e) = hasExpressionInFFT v e
hasExpressionInFFT v (Log e) = hasExpressionInFFT v e
hasExpressionInFFT v (Exp e) = hasExpressionInFFT v e
hasExpressionInFFT v (Erf e) = hasExpressionInFFT v e
hasExpressionInFFT v (Abs e) = hasExpressionInFFT v e
hasExpressionInFFT v (Signum e) = hasExpressionInFFT v e
hasExpressionInFFT v (Scalar e) = hasExpressionInFFT v e
hasExpressionInFFT _ (Var _ _ _ _ Nothing) = False

-- scalarderive gives a derivative of the same type as the original,
-- and will always either be a derivative with respect to a scalar, or
-- with respect to kx.
scalarderive :: Type a => Exprn -> Expression a -> Expression a
scalarderive v e | v == mkExprn e = 1
scalarderive v (Scalar e) = scalar (scalarderive v e)
scalarderive v (Var _ _ _ _ (Just e)) = scalarderive v e
scalarderive _ (Var _ _ _ _ Nothing) = 0
scalarderive v (Sum s _) = pairs2sum $ map dbythis $ sum2pairs s
  where dbythis (f,x) = (f, scalarderive v x)
scalarderive v (Product p i) = pairs2sum $ map dbythis $ product2pairs p
  where dbythis (x,n) = (1, Product p i*toExpression n/x * scalarderive v x)
scalarderive _ (Heaviside _) = error "no scalarderive for Heaviside"
scalarderive v (Cos e) = -sin e * scalarderive v e
scalarderive v (Sin e) = cos e * scalarderive v e
scalarderive v (Exp e) = exp e * scalarderive v e
scalarderive v (Erf e) = 2/sqrt pi*exp (-e**2) * scalarderive v e
scalarderive v (Log e) = scalarderive v e / e
scalarderive _ (Abs _) = error "I didn't think we'd need abs"
scalarderive _ (Signum _) = error "I didn't think we'd need signum"
scalarderive v (Expression e) = scalarderivativeHelper v e

deriveVector :: (Type a, Type b) => Expression b -> Vector a -> Vector a -> Vector b
deriveVector v0 (Vector ddax dday ddaz) (Vector x y z) =
  vector (derive v0 ddax x) (derive v0 dday y) (derive v0 ddaz z)

realspaceGradient :: Vector RealSpace -> Expression RealSpace -> Vector RealSpace
realspaceGradient (Vector vx vy vz) e = vector (derive vx 1 e) (derive vy 1 e) (derive vz 1 e)

derive :: (Type a, Type b) => Expression b -> Expression a -> Expression a -> Expression b
derive v0 dda0 e | Just v <- safeCoerce v0 e,
                   v == e,
                   Just dda <- safeCoerce dda0 v0 = dda
-- The following would treat a scalar derivative of a scalar as
-- scalar, which would make sense.  However, as it turns out, it leads
-- to higher memory use for some reason.  I'm not sure why, but that's
-- why it's disabled for now.

--derive vv@(Scalar v) dda0 (Scalar e)
--  | Just dda <- safeCoerce dda0 vv = dda * Scalar (derive v 1 e)
derive v@(Var _ a b c _) dda0 (Var t _ bb cc (Just e0))
  | Just e <- safeCoerce e0 v,
    Just dda <- safeCoerce dda0 v =
    case isConstant $ derive v dda e of
      Just x -> toExpression x
      Nothing ->
        case isConstant $ derive v 1 e of
          Just x -> toExpression x * dda
          Nothing ->
            if dda == 1 || not (hasExpressionInFFT v e)
            then dda*(Var t ("d" ++ bb ++ "_by_d" ++ a) ("d" ++ bb ++ "_by_d" ++ b)
                      ("\\frac{\\partial "++cc ++"}{\\partial "++c++"}") $ Just $
                      (derive v 1 e))
            else derive v dda e
derive v@(Scalar (Var _ a b c _)) dda0 (Var t _ bb cc (Just e0))
  | Just e <- safeCoerce e0 v,
    Just dda <- safeCoerce dda0 v =
  case isConstant $ derive v dda e of
    Just x -> toExpression x
    Nothing ->
      case isConstant $ derive v 1 e of
        Just x -> toExpression x * dda
        Nothing ->
            if dda == 1 || not (hasExpressionInFFT v e)
            then dda*(Var t ("d" ++ bb ++ "_by_d" ++ a) ("d" ++ bb ++ "_by_d" ++ b)
                      ("\\frac{\\partial "++cc ++"}{\\partial "++c++"}") $ Just $
                      (derive v 1 e))
            else derive v dda e
derive v dda (Var _ _ _ _ (Just e)) = derive v dda e
derive _ _ (Var _ _ _ _ Nothing) = 0
-- In the following line, we avoid dealing with sums that don't
-- include the variable with respect to which we are taking a
-- derivative.
derive v _ (Sum _ s) | not (Set.isSubsetOf (varSet v) s) = 0
derive v dda (Sum s _) = pairs2sum $ map dbythis $ sum2pairs s
  where dbythis (f,x) = (f, derive v dda x)
-- In the following line, we avoid dealing with products that don't
-- include the variable with respect to which we are taking a
-- derivative.
derive v _ (Product _ s) | not (Set.isSubsetOf (varSet v) s) = 0
derive v dda (Product p i) = if False
                             then factorizeSumPairs $ filterNonZero $ map dbythis $ product2pairs p
                             else pairs2sum $ map dbythis $ product2pairs p
  where dbythis (x,n) = (n, derive v (Product p i*dda/x) x)
derive _ _ (Scalar _) = 0 -- FIXME
derive _ _ (Heaviside _) = error "cannot take derivative of Heaviside"
derive v dda (Cos e) = derive v (-dda*sin e) e
derive v dda (Sin e) = derive v (dda*cos e) e
derive v dda (Exp e) = derive v (dda*exp e) e
derive v dda (Log e) = derive v (dda/e) e
derive v dda (Erf e) = derive v (dda*2/sqrt pi*exp (-e**2)) e
derive _ _ (Abs _) = error "I didn't think we'd need abs"
derive _ _ (Signum _) = error "I didn't think we'd need signum"
derive v dda (Expression e) = derivativeHelper v dda e

filterNonZero :: Type a => [(Double, Expression a)] -> [(Double, Expression a)]
filterNonZero = filter nz
  where nz (0,_) = False
        nz (_,x) = x /= 0

hasexpression :: (Type a, Type b) => Expression a -> Expression b -> Bool
hasexpression x e = countexpression x e > 0

hasExprn :: Exprn -> Exprn -> Bool
hasExprn (ES a) (ES b) = hasexpression a b
hasExprn (ES a) (ER b) = hasexpression a b
hasExprn (ES a) (EK b) = hasexpression a b
hasExprn (ER a) (ES b) = hasexpression a b
hasExprn (ER a) (ER b) = hasexpression a b
hasExprn (ER a) (EK b) = hasexpression a b
hasExprn (EK a) (ES b) = hasexpression a b
hasExprn (EK a) (ER b) = hasexpression a b
hasExprn (EK a) (EK b) = hasexpression a b

countexpression :: (Type a, Type b) => Expression a -> Expression b -> Int
countexpression x e = snd $ subAndCount x (s_var "WeAreCounting") e

substitute :: (Type a, Type b) => Expression a -> Expression a -> Expression b -> Expression b
substitute x y e = fst $ subAndCount x y e

substituteE :: Type a => Expression a -> Expression a -> Exprn -> Exprn
substituteE a a' (ES b) = mkExprn $ substitute a a' b
substituteE a a' (ER b) = mkExprn $ substitute a a' b
substituteE a a' (EK b) = mkExprn $ substitute a a' b

countAfterRemoval :: Type a => Expression a -> [Exprn] -> Int
countAfterRemoval v e = Set.size $ mconcat $ map (mapExprn (varsetAfterRemoval v)) e

countAfterRemovalE :: Exprn -> [Exprn] -> Int
countAfterRemovalE a b = mapExprn (\a' -> countAfterRemoval a' b) a

varsetAfterRemoval :: (Type a, Type b) => Expression a -> Expression b -> Set.Set String
-- varsetAfterRemoval v e = varSet (substitute v 2 e) -- This should be equivalent
varsetAfterRemoval x e | mkExprn x == mkExprn e = Set.empty
                       | not (Set.isSubsetOf (varSet x) (varSet e)) = varSet e
varsetAfterRemoval x y
  | EK (Sum xs _) <- mkExprn x,
    EK (Sum es _) <- mkExprn y,
    Just (_, es') <- removeFromMap xs es = varsetAfterRemoval x (map2sum es')
  | ER (Sum xs _) <- mkExprn x,
    ER (Sum es _) <- mkExprn y,
    Just (_, es') <- removeFromMap xs es = varsetAfterRemoval x (map2sum es')
  | ES (Sum xs _) <- mkExprn x,
    ES (Sum es _) <- mkExprn y,
    Just (_, es') <- removeFromMap xs es = varsetAfterRemoval x (map2sum es')
varsetAfterRemoval x y
  | EK (Product xs _) <- mkExprn x,
    EK (Product es _) <- mkExprn y,
    Just (ratio, es') <- removeFromMap xs es,
    abs ratio >= 1 = varsetAfterRemoval x (map2product es')
  | ER (Product xs _) <- mkExprn x,
    ER (Product es _) <- mkExprn y,
    Just (ratio, es') <- removeFromMap xs es,
    abs ratio >= 1 = varsetAfterRemoval x (map2product es')
  | ES (Product xs _) <- mkExprn x,
    ES (Product es _) <- mkExprn y,
    Just (ratio, es') <- removeFromMap xs es,
    abs ratio >= 1 = varsetAfterRemoval x (map2product es')
varsetAfterRemoval x v@(Expression _)
  | EK (Expression (FFT e)) <- mkExprn v = varsetAfterRemoval x e
  | EK (Expression (SetKZeroValue _ e)) <- mkExprn v = varsetAfterRemoval x e
  | ER (Expression (IFFT e)) <- mkExprn v = varsetAfterRemoval x e
  | ES (Expression (Summate e)) <- mkExprn v = varsetAfterRemoval x e
  | otherwise = Set.empty
varsetAfterRemoval x (Sum s _) = Set.unions (map (varsetAfterRemoval x . snd) (sum2pairs s))
varsetAfterRemoval x (Product p _) = Set.unions (map (varsetAfterRemoval x . fst) (product2pairs p))
varsetAfterRemoval x (Heaviside e) = varsetAfterRemoval x e
varsetAfterRemoval x (Cos e) = varsetAfterRemoval x e
varsetAfterRemoval x (Sin e) = varsetAfterRemoval x e
varsetAfterRemoval x (Log e) = varsetAfterRemoval x e
varsetAfterRemoval x (Exp e) = varsetAfterRemoval x e
varsetAfterRemoval x (Erf e) = varsetAfterRemoval x e
varsetAfterRemoval x (Abs e) = varsetAfterRemoval x e
varsetAfterRemoval x (Signum e) = varsetAfterRemoval x e
varsetAfterRemoval x (Var _ _ _ _ (Just e)) = varsetAfterRemoval x e
varsetAfterRemoval _ v@(Var _ _ _ _ Nothing) = varSet v
varsetAfterRemoval x (Scalar e) = varsetAfterRemoval x e


-- removeFromMap is a bit tricky.  It is used to look for a given set
-- of expressions in a Map.Map, and remove them, with some possible
-- "factor".  It is used for removing composite subexpressions from
-- both summations and products.
--
-- This function uses the "do" monad notation, which I generally try
-- to avoid, but I think in this case it helps enough, and the
-- function is complicated enough to start with, that it is
-- worthwhile.
removeFromMap :: Type a => Map.Map (Expression a) Double -> Map.Map (Expression a) Double
                 -> Maybe (Double, Map.Map (Expression a) Double)
removeFromMap xm ym =
  do (xe, xX):_ <- Just $ Map.toList xm
     yX <- Map.lookup xe ym
     let ratio = yX/xX
         filterout ((a,d):rest) emap = do d' <- Map.lookup a emap
                                          if ratio*d == d'
                                            then filterout rest (Map.delete a emap)
                                            else Nothing
         filterout [] emap = Just emap
     ym' <- filterout (Map.toList xm) ym
     Just (ratio, ym')

subAndCount :: (Type a, Type b) => Expression a -> Expression a -> Expression b -> (Expression b, Int)
subAndCount x0 y0 e | Just x <- safeCoerce x0 e,
                      x == e,
                      Just y <- safeCoerce y0 e = (y, 1)
                    | not (Set.isSubsetOf (varSet x0) (varSet e)) = (e, 0) -- quick check
subAndCount x0 y0 e@(Sum es _)
  | Just x@(Sum xs _) <- safeCoerce x0 e,
    Just y <- safeCoerce y0 e,
    Just (ratio, es') <- removeFromMap xs es,
    (e'',n) <- subAndCount x y (map2sum es') = (e'' + toExpression ratio*y, n+1)
subAndCount x0 y0 e@(Product es _)
  | Just x@(Product xs _) <- safeCoerce x0 e,
    Just y <- safeCoerce y0 e,
    Just (ratio, es') <- removeFromMap xs es,
    abs ratio >= 1,
    (e'',n) <- subAndCount x y (map2product es') = (e'' * y**(toExpression ratio), n+1)
subAndCount x y (Expression v) = subAndCountHelper x y v
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
subAndCount x y (Heaviside e)   = (heaviside e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Cos e)   = (cos e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Sin e)   = (sin e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Log e)   = (log e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Exp e)   = (exp e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Erf e)   = (erf e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Abs e)   = (abs e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Signum e) = (signum e', n)
    where (e', n) = subAndCount x y e
subAndCount x y (Var t a b c (Just e)) = (Var t a b c (Just e'), n)
    where (e', n) = subAndCount x y e
subAndCount _ _ v@(Var _ _ _ _ Nothing) = (v, 0)
subAndCount x y (Scalar e) = (scalar e', n)
    where (e', n) = subAndCount x y e


newtype MkBetter a = MB (Maybe (Expression a -> MkBetter a, Expression a))

findRepeatedSubExpression :: Type a => Expression a -> MkBetter a
findRepeatedSubExpression everything = frse everything
  where frse (Expression _) = MB Nothing
        frse (Var _ _ _ _ Nothing) = MB Nothing
        frse (Scalar (Var _ _ _ _ Nothing)) = MB Nothing
        frse e | Just _ <- isConstant e = MB Nothing
        frse e | mytimes > 1 = MB (Just (makebetter e mytimes, e))
          where mytimes = countexpression e everything
        frse x@(Sum s _) | MB (Just (better,_)) <- frs (sum2pairs s) = better x
          where frs ((_,y):ys) = case frse y of
                                   MB Nothing -> frs ys
                                   se -> se
                frs [] = MB Nothing
        frse x@(Product s _) | MB (Just (better,_)) <- frs (product2pairs s) = better x
          where frs ((y,n):ys) = case frse y of
                                   MB Nothing -> frs ys
                                   MB (Just (better,_)) -> better (y ** toExpression n)
                frs [] = MB Nothing
        frse x@(Cos e) | MB (Just (better, _)) <- frse e = better x
        frse x@(Sin e) | MB (Just (better, _)) <- frse e = better x
        frse x@(Log e) | MB (Just (better, _)) <- frse e = better x
        frse x@(Exp e) | MB (Just (better, _)) <- frse e = better x
        frse x@(Erf e) | MB (Just (better, _)) <- frse e = better x
        frse x@(Var _ _ _ _ (Just e)) | MB (Just (better, _)) <- frse e = better x
        frse _ = MB Nothing
        makebetter :: Type a => Expression a -> Int -> Expression a -> MkBetter a
        makebetter e n e' = if n' >= n then MB (Just (makebetter e' n', e'))
                                       else MB (Just (makebetter e n, e))
                 where n' = countexpression e' everything

-- A monoid is (just about) any data type that allows you to combine
-- two elements together to create a third element of the same type.
-- Sets are like this (union combines two sets), and so is Maybe,
-- where we want to short-cut to find the first "Just" element.
searchMonoid :: (Type a, Monoid c) => (forall b. Type b => Expression b -> c)
                    -> Expression a -> c
searchMonoid f x@(Var _ _ _ _ Nothing) = f x
searchMonoid f x@(Var _ _ _ _ (Just e)) = f x `mappend` searchMonoid f e
searchMonoid f x@(Scalar e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Heaviside e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Cos e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Sin e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Exp e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Erf e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Log e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Abs e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Signum e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Product p _) = f x `mappend` mconcat (map (searchMonoid f . fst) $ product2pairs p)
searchMonoid f x@(Sum s _) = f x `mappend` mconcat (map (searchMonoid f . snd) $ sum2pairs s)
searchMonoid f x@(Expression e) = f x `mappend` searchHelper (searchMonoid f) e

findNamedScalars :: Type b => Expression b -> Set.Set String
findNamedScalars = searchMonoid helper
  where helper (Var _ _ b _ (Just e)) | ES _ <- mkExprn e = Set.singleton b
        helper _ = Set.empty

findOrderedInputs :: Type a => Expression a -> [Exprn]
findOrderedInputs e = fst $ helper Set.empty (findInputs e)
  where helper sofar inps
          | x:_ <- filter isok $ Set.toList inps =
            case helper (Set.insert x sofar) (Set.delete x inps) of
              (out, sofar') -> (x:out, sofar')
          | x:_ <- filter lessok $ Set.toList inps =
              case helper sofar (Set.difference (Set.delete x $ findInputsE x) sofar) of
                (out, sofar') ->
                  case helper sofar' (Set.difference inps sofar') of
                    (out', sofar'') -> (out ++ out', sofar'')
          | otherwise = ([], sofar)
            where isok ee = Set.size (Set.difference (findInputsE ee) sofar) == 1 && not (Set.member ee sofar)
                  lessok ee = not (Set.member ee sofar)
        findInputsE (ES ee) = findInputs ee
        findInputsE (EK ee) = findInputs ee
        findInputsE (ER ee) = findInputs ee

findInputs :: Type b => Expression b -> Set.Set Exprn
findInputs = searchMonoid helper
  where helper e | ER (Var _ _ _ _ Nothing) <- mkExprn e,
                   not (Set.member (mkExprn e) grid_description) = Set.insert (mkExprn e) grid_description
        helper e@(Var _ _ _ _ Nothing) = Set.singleton (mkExprn e)
        helper _ = Set.empty
        grid_description = Set.fromList [ES numx, ES numy, ES numz,
                                         ES $ xhat `dot` lat1,
                                         ES $ yhat `dot` lat2,
                                         ES $ zhat `dot` lat3]

findTransforms :: Type b => Expression b -> [Expression KSpace]
findTransforms = Set.toList . searchMonoid helper
  where helper e | EK ee@(Expression (SphericalFourierTransform _ _)) <- mkExprn e = Set.singleton ee
        helper _ = Set.empty

dV :: Expression RealSpace
dV = scalar dVscalar

dVscalar :: Expression Scalar
dVscalar = var "dV" "\\Delta V" $ volume / numx / numy / numz

dr :: Expression KSpace
dr = scalar $ var "dr" "\\Delta r" $ dVscalar ** (1.0/3)

volume :: Type a => Expression a
volume = scalar $ protect "volume" "volume" $ lat1 `dot` (lat2 `cross` lat3)

numx, numy, numz :: Type a => Expression a
numx = s_var "Nx"
numy = s_var "Ny"
numz = s_var "Nz"

lat1, lat2, lat3, rlat1, rlat2, rlat3 :: Vector Scalar
lat1 = Vector (s_var "a1") 0 0
lat2 = Vector 0 (s_var "a2") 0
lat3 = Vector 0 0 (s_var "a3")


rlat1 = (lat2 `cross` lat3) /. (lat1 `dot` (lat2 `cross` lat3))
rlat2 = (lat3 `cross` lat1) /. (lat1 `dot` (lat2 `cross` lat3))
rlat3 = (lat1 `cross` lat2) /. (lat1 `dot` (lat2 `cross` lat3))

nameVector :: Type a => String -> Vector a -> Vector a
nameVector n (Vector x y z) = Vector (nn "x" x) (nn "y" y) (nn "z" z)
  where nn a s = (n++a) === s
\end{code}
