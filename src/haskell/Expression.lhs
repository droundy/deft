\begin{code}
{-# LANGUAGE PatternGuards, Rank2Types #-}

module Expression (Exprn(..),
                   RealSpace(..), r_var, dV,
                   KSpace(..), k_var, imaginary, kx, ky, kz, k, ksqr, setkzero,
                   Scalar(..),
                   ThreeVector, t_var, cross, dot,
                   fft, ifft, integrate, grad, derive,
                   Expression(..), joinFFTs, (===), var, protect,
                   Type(..), Code(..), IsTemp(..),
                   makeHomogeneous, isConstant,
                   setZero, cleanvars, factorize, factorOut,
                   cleanvarsE, initializeE, freeE, newinitializeE, newfreeE,
                   nameE, newdeclareE,
                   sum2pairs, pairs2sum, codeStatementE, newcodeStatementE,
                   product2pairs, pairs2product, product2denominator,
                   hasFFT, hasexpression, hasExprn,
                   searchExpression, searchExpressionDepthFirst,
                   findRepeatedSubExpression, findNamedScalars, findOrderedInputs, findInputs,
                   countexpression, substitute, countAfterRemoval,
                   substituteE, countAfterRemovalE, countVarsE,
                   mapExpressionShortcut, -- just to avoid unused warning
                   countVars, varSet, varSetE)
    where

import Debug.Trace

import qualified Data.Map as Map
import qualified Data.Set as Set
import LatexDouble ( latexDouble )

data RealSpace = IFFT (Expression KSpace)
               deriving ( Eq, Ord, Show )
data KSpace = Delta | -- handy for FFT of homogeneous systems
              Kx | Ky | Kz |
              SetKZeroValue (Expression KSpace) (Expression KSpace) |
              FFT (Expression RealSpace)
            deriving ( Eq, Ord, Show )
data Scalar = Integrate (Expression RealSpace)
            | Component Int (Expression ThreeVector)
            deriving ( Eq, Ord, Show )

data ThreeVector = Permute Int (Expression ThreeVector)
                 | ThreeVector (Expression Scalar) (Expression Scalar) (Expression Scalar)
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

unitvector :: Int -> Expression ThreeVector
unitvector 0 = threeVector 1 0 0
unitvector 1 = threeVector 0 1 0
unitvector 2 = threeVector 0 0 1
unitvector n = unitvector (n `mod` 3)

threeVector :: Expression Scalar -> Expression Scalar -> Expression Scalar -> Expression ThreeVector
threeVector x y z | x == 0 && y == 0 && z == 0 = 0
                  | otherwise = Expression (ThreeVector x y z)

permute :: Int -> Expression ThreeVector -> Expression ThreeVector
permute n e | n `mod` 3 == 0 = e
permute n (Expression (Permute n' e)) = permute (n+n') e
permute n (Expression (ThreeVector x y z)) = permute (n-1) $ Expression (ThreeVector y z x)
permute n e@(Var _ _ _ _ Nothing) = Expression $ Permute n e
permute n (Var t a b c (Just e)) = Var t a b c $ Just $ permute n e
permute _ (Scalar e) = scalar e
permute _ (Cos _) = error "can't take cosine of 3vec"
permute _ (Sin _) = error "can't take sin of 3vec"
permute _ (Exp _) = error "can't take exp of 3vec"
permute _ (Log _) = error "can't take log of 3vec"
permute _ (Abs _) = error "can't take abs of 3vec"
permute _ (Signum _) = error "can't take signum of 3vec"
permute n (Product p _) = pairs2product $ map ff $ product2pairs p
  where ff (e,x) = (permute n e, x)
permute n (Sum s _) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, permute n y)

component :: Type a => Int -> Expression ThreeVector -> Expression a
component n e | n `mod` 3 /= n = component (n `mod` 3) e
component n (Expression (Permute n' e)) = component (n+n') e
component 0 (Expression (ThreeVector x _ _)) = fromScalar x
component 1 (Expression (ThreeVector _ y _)) = fromScalar y
component 2 (Expression (ThreeVector _ _ z)) = fromScalar z
component n (Expression (ThreeVector _ _ _)) = error ("invalid component: "++show n)
component n v@(Var _ _ _ _ Nothing) = scalar $ Expression (Component n v)
component n (Var t _ b c (Just e)) =
  Var t (b++["xyz" !! n]) (b++["xyz" !! n]) (c++"_{"++ "xyz" !! n : "}") $ Just $ component n e
component _ (Scalar e) = scalar e
component _ (Cos _) = error "can't take cosine of 3vec"
component _ (Sin _) = error "can't take sin of 3vec"
component _ (Exp _) = error "can't take exp of 3vec"
component _ (Log _) = error "can't take log of 3vec"
component _ (Abs _) = error "can't take abs of 3vec"
component _ (Signum _) = error "can't take signum of 3vec"
component n (Product p _) = pairs2product $ map ff $ product2pairs p
  where ff (e,x) = (component n e, x)
component n (Sum s _) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, component n y)

cross :: Expression ThreeVector -> Expression ThreeVector -> Expression ThreeVector
cross a b = permute 1 a * permute 2 b - permute 2 a * permute 1 b
dot :: Type a => Expression ThreeVector -> Expression ThreeVector -> Expression a
dot a b = component 0 a * component 0 b + component 1 a * component 1 b + component 2 a * component 2 b
infixl 7 `cross`, `dot`

instance Code ThreeVector where
  codePrec _ (Permute n (Var _ _ v _ Nothing)) = showString (v ++ "[i+"++show n++"]")
  codePrec _ e = error ("Bad permute: " ++ show e)
  latexPrec _ (Permute n (Var _ _ _ v Nothing)) = showString (v ++ "_{+"++show n++"}")
  latexPrec _ e = error ("Bad latex permute: " ++ show e)
instance Type ThreeVector where
  amThreeVector _ = True
  mkExprn = E3
  derivativeHelper v ddv (Permute n e) = derive v (permute (n+2) ddv) e
  derivativeHelper v ddv (ThreeVector x y z) = derive v (component 0 ddv) x +
                                               derive v (component 1 ddv) y +
                                               derive v (component 2 ddv) z
  scalarderivativeHelper v (Permute n e) = permute n (scalarderive v e)
  scalarderivativeHelper v (ThreeVector x y z) = threeVector (scalarderive v x) (scalarderive v y) (scalarderive v z)
  zeroHelper v (Permute n e) = permute n $ setZero v e
  zeroHelper v (ThreeVector x y z) = threeVector (setZero v x) (setZero v y) (setZero v z)
  codeStatementHelper (Var _ _ a _ _) op e =
    unlines ["\tfor (int i=0; i<3; i++) {",
             "\t\t" ++ a ++ op ++ code e ++ ";",
             "\t}"]
  codeStatementHelper _ _ _ = error "bad input to newcodeStatementHelper"
  newcodeStatementHelper (Var _ _ a _ _) op e =
    unlines ["\tfor (int i=0; i<3; i++) {",
             "\t\t" ++ a ++ op ++ newcode e ++ ";",
             "\t}"]
  newcodeStatementHelper _ _ _ = error "bad input to newcodeStatementHelper"
  initialize _ = error "initialize 3vec"
  free e = error $ trace "free error" ("free error " ++ show e)
  newdeclare _ = "Vector"
  newfree e = error $ trace "free error" ("free error " ++ show e)
  toScalar _ = error "how do I make a 3vector homogeneous?"
  mapExpressionHelper' f (Permute n e) = permute n $ f e
  mapExpressionHelper' f (ThreeVector x y z) = threeVector (f x) (f y) (f z)
  subAndCountHelper x y (Permute n e) = case subAndCount x y e of
                                        (e', nn) -> (permute n e', nn)
  subAndCountHelper a b (ThreeVector x y z) =
    case (subAndCount a b x, subAndCount a b y, subAndCount a b z) of
      ((x', nx), (y', ny), (z', nz)) -> (threeVector x' y' z', nx+ny+nz)
  searchHelper f (Permute _ e) = f e
  searchHelper f (ThreeVector x y z) = mconcat [f x, f y, f z]
  safeCoerce a _ = case mkExprn a of
                    E3 a' -> Just a'
                    _ -> Nothing

instance Code RealSpace where
  codePrec _ (IFFT (Var _ _ ksp _ Nothing)) = showString ("ifft(gd, " ++ksp++ ")")
  codePrec _ (IFFT ke) = showString "ifft(gd, " . codePrec 0 ke . showString ")"
  latexPrec _ (IFFT ke) = showString "\\text{ifft}\\left(" . latexPrec 0 ke . showString "\\right)"
instance Type RealSpace where
  amRealSpace _ = True
  mkExprn = ER
  derivativeHelper v ddr (IFFT ke) = derive v (fft ddr) (kinversion ke)
  scalarderivativeHelper v (IFFT ke) = ifft (scalarderive v ke)
  zeroHelper v (IFFT ke) = ifft (setZero v ke)
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
    a ++ op ++ "ifft(gd, " ++ v ++ ");\n"
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
              MB Nothing -> "\t\t" ++ newcode a ++ op ++ newcode x ++ ";"
  initialize (Var IsTemp _ x _ Nothing) = "VectorXd " ++ x ++ "(gd.NxNyNz);"
  initialize _ = error "VectorXd output(gd.NxNyNz);"
  free (Var IsTemp _ x _ Nothing) = x ++ ".resize(0); // Realspace"
  free e = error $ trace "free error" ("free error " ++ show e)
  newdeclare _ = "Vector"
  newinitialize (Var _ _ x _ Nothing) = "Vector " ++ x ++ "(Nx*Ny*Nz);"
  newinitialize _ = error "oops newinitializeE"
  newfree (Var IsTemp _ x _ Nothing) = x ++ ".free(); // Realspace"
  newfree e = error $ trace "free error" ("free error " ++ show e)
  toScalar (IFFT ke) = makeHomogeneous ke
  mapExpressionHelper' f (IFFT ke) = ifft (f ke)
  subAndCountHelper x y (IFFT ke) = case subAndCount x y ke of (ke', n) -> (ifft ke', n)
  searchHelper f (IFFT e) = f e
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
  latexPrec _ Kx = showString "k_{x}"
  latexPrec _ Ky = showString "k_{y}"
  latexPrec _ Kz = showString "k_{z}"
  latexPrec _ Delta = showString "\\delta(k)"
  latexPrec p (SetKZeroValue _ val) = latexPrec p val
  latexPrec _ (FFT r) = showString "\\text{fft}\\left(" . latexPrec 0 r . showString "\\right)"
instance Type KSpace where
  amKSpace _ = True
  mkExprn = EK
  derivativeHelper v ddk (FFT r) = derive v (ifft ddk) r
  derivativeHelper v ddk (SetKZeroValue _ e) = derive v (setkzero 0 ddk) e -- FIXME: how best to handle k=0 derivative?
  derivativeHelper _ _ Kx = 0
  derivativeHelper _ _ Ky = 0
  derivativeHelper _ _ Kz = 0
  derivativeHelper _ _ Delta = 0
  scalarderivativeHelper v (FFT r) = fft (scalarderive v r)
  scalarderivativeHelper v (SetKZeroValue z e) = setkzero (scalarderive v z) (scalarderive v e)
  scalarderivativeHelper _ _ = 0
  zeroHelper v (FFT r) = fft (setZero v r)
  zeroHelper _ Kx = Expression Kx
  zeroHelper _ Ky = Expression Ky
  zeroHelper _ Kz = Expression Kz
  zeroHelper _ Delta = Expression Delta
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
             "\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);",
             codes (1 :: Int) e,
             "\t}"]
      where codes n x = case findRepeatedSubExpression x of
              MB (Just (_,x')) -> "\t\tconst complex t"++ show n ++ " = " ++ code x' ++ ";\n" ++
                                  codes (n+1) (substitute x' (s_var ("t"++show n)) x)
              MB Nothing -> "\t\t" ++ a ++ "[i]" ++ op ++ code x ++ ";"
            setzero = case code $ setKequalToZero e of
                      "0" -> a ++ "[0]" ++ op ++ "0;"
                      k0code -> unlines ["\t{",
                                         "\t\tconst int i = 0;",
                                         "\t\tconst Reciprocal k_i = Reciprocal(0,0,0);",
                                         "\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);",
                                         "\t\t" ++ a ++ "[0]" ++ op ++ k0code ++ ";",
                                         "\t}"]

  codeStatementHelper _ _ _ = error "Illegal input to codeStatementHelper for kspace"


  newcodeStatementHelper (Var _ _ a _ _) op (Expression (FFT (Var _ _ v _ Nothing))) =
    a ++ op ++ "fft(gd, " ++ v ++ ");\n"
  newcodeStatementHelper _ _ (Expression (FFT _)) =
    error "It is a bug to generate newcode for a non-var input to fft"
  newcodeStatementHelper a op (Var _ _ _ _ (Just e)) = newcodeStatementHelper a op e
  newcodeStatementHelper (Var _ _ a _ _) op e =
    unlines [setzero,
             "\tfor (int i=1; i<Nx*Ny*Nz/2; i++) {",
             "\t\tconst int z = i % gd.NzOver2;",
             "\t\tconst int n = (i-z)/gd.NzOver2;",
             "\t\tconst int y = n % gd.Ny;",
             "\t\tconst int xa = (n-y)/gd.Ny;",
             "\t\tconst RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);",
             "\t\tconst Reciprocal k_i = gd.Lat.toReciprocal(rvec);",
             "\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);",
             newcodes (1 :: Int) e,
             "\t}"]
      where newcodes n x = case findRepeatedSubExpression x of
              MB (Just (_,x')) -> "\t\tconst complex t"++ show n ++ " = " ++ newcode x' ++ ";\n" ++
                                  newcodes (n+1) (substitute x' (s_var ("t"++show n)) x)
              MB Nothing -> "\t\t" ++ a ++ "[i]" ++ op ++ newcode x ++ ";"
            setzero = case newcode $ setKequalToZero e of
                      "0" -> a ++ "[0]" ++ op ++ "0;"
                      k0newcode -> unlines ["\t{",
                                         "\t\tconst int i = 0;",
                                         "\t\tconst Reciprocal k_i = Reciprocal(0,0,0);",
                                         "\t\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);",
                                         "\t\t" ++ a ++ "[0]" ++ op ++ k0newcode ++ ";",
                                         "\t}"]

  newcodeStatementHelper _ _ _ = error "Illegal input to newcodeStatementHelper for kspace"
  initialize (Var IsTemp _ x _ Nothing) = "VectorXcd " ++ x ++ "(gd.NxNyNzOver2);"
  initialize _ = error "VectorXcd output(gd.NxNyNzOver2);"
  free (Var IsTemp _ x _ Nothing) = x ++ ".resize(0); // KSpace"
  free _ = error "free error"
  newdeclare _ = "Vector"
  newfree (Var IsTemp _ x _ Nothing) = x ++ ".free(); // KSpace"
  newfree _ = error "free error"
  toScalar Delta = 1
  toScalar Kx = s_var "_kx"
  toScalar Ky = 0
  toScalar Kz = 0
  toScalar (SetKZeroValue val _) = makeHomogeneous val
  toScalar (FFT e) = makeHomogeneous e
  mapExpressionHelper' f (FFT e) = fft (f e)
  mapExpressionHelper' f (SetKZeroValue z v) = setkzero (f z) (f v)
  mapExpressionHelper' _ kk = Expression kk
  subAndCountHelper x y (FFT e) = case subAndCount x y e of (e', n) -> (fft e', n)
  subAndCountHelper x y (SetKZeroValue z e) = (setkzero z' e', n1+n2)
        where (z',n1) = subAndCount x y z
              (e',n2) = subAndCount x y e
  subAndCountHelper _ _ Kx = (kx, 0)
  subAndCountHelper _ _ Ky = (ky, 0)
  subAndCountHelper _ _ Kz = (kz, 0)
  subAndCountHelper _ _ Delta = error "unhandled case: Delta"
  searchHelper f (FFT e) = f e
  searchHelper f (SetKZeroValue _ e) = f e
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

cleanvarsE :: Exprn -> Exprn
cleanvarsE (ES e) = ES $ cleanvars e
cleanvarsE (EK e) = EK $ cleanvars e
cleanvarsE (ER e) = ER $ cleanvars e
cleanvarsE (E3 e) = E3 $ cleanvars e

mapExpression' :: Type a => (forall b. Type b => Expression b -> Expression b) -> Expression a -> Expression a
mapExpression' f (Var IsTemp a b c (Just e)) = f $ Var IsTemp a b c $ Just $ mapExpression' f e
mapExpression' f e@(Var CannotBeFreed _ _ _ (Just _)) = f e
mapExpression' f (Var tt c v t Nothing) = f $ Var tt c v t Nothing
mapExpression' f (Scalar e) = f $ Scalar (mapExpression' f e)
mapExpression' f (Cos e) = f $ cos (mapExpression' f e)
mapExpression' f (Sin e) = f $ sin (mapExpression' f e)
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
mapExpressionShortcut f (Cos e) = cos (mapExpressionShortcut f e)
mapExpressionShortcut f (Sin e) = sin (mapExpressionShortcut f e)
mapExpressionShortcut f (Exp e) = exp (mapExpressionShortcut f e)
mapExpressionShortcut f (Log e) = log (mapExpressionShortcut f e)
mapExpressionShortcut f (Abs e) = abs (mapExpressionShortcut f e)
mapExpressionShortcut f (Signum e) = signum (mapExpressionShortcut f e)
mapExpressionShortcut f (Product p _) = pairs2product $ map ff $ product2pairs p
  where ff (e,n) = (mapExpressionShortcut f e, n)
mapExpressionShortcut f (Sum s _) = pairs2sum $ map ff $ sum2pairs s
  where ff (x,y) = (x, mapExpressionShortcut f y)
mapExpressionShortcut f (Expression x) = mapExpressionHelper' (mapExpressionShortcut f) x

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
searchExpression i f (Cos e) = searchExpression i f e
searchExpression i f (Sin e) = searchExpression i f e
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
searchExpressionDepthFirst i f x@(Cos e) = searchExpressionDepthFirst i f e `mor` f x
searchExpressionDepthFirst i f x@(Sin e) = searchExpressionDepthFirst i f e `mor` f x
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
  where helper (Var IsTemp b c d (Just e)) | ES _ <- mkExprn e = Var IsTemp b c d (Just e)
                                           | otherwise = e
        helper e = e

isEven :: (Type a) => Exprn -> Expression a -> Double
isEven v e | v == mkExprn e = -1
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
isEven v e = case mkExprn e of
             EK (Expression (FFT r)) -> isEven v r
             EK _ -> 1
             ER (Expression (IFFT ks)) -> isEven v ks
             ES (Expression (Integrate x)) -> isEven v x
             _ -> 1 -- Expression _.  Technically, it might be good to recurse into this

-- expand does a few terms in a taylor expansion, and is intended only
-- to be a tool for setting a variable (k, in particular) to zero.
expand :: Type a => Expression a -> Expression a -> Expression a
expand v (Var a b c d (Just e)) = Var a b c d (Just $ expand v e)
expand _ e@(Var _ _ _ _ Nothing) = e
expand _ (Scalar e) = Scalar e
expand v (Cos e) = if setZero (mkExprn v) e == 0
                   then 1 - e'**2/2 + e'**4/4/3/2
                   else cos e'
     where e' = expand v e
expand v (Sin e) = if setZero (mkExprn v) e == 0
                   then e' - e'**3/3/2
                   else sin e'
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
         then error ("L'Hopital's rule failure:\n" ++
                     latex n ++ "\n /\n  " ++ latex d ++ "\n\n\n"
                     ++ latex (Product p i) ++ "\n\n\n" ++ latex zn)
         else setZero v (scalarderive v n / scalarderive v d)
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
  codePrec _ (Component n (Var _ _ v _ Nothing)) = showString (v ++ "[" ++ show n ++ "]")
  codePrec _ (Component _ e) = error ("bad component: " ++ show e)
  latexPrec _ (Integrate r) = showString "\\int " . latexPrec 0 r
  latexPrec _ (Component n (Var _ _ _ v Nothing)) = showString (v ++ "_{" ++ "xyz" !! n : "}")
  latexPrec _ (Component _ e) = error ("bad component: " ++ show e)
instance Type Scalar where
  s_var ("complex(0,1)") = Var CannotBeFreed "complex(0,1)" "complex(0,1)" "i" Nothing
  s_var v@['d',_] = Var CannotBeFreed v v v Nothing -- for differentials
  s_var vv@(a:v@(_:_)) = Var CannotBeFreed vv vv (a : '_' : '{' : v ++ "}") Nothing
  s_var v = Var CannotBeFreed v v v Nothing
  s_tex vv tex = Var CannotBeFreed vv vv tex Nothing
  amScalar _ = True
  mkExprn = ES
  derivativeHelper v dds (Integrate e) = derive v (scalar $ dds*dV) e
  derivativeHelper v dds (Component n e) = derive v (unitvector n * scalar dds) e
  scalarderivativeHelper v (Integrate e) = integrate (scalarderive v e)
  scalarderivativeHelper v (Component n e) = component n (scalarderive v e)
  zeroHelper v (Integrate e) = integrate (setZero v e)
  zeroHelper v (Component n e) = component n (setZero v e)
  codeStatementHelper a " = " (Var _ _ _ _ (Just e)) = codeStatementHelper a " = " e
  codeStatementHelper a " = " (Expression (Integrate e)) =
    code a ++ " = 0;\n\tfor (int i=0; i<gd.NxNyNz; i++) {\n\t\t" ++
    code a ++ " += " ++ code (e * s_var "gd.dvolume") ++
    ";\n\t}\n"
  codeStatementHelper _ op (Expression (Integrate _)) = error ("Haven't implemented "++op++" for integrate...")
  codeStatementHelper a op e = code a ++ op ++ code e ++ ";"

  newcodeStatementHelper a " = " (Var _ _ _ _ (Just e)) = newcodeStatementHelper a " = " e
  newcodeStatementHelper a " = " (Expression (Integrate e)) =
    newcode a ++ " = 0;\n\tfor (int i=0; i<Nx*Ny*Nz; i++) {\n\t\t" ++
    newcode a ++ " += " ++ newcode (e * dV) ++
    ";\n\t}\n"
  newcodeStatementHelper _ op (Expression (Integrate _)) = error ("Haven't implemented "++op++" for integrate...")
  newcodeStatementHelper a op e = newcode a ++ op ++ newcode e ++ ";"
  initialize (Var _ _ x _ Nothing) = "double " ++ x ++ " = 0;\n"
  initialize v = error ("bug in initialize Scalar: "++show v)
  newdeclare _ = "double"
  toScalar (Integrate r) = makeHomogeneous r
  toScalar (Component n e) = component n e
  fromScalar = id
  mapExpressionHelper' f (Integrate e) = integrate (f e)
  mapExpressionHelper' f (Component n e) = component n (f e)
  subAndCountHelper x y (Integrate e) = case subAndCount x y e of (e', n) -> (integrate e', n)
  subAndCountHelper x y (Component nn e) = case subAndCount x y e of (e', n) -> (component nn e', n)
  searchHelper f (Integrate e) = f e
  searchHelper f (Component _ e) = f e
  safeCoerce a _ = case mkExprn a of
                    ES a' -> Just a'
                    _ -> Nothing

scalar :: Type a => Expression Scalar -> Expression a
scalar (Scalar e) = scalar e
scalar e | Just c <- isConstant e = toExpression c
scalar e = fromScalar e

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

t_var :: String -> Expression ThreeVector
t_var v@([_]) = Var CannotBeFreed (v++"[i]") v ("\\vec{"++v++"}") Nothing
t_var (a:v) = Var CannotBeFreed (a:v++"[i]") (a:v) ("\\vec{"++[a]++"_{"++v++"}}") Nothing
t_var "" = error "t_var needs non-empty string"

imaginary :: Expression KSpace
imaginary = Var CannotBeFreed "complex(0,1)" "complex(0,1)" "i" Nothing

infix 4 ===

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

-- An expression statement holds a mathematical expression, which is
-- be any one of several different types: RealSpace (as in n(\vec{r})),
-- KSpace (as in w(\vec{k})), or Scalar (as in an ordinary scalar value
-- like k_BT.
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

-- An "E" type is *any* type of expression.  This is useful when we
-- want to specify subexpressions, for instance, when we might not
-- want to care in advance which sort of subexpression it is.  Also
-- useful for comparing two expressions that might not be the same
-- type.
data Exprn = EK (Expression KSpace) | ER (Expression RealSpace) | ES (Expression Scalar)
           | E3 (Expression ThreeVector)
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
  codePrec p (ES e) = codePrec p e
  codePrec p (EK e) = codePrec p e
  codePrec p (ER e) = codePrec p e
  codePrec p (E3 e) = codePrec p e
  latexPrec p (ES e) = latexPrec p e
  latexPrec p (EK e) = latexPrec p e
  latexPrec p (ER e) = latexPrec p e
  latexPrec p (E3 e) = latexPrec p e
instance (Type a, Code a) => Code (Expression a) where
  codePrec _ (Var _ c _ _ Nothing) = showString c
  codePrec p (Var _ _ _ _ (Just e)) = codePrec p e
  codePrec p (Scalar x) = codePrec p x
  codePrec p (Expression x) = codePrec p x
  codePrec _ (Cos x) = showString "cos(" . codePrec 0 x . showString ")"
  codePrec _ (Sin x) = showString "sin(" . codePrec 0 x . showString ")"
  codePrec _ (Exp x) = showString "exp(" . codePrec 0 x . showString ")"
  codePrec _ (Log x) = showString "log(" . codePrec 0 x . showString ")"
  codePrec _ (Abs x) = showString "fabs(" . codePrec 0 x . showString ")"
  codePrec _ (Signum _) = undefined
  codePrec _ (Product p i) | Product p i == 1 = showString "1"
  codePrec pree (Product p _) = showParen (pree > 7) $
                           if den == 1
                           then codesimple num
                           else codesimple num . showString "/" . codePrec 8 den
    where codesimple [] = showString "1"
          codesimple [(a,n)] = codee a n
          codesimple [(a,n),(b,m)] = codee a n . showString "*" . codee b m
          codesimple ((a,n):es) = codee a n . showString "*" . codesimple es
          num = product2numerator_pairs p
          den = product2denominator p
          codee _ 0 = showString "1" -- this shouldn't happen...
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
  codePrec _ (Sum s i) | Sum s i == 0 = showString "0"
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
  latexPrec _ (Cos x) = showString "\\cos(" . latexPrec 0 x . showString ")"
  latexPrec _ (Sin x) = showString "\\sin(" . latexPrec 0 x . showString ")"
  latexPrec _ (Exp x) = showString "\\exp\\left(" . latexPrec 0 x . showString "\\right)"
  latexPrec _ (Log x) = showString "\\log(" . latexPrec 0 x . showString ")"
  latexPrec _ (Abs x) = showString "\\left|" . latexPrec 0 x . showString "\\right|"
  latexPrec _ (Signum _) = undefined
  latexPrec p (Product x _) | Map.size x == 1 && product2denominator x == 1 =
    case product2pairs x of
      [(_,0)] -> showString "1" -- this shouldn't happen...
      [(_, n)] | n < 0 -> error "shouldn't have negative power here"
      [(e, 1)] ->   latexPrec p e
      [(e, 0.5)] -> showString "\\sqrt{" . latexPrec 0 e . showString "}"
      [(e, n)] -> latexPrec 8 e . showString ("^{" ++ latexDouble n ++ "}")
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
  amThreeVector :: Expression a -> Bool
  amThreeVector _ = False
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
  newinitialize e@(Var _ _ x _ _) = newdeclare e ++ " " ++ x
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
initializeE (E3 e) = initialize e
freeE :: Exprn -> String
freeE (ES e) = free e
freeE (ER e) = free e
freeE (EK e) = free e
freeE (E3 e) = free e

newinitializeE :: Exprn -> String
newinitializeE (ES e) = newinitialize e
newinitializeE (EK e) = newinitialize e
newinitializeE (ER e) = newinitialize e
newinitializeE (E3 e) = newinitialize e

newdeclareE :: Exprn -> String
newdeclareE (ES e) = newdeclare e
newdeclareE (EK e) = newdeclare e
newdeclareE (ER e) = newdeclare e
newdeclareE (E3 e) = newdeclare e

nameE :: Exprn -> String
nameE (ES (Var _ _ v _ Nothing)) = v
nameE (EK (Var _ _ v _ Nothing)) = v
nameE (ER (Var _ _ v _ Nothing)) = v
nameE (E3 (Var _ _ v _ Nothing)) = v
nameE e = show e

newfreeE :: Exprn -> String
newfreeE (ES e) = newfree e
newfreeE (ER e) = newfree e
newfreeE (EK e) = newfree e
newfreeE (E3 e) = newfree e

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
countVars s = Set.size $ varSet s

countVarsE :: Exprn -> Int
countVarsE = Set.size . varSetE

varSetE :: Exprn -> Set.Set String
varSetE (ES e) = varSet e
varSetE (EK e) = varSet e
varSetE (ER e) = varSet e
varSetE (E3 e) = varSet e

varSet :: Type a => Expression a -> Set.Set String
varSet e@(Expression _) = case mkExprn e of
                          EK (Expression (FFT e')) -> varSet e'
                          EK (Expression (SetKZeroValue _ e')) -> varSet e'
                          ER (Expression (IFFT e')) -> varSet e'
                          ES (Expression (Integrate e')) -> varSet e'
                          _ -> Set.empty
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

hasFFT :: Type a => Expression a -> Bool
hasFFT e@(Expression _) = case mkExprn e of
  EK (Expression (FFT _)) -> True
  EK (Expression (SetKZeroValue z e')) -> hasFFT e' || hasFFT z
  ER (Expression (IFFT _)) -> True
  ES (Expression (Integrate _)) -> True -- a bit weird... rename this function?
  _ -> False
hasFFT (Var _ _ _ _ (Just e)) = hasFFT e
hasFFT (Sum s _) = or $ map (hasFFT . snd) (sum2pairs s)
hasFFT (Product p _) = or $ map (hasFFT . fst) (product2pairs p)
hasFFT (Sin e) = hasFFT e
hasFFT (Cos e) = hasFFT e
hasFFT (Log e) = hasFFT e
hasFFT (Exp e) = hasFFT e
hasFFT (Abs e) = hasFFT e
hasFFT (Signum e) = hasFFT e
hasFFT (Var _ _ _ _ Nothing) = False
hasFFT (Scalar e) = hasFFT e

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
factorOut e (Sum s _) | [(_,xx)] <- sum2pairs s = factorOut e xx
factorOut e (Product p _) = Map.lookup e p
factorOut _ _ = Nothing

allFactors :: Type a => Expression a -> Set.Set (Expression a)
allFactors (Product p _) = Map.keysSet p
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
                                      ES (Expression (Integrate e')) -> hasExpressionInFFT v e'
                                      EK (Expression Kx) -> False
                                      EK (Expression Ky) -> False
                                      EK (Expression Kz) -> False
                                      _ -> error "inexhaustive pattern in hasExpressionInFFT"
hasExpressionInFFT v (Var _ _ _ _ (Just e)) = hasExpressionInFFT v e
hasExpressionInFFT v (Sum s _) = or $ map (hasExpressionInFFT v . snd) (sum2pairs s)
hasExpressionInFFT v (Product p _) = or $ map (hasExpressionInFFT v . fst) (product2pairs p)
hasExpressionInFFT v (Sin e) = hasExpressionInFFT v e
hasExpressionInFFT v (Cos e) = hasExpressionInFFT v e
hasExpressionInFFT v (Log e) = hasExpressionInFFT v e
hasExpressionInFFT v (Exp e) = hasExpressionInFFT v e
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
scalarderive v (Cos e) = -sin e * scalarderive v e
scalarderive v (Sin e) = cos e * scalarderive v e
scalarderive v (Exp e) = exp e * scalarderive v e
scalarderive v (Log e) = scalarderive v e / e
scalarderive _ (Abs _) = error "I didn't think we'd need abs"
scalarderive _ (Signum _) = error "I didn't think we'd need signum"
scalarderive v (Expression e) = scalarderivativeHelper v e


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
derive v@(Var _ a b c _) dda0 (Var t aa bb cc (Just e0))
  | Just e <- safeCoerce e0 v,
    Just dda <- safeCoerce dda0 v =
    case isConstant $ derive v dda e of
      Just x -> toExpression x
      Nothing ->
        case isConstant $ derive v 1 e of
          Just x -> toExpression x * dda
          Nothing ->
            if dda == 1 || not (hasExpressionInFFT v e)
            then dda*(Var t ("d" ++ aa ++ "_by_d" ++ a) ("d" ++ bb ++ "_by_d" ++ b)
                      ("\\frac{\\partial "++cc ++"}{\\partial "++c++"}") $ Just $
                      (derive v 1 e))
            else derive v dda e
derive v@(Scalar (Var _ a b c _)) dda0 (Var t aa bb cc (Just e0))
  | Just e <- safeCoerce e0 v,
    Just dda <- safeCoerce dda0 v =
  case isConstant $ derive v dda e of
    Just x -> toExpression x
    Nothing ->
      case isConstant $ derive v 1 e of
        Just x -> toExpression x * dda
        Nothing ->
            if dda == 1 || not (hasExpressionInFFT v e)
            then dda*(Var t ("d" ++ aa ++ "_by_d" ++ a) ("d" ++ bb ++ "_by_d" ++ b)
                      ("\\frac{\\partial "++cc ++"}{\\partial "++c++"}") $ Just $
                      (derive v 1 e))
            else derive v dda e
derive v dda (Var _ _ _ _ (Just e)) = derive v dda e
derive _ _ (Var _ _ _ _ Nothing) = 0
derive v dda (Sum s _) = pairs2sum $ map dbythis $ sum2pairs s
  where dbythis (f,x) = (f, derive v dda x)
derive v dda (Product p i) = pairs2sum $ map dbythis $ product2pairs p
  where dbythis (x,n) = (1, derive v (Product p i*toExpression n*dda/x) x)
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

hasExprn :: Exprn -> Exprn -> Bool
hasExprn (ES a) (ES b) = hasexpression a b
hasExprn (ES a) (ER b) = hasexpression a b
hasExprn (ES a) (EK b) = hasexpression a b
hasExprn (ES a) (E3 b) = hasexpression a b
hasExprn (ER a) (ES b) = hasexpression a b
hasExprn (ER a) (ER b) = hasexpression a b
hasExprn (ER a) (EK b) = hasexpression a b
hasExprn (ER a) (E3 b) = hasexpression a b
hasExprn (EK a) (ES b) = hasexpression a b
hasExprn (EK a) (ER b) = hasexpression a b
hasExprn (EK a) (EK b) = hasexpression a b
hasExprn (EK a) (E3 b) = hasexpression a b
hasExprn (E3 a) (ES b) = hasexpression a b
hasExprn (E3 a) (ER b) = hasexpression a b
hasExprn (E3 a) (EK b) = hasexpression a b
hasExprn (E3 a) (E3 b) = hasexpression a b

countexpression :: (Type a, Type b) => Expression a -> Expression b -> Int
countexpression x e = snd $ subAndCount x (s_var "WeAreCounting") e

substitute :: (Type a, Type b) => Expression a -> Expression a -> Expression b -> Expression b
substitute x y e = fst $ subAndCount x y e

substituteE :: Type a => Expression a -> Expression a -> Exprn -> Exprn
substituteE a a' (ES b) = mkExprn $ substitute a a' b
substituteE a a' (ER b) = mkExprn $ substitute a a' b
substituteE a a' (EK b) = mkExprn $ substitute a a' b
substituteE a a' (E3 b) = mkExprn $ substitute a a' b

countAfterRemoval :: (Type a, Type b) => Expression a -> Expression b -> Int
countAfterRemoval v e = Set.size (varsetAfterRemoval v e)

countAfterRemovalE :: Type a => Exprn -> Expression a -> Int
countAfterRemovalE (EK a) b = countAfterRemoval a b
countAfterRemovalE (ER a) b = countAfterRemoval a b
countAfterRemovalE (ES a) b = countAfterRemoval a b
countAfterRemovalE (E3 a) b = countAfterRemoval a b

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
  | ES (Expression (Integrate e)) <- mkExprn v = varsetAfterRemoval x e
  | otherwise = Set.empty
varsetAfterRemoval x (Sum s _) = Set.unions (map (varsetAfterRemoval x . snd) (sum2pairs s))
varsetAfterRemoval x (Product p _) = Set.unions (map (varsetAfterRemoval x . fst) (product2pairs p))
varsetAfterRemoval x (Cos e) = varsetAfterRemoval x e
varsetAfterRemoval x (Sin e) = varsetAfterRemoval x e
varsetAfterRemoval x (Log e) = varsetAfterRemoval x e
varsetAfterRemoval x (Exp e) = varsetAfterRemoval x e
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
subAndCount x y (Scalar e) = (scalar e', n)
    where (e', n) = subAndCount x y e


newtype MkBetter a = MB (Maybe (Expression a -> MkBetter a, Expression a))

findRepeatedSubExpression :: Type a => Expression a -> MkBetter a
findRepeatedSubExpression everything = frse everything
  where frse x@(Sum s _) | Map.size s > 1,
                           MB (Just (better,_)) <- frs (sum2pairs s) = better x
          where frs ((_,y):ys) = case frse y of
                                   MB Nothing -> frs ys
                                   se -> se
                frs [] = MB Nothing
        frse x@(Product s _) | Map.size s > 1,
                               MB (Just (better,_)) <- frs (product2pairs s) = better x
          where frs ((y,n):ys) = case frse y of
                                   MB Nothing -> frs ys
                                   MB (Just (better,_)) -> better (y ** toExpression n)
                frs [] = MB Nothing
        frse x@(Cos e) | MB (Just (better, _)) <- frse e = better x
        frse x@(Sin e) | MB (Just (better, _)) <- frse e = better x
        frse x@(Log e) | MB (Just (better, _)) <- frse e = better x
        frse x@(Exp e) | MB (Just (better, _)) <- frse e = better x
        frse (Expression _) = MB Nothing
        frse (Var _ _ _ _ Nothing) = MB Nothing
        frse x@(Var _ _ _ _ (Just e)) | MB (Just (better, _)) <- frse e = better x
        frse (Scalar (Var _ _ _ _ Nothing)) = MB Nothing
        frse e = if mytimes > 1
                 then MB (Just (makebetter e mytimes, e))
                 else MB Nothing
           where mytimes = countexpression e everything
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
searchMonoid f x@(Cos e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Sin e) = f x `mappend` searchMonoid f e
searchMonoid f x@(Exp e) = f x `mappend` searchMonoid f e
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
        findInputsE (E3 ee) = findInputs ee

findInputs :: Type b => Expression b -> Set.Set Exprn
findInputs = searchMonoid helper
  where helper e | ER (Var _ _ _ _ Nothing) <- mkExprn e,
                   not (Set.member (mkExprn e) grid_description) = Set.insert (mkExprn e) grid_description
        helper e@(Var _ _ _ _ Nothing) = Set.singleton (mkExprn e)
        helper _ = Set.empty
        grid_description = Set.fromList [ES numx, ES numy, ES numz,
                                         E3 lat1, E3 lat2, E3 lat3]

dV, volume :: Type a => Expression a
dV = protect "dV" "\\Delta V" $ volume / numx / numy / numz
volume = protect "volume" "volume" $ lat1 `dot` (lat2 `cross` lat3)

numx, numy, numz :: Type a => Expression a
numx = s_var "Nx"
numy = s_var "Ny"
numz = s_var "Nz"

lat1, lat2, lat3 :: Expression ThreeVector
lat1 = Var CannotBeFreed "lat1" "lat1" "\\vec{a_1}" Nothing
lat2 = Var CannotBeFreed "lat2" "lat2" "\\vec{a_2}" Nothing
lat3 = Var CannotBeFreed "lat3" "lat3" "\\vec{a_3}" Nothing
\end{code}
