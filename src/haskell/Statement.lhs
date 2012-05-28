The statement type allows us to string together a sequence of
definitions and assignments.

\begin{code}
{-# LANGUAGE GADTs, PatternGuards #-}

module Statement ( Statement(..),
                   ToDo(..),
                   codeStatement, 
                   codeStatements,
                   latexStatements,
                   latexSimp,
                   simp2,
                   findToDo,
                   countFFT,
                   checkDup,
                   peakMem,
                   reuseVar,
                   freeVectors)
    where

import Expression
import Data.List ( nubBy, partition, (\\) )
-- import qualified Data.Set as Set

data Statement = AssignR (Expression RealSpace) (Expression RealSpace)
               | AssignK (Expression KSpace) (Expression KSpace)
               | AssignS (Expression Scalar) (Expression Scalar)
               | InitializeR (Expression RealSpace)
               | InitializeK (Expression KSpace)
               | InitializeS (Expression Scalar)
               | FreeR (Expression RealSpace)
               | FreeK (Expression KSpace)

instance Show Statement where
  showsPrec = showsS
showsS :: Int -> Statement -> ShowS
showsS _ (AssignR x y) = showsPrec 0 x . showString " := " . showsPrec 0 y
showsS _ (AssignK x y) = showsPrec 0 x . showString " := " . showsPrec 0 y
showsS _ (AssignS x y) = showsPrec 0 x . showString " := " . showsPrec 0 y
showsS _ (InitializeR x) = showString "Initialize " . showsPrec 0 x
showsS _ (InitializeK x) = showString "Initialize " . showsPrec 0 x
showsS _ (InitializeS x) = showString "Initialize " . showsPrec 0 x
showsS _ (FreeR x) = showString "Free " . showsPrec 0 x
showsS _ (FreeK x) = showString "Free " . showsPrec 0 x

instance Code Statement where
  codePrec = codeS
  latexPrec = latexS
codeS :: Int -> Statement -> ShowS
codeS _ (AssignR x y) = showString (prefix y) . codePrec 0 x . showString " = " . codePrec 0 y . showString ";" . showString (postfix y)
codeS _ (AssignK x y) = showString (prefix y) . codePrec 0 x . showString " = " . codePrec 0 y . showString ";" . showString (postfix y)
codeS _ (AssignS x y) = showString (prefix y) . codePrec 0 x . showString " = " . codePrec 0 y . showString ";" . showString (postfix y)
codeS _ (InitializeR e) = showString (initialize e)
codeS _ (InitializeK e) = showString (initialize e)
codeS _ (InitializeS e) = showString (initialize e)
codeS _ (FreeR e) = showString (free e)
codeS _ (FreeK e) = showString (free e)

latexS :: Int -> Statement -> ShowS
latexS _ (AssignR x y) = latexPrec 0 x . showString " = " . latexPrec 0 (cleanvars y)
latexS _ (AssignK x y) = latexPrec 0 x . showString " = " . latexPrec 0 (cleanvars y)
latexS _ (AssignS x y) = latexPrec 0 x . showString " = " . latexPrec 0 (cleanvars y)
latexS _ (InitializeR e) = showString (initialize e)
latexS _ (InitializeK e) = showString (initialize e)
latexS _ (InitializeS e) = showString (initialize e)
latexS _ (FreeR e) = showString (free e)
latexS _ (FreeK e) = showString (free e)

latexStatements :: [Statement] -> String
latexStatements x = unlines $ map (\e -> "\n\\begin{dmath}\n" ++ latex e ++ "\n\\end{dmath}") x

latexSimp :: (Type a) => Expression a -> String
latexSimp e = "\\documentclass{article}\n\\usepackage{amsmath}\n\\usepackage{breqn}\n\\begin{document}\n\n" ++
              latexStatements ( freeVectors $ f sts) ++ "\n\\begin{dmath}\n" ++ latex (cleanvars e') ++ "\n\\end{dmath}" ++
              "\n\\end{document}"
    where f (x:xs) = x:(f xs)
          f [] = []
          (sts,e') = simp2 e

codeStatements :: [Statement] -> String
codeStatements (InitializeS v : AssignS v' e : ss)
  | v == v' = "\tdouble " ++ codeStatementHelper v " = " e ++ "\n" ++
              codeStatements ss
codeStatements (s:ss) = codeStatement s ++ "\n" ++ codeStatements ss
codeStatements [] = ""

codeStatement :: Statement -> String
codeStatement (AssignR x y) = "\t" ++ codeStatementHelper x " = " y
codeStatement (AssignK x y) = "\t" ++ codeStatementHelper x " = " y
codeStatement (AssignS x y) = "\t" ++ codeStatementHelper x " = " y
codeStatement (InitializeR e) = "\t" ++ initialize e
codeStatement (InitializeK e) = "\t" ++ initialize e
codeStatement (InitializeS e) = "\t" ++ initialize e
codeStatement (FreeR e) = "\t" ++ free e
codeStatement (FreeK e) = "\t" ++ free e

substituteS :: Type a => Expression a -> Expression a -> Statement -> Statement
substituteS x y (AssignR s e) = AssignR s (substitute x y e)
substituteS x y (AssignK s e) = AssignK s (substitute x y e)
substituteS x y (AssignS s e) = AssignS s (substitute x y e)
substituteS x y (InitializeR e) = InitializeR (substitute x y e)
substituteS x y (InitializeK e) = InitializeK (substitute x y e)
substituteS x y (InitializeS e) = InitializeS (substitute x y e)
substituteS x y (FreeR e) = FreeR (substitute x y e)
substituteS x y (FreeK e) = FreeK (substitute x y e)

countFFT :: [Statement] -> Int
countFFT = sum . map helper
    where helper (AssignR _ (Expression (IFFT _))) = 1
          helper (AssignK _ (Expression (FFT _))) = 1
          helper (AssignR _ (Var _ _ _ _ (Just (Expression (IFFT _))))) = 1
          helper (AssignK _ (Var _ _ _ _ (Just (Expression (FFT _))))) = 1
          helper _ = 0

peakMem :: [Statement] -> Int
peakMem = maximum . (helper 0) . freeVectors
    where helper n (x:xs) | (InitializeR _) <- x = (n+1) : (helper (n+1) xs)
                          | (InitializeK _) <- x = (n+1) : (helper (n+1) xs)
                          | (FreeR _) <- x = (n-1) : (helper (n-1) xs)
                          | (FreeK _) <- x = (n-1) : (helper (n-1) xs)
                          | otherwise = n : helper n xs
          helper n [] = [n]

checkDup :: [Statement] -> [Statement]
checkDup = nubBy sameInit
    where sameInit (InitializeR x) (InitializeR y) = x == y
          sameInit (InitializeK x) (InitializeK y) = x == y
          sameInit _ _ = False

hasE :: Type a => Expression a -> [Statement] -> Bool
hasE _ [] = False
hasE e (AssignK _ a:_) | hasexpression e a = True
hasE e (AssignR _ a:_) | hasexpression e a = True
hasE e (AssignS _ a:_) | hasexpression e a = True
hasE e (FreeK a:_) | Same <- compareExpressions a e = True
hasE e (FreeR a:_) | Same <- compareExpressions a e = True
hasE e (_:rest) = hasE e rest

-- FIXME: This freeVectors uses O(n^2) calls to hasexpression, which
-- is worse than before, and is worse than we should be able to do!
-- :( On the plus side, it now preserves variable names, and reuseVar
-- works again.
freeVectors :: [Statement] -> [Statement]
freeVectors = freeHelper [] []
    where freeHelper :: [Expression RealSpace] -> [Expression KSpace] -> [Statement] -> [Statement]
          freeHelper rs ks (xnn@(FreeR v):xs) = [xnn] ++ freeHelper (rs \\ [v]) ks xs
          freeHelper rs ks (xnn@(FreeK v):xs) = [xnn] ++ freeHelper rs (ks \\ [v]) xs
          freeHelper rs ks (xnn@(InitializeR v@(Var IsTemp _ _ _ Nothing)):xs) = [xnn] ++ freeHelper (v:rs) ks xs
          freeHelper _ _ (InitializeR v:_) = error ("crazy InitializeR: " ++ show v)
          freeHelper rs ks (xnn@(InitializeK v@(Var IsTemp _ _ _ Nothing)):xs) = [xnn] ++ freeHelper rs (v:ks) xs
          freeHelper _ _ (InitializeK v:_) = error ("crazy InitializeK: " ++ show v)
          freeHelper rs ks (xnn@(InitializeS _):xs) = [xnn] ++ freeHelper rs ks xs
          freeHelper rs ks (xnn@(AssignR _ _):xs) = xnn : freeme ++ freeHelper rs' ks' xs
            where (rs', rsfree) = partition (`hasE` xs) rs
                  (ks', ksfree) = partition (`hasE` xs) ks
                  freeme = map FreeR rsfree ++ map FreeK ksfree
          freeHelper rs ks (xnn@(AssignK _ _):xs) = xnn : freeme ++ freeHelper rs' ks' xs
            where (rs', rsfree) = partition (`hasE` xs) rs
                  (ks', ksfree) = partition (`hasE` xs) ks
                  freeme = map FreeR rsfree ++ map FreeK ksfree
          freeHelper rs ks (xnn@(AssignS _ _):xs) = xnn : freeme ++ freeHelper rs' ks' xs
            where (rs', rsfree) = partition (`hasE` xs) rs
                  (ks', ksfree) = partition (`hasE` xs) ks
                  freeme = map FreeR rsfree ++ map FreeK ksfree
          freeHelper _ _ [] = []

reuseVar :: [Statement] -> [Statement]
reuseVar ((InitializeR iivar@(Var IsTemp _ _ _ Nothing)):(AssignR n e):(FreeR ffvar@(Var IsTemp _ _ _ Nothing)):xs)
    | iivar == n = (AssignR ffvar e) : reuseVar (map (substituteS iivar ffvar) xs)
    | otherwise = error "RS initialize error: "
reuseVar ((InitializeK iivar@(Var IsTemp _ _ _ Nothing)):(AssignK n e):(FreeK ffvar@(Var IsTemp _ _ _ Nothing)):xs)
    | iivar == n = (AssignK ffvar e) : reuseVar (map (substituteS iivar ffvar) xs)
    | otherwise = error "initialize error"
reuseVar (x:xs) = x : (reuseVar xs)
reuseVar [] = []

\end{code}

\begin{code}

data ToDo = DoK (Expression KSpace) | DoR (Expression RealSpace) | DoS (Expression Scalar) | DoNothing
            deriving (Eq, Show)

numtotake :: Int
numtotake = 20

subsq :: [a] -> [[a]]
subsq xs = map (:[]) xs ++ rest xs
    where rest (y:ys@(_:_)) = map (:[y]) ys ++ rest ys
          rest _ = []

findToDo :: (Type a, Type b) => Expression a -> Expression b -> ToDo
findToDo everything (Expression e)
    | Same <- isKSpace (Expression e), FFT e' <- e = findToDo everything e'
    | Same <- isRealSpace (Expression e), IFFT e' <- e = findToDo everything e'
    | Same <- isScalar (Expression e), Integrate e' <- e = case findToDo everything e' of
                                                             DoNothing -> if not (hasFFT e')
                                                                          then DoS $ Expression e
                                                                          else DoNothing
                                                             dothis -> dothis
    | otherwise = if hasFFT (Expression e)
                  then error ("FFT not detected in findToDo: " ++ show e)
                  else DoNothing -- error ("Missed Expression type in findToDo: " ++ show e)
-- findToDo _ (Sum i _) | Set.size i < 2 = DoNothing
findToDo everything (Sum s i)
    | Same <- isRealSpace (Sum s i), todo:_ <- filter simplifiable subes = DoR todo
    | Same <- isKSpace (Sum s i), todo:_ <- filter simplifiable subes = DoK todo
    where acceptables = take numtotake $ filter (\(_,e) -> {- not (hasK e) && -} countVars e > 0 && not (hasFFT e)) $ sum2pairs s
          subes = map pairs2sum $ take numtotake $ subsq acceptables
          simplifiable sube = countVarssube > 1 && countVarssube < 3 && ithelps sube {- && (countexpression sube everything) > 1 -}
              where countVarssube = countVars sube
          oldnum = countVars everything
          ithelps :: Type a => Expression a -> Bool
          ithelps e | Same <- isRealSpace e = countVars (substitute e (r_var "rtempWhatIfISubstituteThis") everything) < oldnum
                    | Same <- isKSpace e = countVars (substitute e (k_var "ktempWhatIfISubstituteThis") everything) < oldnum
                    | otherwise = False
findToDo everything (Sum s _) = case filter (/= DoNothing) $ map sub $ sum2pairs s of
                                [] -> DoNothing
                                dothis:_ -> dothis
    where sub (_,e) = findToDo everything e
findToDo _ (Product p _) | [] <- product2pairs p = DoNothing
-- findToDo _ (Product i _) | False && Set.size i < 2 = DoNothing
findToDo everything (Product s i)
    | Same <- isRealSpace (Product s i), todo:_ <- filter simplifiable subes = DoR todo
    | Same <- isKSpace (Product s i), todo:_ <- filter simplifiable subes = DoK todo
    where acceptables = take numtotake $ filter (\(e,_) -> not (hasK e) && countVars e > 0 && not (hasFFT e)) $ product2pairs s
          subes = map pairs2product $ take numtotake $ subsq acceptables
          simplifiable sube = countVarssube > 1 && countVarssube < 3 && ithelps sube {- && (countexpression sube everything) > 1 -}
              where countVarssube = countVars sube
          oldnum = countVars everything
          ithelps :: Type a => Expression a -> Bool
          ithelps e | Same <- isRealSpace e = countVars (substitute e (r_var "rtempWhatIfISubstituteThis") everything) < oldnum
                    | Same <- isKSpace e = countVars (substitute e (k_var "ktempWhatIfISubstituteThis") everything) < oldnum
                    | otherwise = False
findToDo everything (Product p _) =
    if iszero (product2denominator p)
    then case filter notk $ filter (/= DoNothing) $ map sub $ product2pairs p of
           [] -> DoNothing
           dothis:_ -> dothis
    else case filter (/= DoNothing) $ map sub $ product2pairs p of
           [] -> DoNothing
           dothis:_ -> dothis
    where sub (e, _) = findToDo everything e
          iszero e = hasK e && setZero kx (setZero ky (setZero kz e)) == 0
          notk (DoK e) = not (hasK e)
          notk _ = True
findToDo x (Cos e) = findToDo x e
findToDo x (Sin e) = findToDo x e
findToDo x (Log e) = findToDo x e
findToDo x (Exp e) = findToDo x e
findToDo x (Abs e) = findToDo x e
findToDo x (Signum e) = findToDo x e
findToDo x (Var t a b c (Just e)) =
  case findToDo x e of
    DoR e' -> case compareExpressions e e' of
                Same -> DoR $ Var t a b c (Just e)
                Different -> DoR e'
    DoK e' -> case compareExpressions e e' of
                Same -> DoK $ Var t a b c (Just e)
                Different -> DoK e'
    DoS e' -> case compareExpressions e e' of
                Same -> DoS $ Var t a b c (Just e)
                Different -> DoS e'
    DoNothing -> case isScalar e of
                 Different -> DoNothing
                 Same -> if hasFFT e then DoNothing else DoS $ Var t a b c (Just e)
--findToDo _ (Var _ _ _ (Just e)) = findToDo e
findToDo _ (Var _ _ _ _ Nothing) = DoNothing
findToDo everything (Scalar e) = findToDo everything e


findFFTtodo :: (Type a, Type b) => Expression a -> Expression b -> ToDo
findFFTtodo everything (Expression e)
    | Same <- isKSpace (Expression e), FFT (Var _ _ _ _ Nothing) <- e = DoK (Expression e)
    | Same <- isKSpace (Expression e), FFT e' <- e = case findFFTtodo everything e' of
                                                       DoNothing -> DoR $ e'
                                                       dothis -> dothis
    | Same <- isRealSpace (Expression e), IFFT (Var _ _ _ _ Nothing) <- e = DoR (Expression e)
    | Same <- isRealSpace (Expression e), IFFT e' <- e = case findFFTtodo everything e' of
                                                           DoNothing -> DoK $ e'
                                                           dothis -> dothis
    | Same <- isScalar (Expression e), Integrate e' <- e = case findFFTtodo everything e' of
                                                             DoNothing -> DoS $ Expression e
                                                             dothis -> dothis
    | otherwise = if hasFFT (Expression e)
                  then error ("FFT not detected in : " ++ show e)
                  else DoNothing -- error ("Missed Expression type in findFFTtodo: " ++ show e)
findFFTtodo everything (Sum s i) = case filter (/= DoNothing) $ map sub $ sum2pairs s of
                                [] -> if hasFFT (Sum s i)
                                      then error ("FFT not detected in : " ++ show (Sum s i))
                                      else DoNothing
                                dothis:_ -> dothis
    where sub (_,e) = findFFTtodo everything e
findFFTtodo everything (Product p i) = case filter (/= DoNothing) $ map sub $ product2pairs p of
                         [] -> if hasFFT (Product p i)
                               then error ("FFT not detected in : " ++ show (Product p i))
                               else DoNothing
                         dothis:_ -> dothis
    where sub (e, _) = findFFTtodo everything e
findFFTtodo x (Cos e) = findFFTtodo x e
findFFTtodo x (Sin e) = findFFTtodo x e
findFFTtodo x (Log e) = findFFTtodo x e
findFFTtodo x (Exp e) = findFFTtodo x e
findFFTtodo x (Abs e) = findFFTtodo x e
findFFTtodo x (Signum e) = findFFTtodo x e
findFFTtodo x (Var t a b c (Just e)) =
  case findFFTtodo x e of
    DoR e' -> case compareExpressions e e' of
                Same -> DoR $ Var t a b c (Just e)
                Different -> DoR e'
    DoK e' -> case compareExpressions e e' of
                Same -> DoK $ Var t a b c (Just e)
                Different -> DoK e'
    DoS e' -> case compareExpressions e e' of
                Same -> DoS $ Var t a b c (Just e)
                Different -> DoS e'
    DoNothing -> DoNothing
--findFFTtodo _ (Var _ _ _ (Just e)) = findFFTtodo e
findFFTtodo _ (Var _ _ _ _ Nothing) = DoNothing
findFFTtodo everything (Scalar e) = findFFTtodo everything e


simp2 :: Type a => Expression a -> ([Statement], Expression a)
simp2 = simp2helper (0 :: Int) [] -- . cleanvars
    where simp2helper n sts e = case findToDo e e of
                                  DoK ke -> simp2helper (n+1) (sts++[inike, setke]) e'
                                      where v = case ke of --Var _ x@('r':_) t _ -> (x, t) Why we do this ???
                                                           Var _ xi x t _ -> Var IsTemp xi x t Nothing
                                                           _ -> Var IsTemp ("ktemp_" ++ show n++"[i]")
                                                                    ("ktemp_"++show n) ("ktemp_{" ++ show n ++ "}") Nothing
                                            inike = InitializeK v
                                            setke = AssignK v ke
                                            e'  = substitute ke v e
                                  DoR re -> simp2helper (n+1) (sts++[inire, setre]) e'
                                      where v = case re of --Var _ x@('k':_) t _ -> (x,t)
                                                           Var _ xi x t _ -> Var IsTemp xi x t Nothing
                                                           _ -> Var IsTemp ("rtemp_" ++ show n++"[i]")
                                                                    ("rtemp_"++show n) ("rtemp_{" ++ show n ++ "}") Nothing
                                            inire = InitializeR v
                                            setre = AssignR v re
                                            e'  = substitute re v e
                                  DoS re -> simp2helper (n+1) (sts++[inire, setre]) e'
                                      where v = case re of Var _ _ x t _ -> Var CannotBeFreed x x t Nothing
                                                           _ -> Var CannotBeFreed ("s" ++ show n) ("s" ++ show n) ("s_{" ++ show n ++ "}") Nothing
                                            inire = InitializeS v
                                            setre = AssignS v re
                                            e'  = substitute re v e
                                  DoNothing -> case findFFTtodo e e of
                                                 DoK ke -> simp2helper (n+1) (sts++[inike, setke]) e'
                                                     where v = case ke of
                                                                 Var _ xi x t _ -> Var IsTemp xi x t Nothing
                                                                 _ -> Var IsTemp ("ktemp_" ++ show n++"[i]")
                                                                          ("ktemp_"++show n) ("ktemp_{" ++ show n ++ "}") Nothing
                                                           inike = InitializeK v
                                                           setke = AssignK v ke
                                                           e'  = substitute ke v e
                                                 DoR re -> simp2helper (n+1) (sts++[inire, setre]) e'
                                                     where v = case re of
                                                                 Var _ xi x t _ -> Var IsTemp xi x t Nothing
                                                                 _ -> Var IsTemp ("rtemp_" ++ show n++"[i]")
                                                                      ("rtemp_"++show n) ("rtemp_{" ++ show n ++ "}") Nothing
                                                           inire = InitializeR v
                                                           setre = AssignR v re
                                                           e'  = substitute re v e
                                                 DoS re -> simp2helper (n+1) (sts++[inire, setre]) e'
                                                     where v = case re of Var _ _ x t _ -> Var CannotBeFreed x x t Nothing
                                                                          _ -> Var CannotBeFreed ("s" ++ show n) ("s" ++ show n) ("s_{" ++ show n ++ "}") Nothing
                                                           inire = InitializeS v
                                                           setre = AssignS v re
                                                           e'  = substitute re v e
                                                 DoNothing -> (sts, e)

\end{code}
