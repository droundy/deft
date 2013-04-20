The statement type allows us to string together a sequence of
definitions and assignments.

\begin{code}
{-# LANGUAGE PatternGuards #-}

module Statement ( Statement(..),
                   codeStatements,
                   newcodeStatements,
                   latexStatements,
                   countFFT,
                   checkDup,
                   peakMem,
                   reuseVar,
                   freeVectors)
    where

import Expression
import Data.List ( nubBy, partition, (\\) )

data Statement = Assign Exprn Exprn
               | Initialize Exprn
               | Free Exprn

instance Show Statement where
  showsPrec _ (Assign x y) = showsPrec 0 x . showString " := " . showsPrec 0 y
  showsPrec _ (Initialize x) = showString "Initialize " . showsPrec 0 x
  showsPrec _ (Free x) = showString "Free " . showsPrec 0 x

instance Code Statement where
  codePrec _ (Assign x y) = showString ("\t" ++ codeStatementE x " = " y)
  codePrec _ (Initialize e) = showString ("\t" ++ initializeE e)
  codePrec _ (Free e) = showString ("\t" ++ freeE e)
  newcodePrec _ (Assign x y) = showString ("\t" ++ newcodeStatementE x " = " y)
  newcodePrec _ (Initialize e) = showString ("\t" ++ newinitializeE e)
  newcodePrec _ (Free e) = showString ("\t" ++ newfreeE e)
  latexPrec _ (Assign x y) = latexPrec 0 x . showString " = " . mapExprn (latexPrec 0 . cleanvars) y
  latexPrec _ (Initialize e) = showString (initializeE e)
  latexPrec _ (Free e) = showString (freeE e)

latexStatements :: [Statement] -> String
latexStatements x = unlines $ map (\e -> "\n\\begin{dmath}\n" ++ latex e ++ "\n\\end{dmath}") x

codeStatements :: [Statement] -> String
codeStatements (Initialize v : Assign v' (ES e) : ss)
  | v == v' = "\tdouble " ++ code (Assign v (ES e)) ++ "\n" ++ codeStatements ss
codeStatements (s:ss) = code s ++ "\n" ++ codeStatements ss
codeStatements [] = ""

newcodeStatements :: [Statement] -> String
newcodeStatements (Initialize v : Assign v' (ES e) : ss)
  | v == v' = "\tdouble " ++ newcode (Assign v (ES e)) ++ "\n" ++ newcodeStatements ss
newcodeStatements (s:ss) = newcode s ++ "\n" ++ newcodeStatements ss
newcodeStatements [] = ""

substituteS :: Type a => Expression a -> Expression a -> Statement -> Statement
substituteS x y (Assign s e) = Assign s (substituteE x y e)
substituteS x y (Initialize e) = Initialize (substituteE x y e)
substituteS x y (Free e) = Free (substituteE x y e)

countFFT :: [Statement] -> Int
countFFT = sum . map helper
    where helper (Assign _ (ER (Expression (IFFT _)))) = 1
          helper (Assign _ (EK (Expression (FFT _)))) = 1
          helper (Assign _ (ER (Var _ _ _ _ (Just (Expression (IFFT _)))))) = 1
          helper (Assign _ (EK (Var _ _ _ _ (Just (Expression (FFT _)))))) = 1
          helper _ = 0

peakMem :: [Statement] -> Int
peakMem = maximum . (helper 0) . freeVectors
    where helper n (x:xs) | (Initialize (EK _)) <- x = (n+1) : (helper (n+1) xs)
                          | (Initialize (ER _)) <- x = (n+1) : (helper (n+1) xs)
                          | (Free (ER _)) <- x = (n-1) : (helper (n-1) xs)
                          | (Free (EK _)) <- x = (n-1) : (helper (n-1) xs)
                          | otherwise = n : helper n xs
          helper n [] = [n]

checkDup :: [Statement] -> [Statement]
checkDup = nubBy sameInit
    where sameInit (Initialize x) (Initialize y) = x == y
          sameInit _ _ = False

hasE :: Exprn -> [Statement] -> Bool
hasE _ [] = False
hasE e (Assign _ a:_) | hasExprn e a = True
hasE e (Free a:_) | e == a = True
hasE e (_:rest) = hasE e rest

-- FIXME: This freeVectors uses O(n^2) calls to hasexpression, which
-- is worse than before, and is worse than we should be able to do!
-- :( On the plus side, it now preserves variable names, and reuseVar
-- works again.
freeVectors :: [Statement] -> [Statement]
freeVectors = freeHelper []
    where freeHelper :: [Exprn] -> [Statement] -> [Statement]
          freeHelper vs (xnn@(Free v):xs) = [xnn] ++ freeHelper (vs \\ [v]) xs
          freeHelper vs (xnn@(Initialize v):xs) | istemp v = [xnn] ++ freeHelper (v:vs) xs
                                                 | otherwise = [xnn] ++ freeHelper vs xs
          freeHelper vs (xnn@(Assign _ _):xs) = xnn : freeme ++ freeHelper vs' xs
            where (vs', vsfree) = partition (`hasE` xs) vs
                  freeme = map Free vsfree
          freeHelper _ [] = []
          istemp (ES (Var IsTemp _ _ _ Nothing)) = True
          istemp (EK (Var IsTemp _ _ _ Nothing)) = True
          istemp (ER (Var IsTemp _ _ _ Nothing)) = True
          istemp _ = False

reuseVar :: [Statement] -> [Statement]
reuseVar ((Initialize (ER iivar@(Var IsTemp _ _ _ Nothing))) :
          (Assign n e):(Free (ER ffvar@(Var IsTemp _ _ _ Nothing))) : xs)
    | ER iivar == n = (Assign (ER ffvar) e) : reuseVar (map (substituteS iivar ffvar) xs)
    | otherwise = error "RS initialize error: "
reuseVar ((Initialize (EK iivar@(Var IsTemp _ _ _ Nothing))) :
          (Assign n e):(Free (EK ffvar@(Var IsTemp _ _ _ Nothing))) : xs)
    | EK iivar == n = (Assign (EK ffvar) e) : reuseVar (map (substituteS iivar ffvar) xs)
    | otherwise = error "KS initialize error: "
reuseVar (x:xs) = x : (reuseVar xs)
reuseVar [] = []
\end{code}
