The statement type allows us to string together a sequence of
definitions and assignments.

\begin{code}
{-# LANGUAGE PatternGuards #-}

module Statement ( Statement(..),
                   codeStatements,
                   newcodeStatements,
                   latexStatements,
                   latexSimp,
                   simp2,
                   findToDo, findNamedSubexpression,
                   findNamedScalar, findFFTtodo, findFFTinputtodo,
                   countFFT,
                   checkDup,
                   peakMem,
                   reuseVar,
                   freeVectors)
    where

import Expression
import Data.List ( nubBy, partition, (\\) )
import qualified Data.Set as Set

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
  latexPrec _ (Assign x y) = latexPrec 0 x . showString " = " . latexPrec 0 (cleanvarsE y)
  latexPrec _ (Initialize e) = showString (initializeE e)
  latexPrec _ (Free e) = showString (freeE e)

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


numtotake :: Int
numtotake = 20000

subsq :: [a] -> [[a]]
subsq xs = -- map (:[]) xs ++
           rest xs
    where rest (y:ys@(_:_)) = map (:[y]) ys ++ rest ys
          rest _ = []

findNamedScalar :: Type b => Expression b -> Maybe Exprn
findNamedScalar = searchExpressionDepthFirst Set.empty helper
  where helper e@(Var _ _ _ _ (Just _)) | ES _ <- mkExprn e = Just $ mkExprn e
        helper _ = Nothing

findNamedSubexpression :: Type b => Expression b -> Maybe Exprn
findNamedSubexpression = searchExpressionDepthFirst Set.empty helper
  where helper e@(Var _ _ _ _ (Just _)) = Just $ mkExprn e
        helper _ = Nothing

findToDo :: (Type a, Type b) => Set.Set String -> Expression a -> Expression b -> Maybe Exprn
findToDo i everything = searchExpression i helper
  where helper (Sum _ ii) | Set.size ii < 2 = Nothing
        helper (Product _ ii) | Set.size ii < 2 = Nothing
        helper (Sum s _) | todo:_ <- filter simplifiable subes = Just todo
          where subes = map (mkExprn . pairs2sum) $ take numtotake $ subsq $
                        take numtotake $ filter (\(_,e) -> countVars e > 0 && not (hasFFT e)) $ sum2pairs s
                simplifiable sube = countVarssube > 1 && countVarssube < 3 && ithelps sube
                  where countVarssube = countVarsE sube
                oldnum = countVars everything
                ithelps e = countAfterRemovalE e everything + 1 < oldnum
        helper (Product p _) | todo:_ <- filter simplifiable subes = Just todo
          where subes = map (mkExprn . pairs2product) $ take numtotake $ subsq $
                        take numtotake $ filter (\(e,_) -> countVars e > 0 &&
                                                           not (hasFFT e)) $ product2pairs p
                simplifiable sube = countVarssube > 1 && countVarssube < 3 && ithelps sube
                  where countVarssube = countVarsE sube
                oldnum = countVars everything
                ithelps e = countAfterRemovalE e everything + 1 < oldnum
        helper _ = Nothing

findFFTtodo :: (Type a, Type b) => Expression a -> Expression b -> Maybe Exprn
findFFTtodo _ = searchExpression Set.empty helper
  where helper e@(Expression _)
          | EK (Expression (FFT (Var _ _ _ _ Nothing))) <- mkExprn e = Just $ mkExprn e
          | ER (Expression (IFFT (Var _ _ _ _ Nothing))) <- mkExprn e = Just $ mkExprn e
        helper _ = Nothing

findFFTinputtodo :: (Type a, Type b) => Set.Set String -> Expression a -> Expression b -> Maybe Exprn
findFFTinputtodo i _ = searchExpressionDepthFirst i helper
  where helper e@(Expression _)
          | EK (Expression (FFT e')) <- mkExprn e, not (hasFFT e') = Just $ ER e'
          | ER (Expression (IFFT e')) <- mkExprn e, not (hasFFT e') = Just $ EK e'
          | ES (Expression (Integrate e')) <- mkExprn e, not (hasFFT e') = Just $ mkExprn e
        helper _ = Nothing

simp2 :: Type a => Expression a -> ([Statement], Expression a)
simp2 eee = case scalarhelper [] 0 eee of
            (a,_,b) -> (a,b)
    where -- First, we want to evaluate any purely scalar expressions
          -- that we can, just to simplify things!
          scalarhelper :: Type a => [Statement] -> Int -> Expression a -> ([Statement], Int, Expression a)
          scalarhelper sts n e =
            case findNamedScalar e of
              Just (ES s@(Var _ _ x t (Just e'))) ->
                case simp2helper Set.empty n [] e e' of
                  ([],_,_) -> scalarhelper (sts++[Initialize (ES v), Assign (ES v) (ES s)]) n (substitute s v e)
                              where v = Var CannotBeFreed x x t Nothing :: Expression Scalar
                  (sts',n',e'') -> scalarhelper (sts++sts') n' e''
              Nothing -> simp2helper Set.empty (n :: Int) sts e e
              _ -> error "bad result in scalarhelper"
          -- Then we go looking for memory to save or ffts to evaluate...
          simp2helper :: (Type a, Type b) => Set.Set String -> Int -> [Statement] -> Expression a
                         -> Expression b
                         -> ([Statement], Int, Expression a)
          simp2helper i n sts everything e =
            if Set.size i == 0
            then handletodos [findToDo i everything e, findFFTtodo everything e, findFFTinputtodo i everything e]
            else handletodos [findToDo i everything e, findFFTtodo everything e,
                             findFFTinputtodo i everything e, findFFTinputtodo Set.empty everything e]
            where handletodos [] = (sts, n, everything)
                  handletodos (todo:ts) =
                    case todo of
                      Just (EK ke) -> simp2helper (varSet v) (n+1)
                                      (sts++[Initialize (EK v), Assign (EK v) (EK ke)])
                                      (substitute ke v everything) (substitute ke v e)
                        where v :: Expression KSpace
                              v = case ke of
                                Var _ xi x t _ -> Var IsTemp xi x t Nothing
                                _ -> Var IsTemp ("ktemp" ++ show n++"[i]")
                                                ("ktemp"++show n)
                                                ("\\tilde{f}_{" ++ show n ++ "}") Nothing
                      Just (ER re) -> simp2helper (varSet v) (n+1)
                                      (sts++[Initialize (ER v), Assign (ER v) (ER re)])
                                      (substitute re v everything) (substitute re v e)
                        where v :: Expression RealSpace
                              v = case re of
                                Var _ xi x t _ -> Var IsTemp xi x t Nothing
                                _ -> Var IsTemp ("rtemp" ++ show n++"[i]")
                                                ("rtemp"++show n)
                                                ("f_{" ++ show n ++ "}") Nothing
                      Just (ES se) -> simp2helper (varSet v) (n+1)
                                      (sts++[Initialize (ES v), Assign (ES v) (ES se)])
                                      (substitute se v everything) (substitute se v e)
                        where v :: Expression Scalar
                              v = case se of
                                Var _ xi x t _ -> Var IsTemp xi x t Nothing
                                _ -> Var CannotBeFreed ("s" ++ show n)
                                                       ("s"++show n)
                                                       ("s_{" ++ show n ++ "}") Nothing
                      Just (E3 ve) -> simp2helper (varSet v) (n+1)
                                      (sts++[Initialize (E3 v), Assign (E3 v) (E3 ve)])
                                      (substitute ve v everything) (substitute ve v e)
                        where v :: Expression ThreeVector
                              v = case ve of
                                Var _ xi x t _ -> Var IsTemp xi x t Nothing
                                _ -> Var CannotBeFreed ("v" ++ show n)
                                                       ("v"++show n)
                                                       ("\\vec{v_{" ++ show n ++ "}}") Nothing
                      Nothing -> handletodos ts

\end{code}
