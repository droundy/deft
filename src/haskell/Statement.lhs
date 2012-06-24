The statement type allows us to string together a sequence of
definitions and assignments.

\begin{code}
{-# LANGUAGE PatternGuards #-}

module Statement ( Statement(..),
                   codeStatements,
                   latexStatements,
                   latexSimp,
                   simp2,
                   findToDo, findNamedSubexpression, findNiceScalar,
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

\begin{code}
numtotake :: Int
numtotake = 20000

subsq :: [a] -> [[a]]
subsq xs = -- map (:[]) xs ++
           rest xs
    where rest (y:ys@(_:_)) = map (:[y]) ys ++ rest ys
          rest _ = []

findNamedScalar :: Type b => Expression b -> Maybe Exprn
findNamedScalar e@(Expression _)
  | EK (Expression (FFT e')) <- mkExprn e = findNamedScalar e'
  | ER (Expression (IFFT e')) <- mkExprn e = findNamedScalar e'
  | ES (Expression (Integrate e')) <- mkExprn e = findNamedScalar e'
  | otherwise = Nothing
findNamedScalar (Sum s _) = case filter (/= Nothing) $ map sub $ sum2pairs s of
                                [] -> Nothing
                                dothis:_ -> dothis
    where sub (_,e) = findNamedScalar e
findNamedScalar (Product p _) = case filter (/= Nothing) $ map sub $ product2pairs p of
                                               [] -> Nothing
                                               dothis:_ -> dothis
    where sub (e, _) = findNamedScalar e
findNamedScalar (Cos e) = findNamedScalar e
findNamedScalar (Sin e) = findNamedScalar e
findNamedScalar (Log e) = findNamedScalar e
findNamedScalar (Exp e) = findNamedScalar e
findNamedScalar (Abs e) = findNamedScalar e
findNamedScalar (Signum e) = findNamedScalar e
findNamedScalar (Var t a b c (Just e)) =
  case findNamedScalar e of
    Just (ES e') -> if ES e' == mkExprn e
                    then Just $ ES $ Var t a b c (Just e')
                    else Just $ ES e'
    Nothing -> if amScalar e
               then Just $ mkExprn $ Var t a b c (Just e)
               else Nothing
    _ -> error "impossible case in findNamedScalar"
findNamedScalar (Var _ _ _ _ Nothing) = Nothing
findNamedScalar (Scalar e) = findNamedScalar e

findNiceScalar :: Type b => Expression b -> Maybe Exprn
findNiceScalar e@(Expression _)
  | EK (Expression (FFT e')) <- mkExprn e = findNiceScalar e'
  | ER (Expression (IFFT e')) <- mkExprn e = findNiceScalar e'
  | ES (Expression (Integrate e')) <- mkExprn e = if hasFFT e'
                                                  then findNiceScalar e'
                                                  else Just $ mkExprn e
  | otherwise = Nothing
findNiceScalar (Sum s _) = case filter (/= Nothing) $ map sub $ sum2pairs s of
                             [] -> Nothing
                             dothis:_ -> dothis
    where sub (_,e) = findNiceScalar e
findNiceScalar (Product p _) = case filter (/= Nothing) $ map sub $ product2pairs p of
                                 [] -> Nothing
                                 dothis:_ -> dothis
    where sub (e, _) = findNiceScalar e
findNiceScalar (Cos e) = findNiceScalar e
findNiceScalar (Sin e) = findNiceScalar e
findNiceScalar (Log e) = findNiceScalar e
findNiceScalar (Exp e) = findNiceScalar e
findNiceScalar (Abs e) = findNiceScalar e
findNiceScalar (Signum e) = findNiceScalar e
findNiceScalar (Var _ _ _ _ (Just e)) = findNiceScalar e
findNiceScalar (Var _ _ _ _ Nothing) = Nothing
findNiceScalar (Scalar (Var _ _ _ _ Nothing)) = Nothing
findNiceScalar (Scalar e) | not (hasFFT e) = Just $ ES e
                          | otherwise = Nothing

findNamedSubexpression :: Type b => Expression b -> Maybe Exprn
findNamedSubexpression e@(Expression _)
  | EK (Expression (FFT e')) <- mkExprn e = findNamedSubexpression e'
  | ER (Expression (IFFT e')) <- mkExprn e = findNamedSubexpression e'
  | ES (Expression (Integrate e')) <- mkExprn e = findNamedSubexpression e'
  | otherwise = Nothing
findNamedSubexpression (Sum s _) = case filter (/= Nothing) $ map sub $ sum2pairs s of
                                     [] -> Nothing
                                     dothis:_ -> dothis
    where sub (_,e) = findNamedSubexpression e
findNamedSubexpression (Product p _) = case filter (/= Nothing) $ map sub $ product2pairs p of
                                         [] -> Nothing
                                         dothis:_ -> dothis
    where sub (e, _) = findNamedSubexpression e
findNamedSubexpression (Cos e) = findNamedSubexpression e
findNamedSubexpression (Sin e) = findNamedSubexpression e
findNamedSubexpression (Log e) = findNamedSubexpression e
findNamedSubexpression (Exp e) = findNamedSubexpression e
findNamedSubexpression (Abs e) = findNamedSubexpression e
findNamedSubexpression (Signum e) = findNamedSubexpression e
findNamedSubexpression e@(Var _ _ _ _ (Just e'))
  | Just ee <- findNamedSubexpression e' = Just ee
  | otherwise = Just (mkExprn e)
  | otherwise = error "impossible case in findNamedSubexpression"
findNamedSubexpression (Var _ _ _ _ Nothing) = Nothing
findNamedSubexpression (Scalar e) = findNamedSubexpression e

findToDo :: (Type a, Type b) => Set.Set String -> Expression a -> Expression b -> Maybe Exprn
findToDo _ _ (Var _ _ _ _ Nothing) = Nothing
findToDo i _ e | Set.size i > 0 && not (Set.isSubsetOf i (varSet e)) = Nothing
               | countVars e == 0 = Nothing
--findToDo i everything e
--  | Set.isSubsetOf i (varSet e) && countVars e == 1 && not (hasFFT e) &&
--                 countAfterRemoval e everything + 1 < countVars everything = Just $ mkExprn e
findToDo i everything e@(Expression _)
    | EK (Expression (FFT e')) <- mkExprn e = findToDo i everything e'
    | ER (Expression (IFFT e')) <- mkExprn e = findToDo i everything e'
    | ES (Expression (Integrate e')) <- mkExprn e = if hasFFT e'
                                                    then findToDo i everything e'
                                                    else Just $ ES $ Expression (Integrate e')
    | otherwise = Nothing
findToDo _ _ (Sum _ i) | Set.size i < 2 = Nothing
findToDo _ everything (Sum s _) | todo:_ <- filter simplifiable subes = Just todo
    where subes = map (mkExprn . pairs2sum) $ take numtotake $ subsq $
                  take numtotake $ filter (\(_,e) -> countVars e > 0 && not (hasFFT e)) $ sum2pairs s
          simplifiable sube = countVarssube > 1 && countVarssube < 3 && ithelps sube
              where countVarssube = countVarsE sube
          oldnum = countVars everything
          ithelps e = countAfterRemovalE e everything + 1 < oldnum
findToDo i everything (Sum s _) = case filter (/= Nothing) $ map sub $ sum2pairs s of
                                [] -> Nothing
                                dothis:_ -> dothis
    where sub (_,e) = findToDo i everything e
findToDo _ _ (Product _ i) | Set.size i < 2 = Nothing
findToDo _ everything (Product p _) | todo:_ <- filter simplifiable subes = Just todo
    where subes = map (mkExprn . pairs2product) $ take numtotake $ subsq $
                  take numtotake $ filter (\(e,_) -> countVars e > 0 &&
                                                     kok (mkExprn e) &&
                                                     not (hasFFT e)) $ product2pairs p
          simplifiable sube = countVarssube > 1 && countVarssube < 3 && ithelps sube
              where countVarssube = countVarsE sube
          oldnum = countVars everything
          ithelps e = countAfterRemovalE e everything + 1 < oldnum
          zerodenom = iszero (product2denominator p)
          kok (EK kk) = if zerodenom -- kok is true if we don't have to worry about a divide by
                        then not (hasK kk) -- zero in the limit when when k == 0
                        else True
          kok _ = True
          iszero e = hasK e && setZero (EK kx) (setZero (EK ky) (setZero (EK kz) e)) == 0
findToDo i everything (Product p _) =
    if iszero (product2denominator p)
    then case filter notk $ filter (/= Nothing) $ map sub $ product2pairs p of
           [] -> Nothing
           dothis:_ -> dothis
    else case filter (/= Nothing) $ map sub $ product2pairs p of
           [] -> Nothing
           dothis:_ -> dothis
    where sub (e, _) = findToDo i everything e
          iszero e = hasK e && setZero (EK kx) (setZero (EK ky) (setZero (EK kz) e)) == 0
          notk (Just (EK e)) = not (hasK e)
          notk _ = True
findToDo i x (Cos e) = findToDo i x e
findToDo i x (Sin e) = findToDo i x e
findToDo i x (Log e) = findToDo i x e
findToDo i x (Exp e) = findToDo i x e
findToDo i x (Abs e) = findToDo i x e
findToDo i x (Signum e) = findToDo i x e
findToDo i x (Var t a b c (Just e)) =
  case findToDo i x e of
    Just e' -> if e' == mkExprn e
               then case e' of
                    ES e'' -> Just $ ES (Var t a b c (Just e''))
                    EK e'' -> Just $ EK (Var t a b c (Just e''))
                    ER e'' -> Just $ ER (Var t a b c (Just e''))
               else Just e'
    Nothing -> if amScalar e && hasFFT e
               then Just $ mkExprn $ Var t a b c (Just e)
               else Nothing
findToDo i everything (Scalar e) = findToDo i everything e


findFFTtodo :: (Type a, Type b) => Expression a -> Expression b -> Maybe Exprn
findFFTtodo everything e@(Expression _)
    | EK (Expression (FFT (Var _ _ _ _ Nothing))) <- mkExprn e = Just $ mkExprn e
    | EK (Expression (FFT e')) <- mkExprn e = findFFTtodo everything e'
    | ER (Expression (IFFT (Var _ _ _ _ Nothing))) <- mkExprn e = Just $ mkExprn e
    | ER (Expression (IFFT e')) <- mkExprn e = findFFTtodo everything e'
    | ES (Expression (Integrate e')) <- mkExprn e = findFFTtodo everything e'
    | otherwise = Nothing
findFFTtodo everything (Sum s _) = case filter (/= Nothing) $ map sub $ sum2pairs s of
                                [] -> Nothing
                                dothis:_ -> dothis
    where sub (_,e) = findFFTtodo everything e
findFFTtodo everything (Product p _) = case filter (/= Nothing) $ map sub $ product2pairs p of
                         [] -> Nothing
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
    Just e' -> if e' == mkExprn e
               then case e' of
                    ES e'' -> Just $ ES (Var t a b c (Just e''))
                    EK e'' -> Just $ EK (Var t a b c (Just e''))
                    ER e'' -> Just $ ER (Var t a b c (Just e''))
               else Just e'
    Nothing -> Nothing
findFFTtodo _ (Var _ _ _ _ Nothing) = Nothing
findFFTtodo everything (Scalar e) = findFFTtodo everything e


findFFTinputtodo :: (Type a, Type b) => Set.Set String -> Expression a -> Expression b -> Maybe Exprn
findFFTinputtodo i everything e@(Expression _)
    | EK (Expression (FFT e')) <- mkExprn e = if hasFFT e'
                                              then findFFTinputtodo i everything e'
                                              else Just $ ER e'
    | ER (Expression (IFFT e')) <- mkExprn e = if hasFFT e'
                                               then findFFTinputtodo i everything e'
                                               else Just $ EK e'
    | ES (Expression (Integrate e')) <- mkExprn e = if hasFFT e'
                                                    then findFFTinputtodo i everything e'
                                                    else Just $ ES $ Expression $ Integrate e'
    | otherwise = Nothing
findFFTinputtodo i everything (Product p _) = case filter (/= Nothing) $ map sub $ product2pairs p of
                                                [] -> Nothing
                                                dothis:_ -> dothis
    where sub (e, _) = findFFTinputtodo i everything e
findFFTinputtodo i x (Cos e) = findFFTinputtodo i x e
findFFTinputtodo i x (Sin e) = findFFTinputtodo i x e
findFFTinputtodo i x (Log e) = findFFTinputtodo i x e
findFFTinputtodo i x (Exp e) = findFFTinputtodo i x e
findFFTinputtodo i x (Abs e) = findFFTinputtodo i x e
findFFTinputtodo i x (Signum e) = findFFTinputtodo i x e
findFFTinputtodo i x (Var t a b c (Just e)) =
  case findFFTinputtodo i x e of
    Just e' -> if e' == mkExprn e
               then case e' of
                    ES e'' -> Just $ ES (Var t a b c (Just e''))
                    EK e'' -> Just $ EK (Var t a b c (Just e''))
                    ER e'' -> Just $ ER (Var t a b c (Just e''))
               else Just e'
    Nothing -> Nothing
findFFTinputtodo _ _ (Var _ _ _ _ Nothing) = Nothing
findFFTinputtodo i everything (Scalar e) = findFFTinputtodo i everything e
-- If possible, we want to only look for ffts that are in a sum
-- including the variable we most recently evaluated, since that's
-- where we're most likely to find more memory-saving improvements.
findFFTinputtodo i _ e | Set.size i > 0 && not (Set.isSubsetOf i (varSet e)) = Nothing
findFFTinputtodo i everything (Sum s _) = case filter (/= Nothing) $ map sub $ sum2pairs s of
                                            [] -> Nothing
                                            dothis:_ -> dothis
    where sub (_,e) = findFFTinputtodo i everything e


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
              Nothing ->
                case findNiceScalar e of
                  Just (ES s) -> scalarhelper (sts++[Initialize (ES v), Assign (ES v) (ES s)]) (n+1) (substitute s v e)
                    where v = Var CannotBeFreed ("s"++show n) ("s"++show n)
                                                ("s_{"++show n++"}") Nothing :: Expression Scalar
                  Nothing -> simp2helper Set.empty (n :: Int) sts e e
                  dothis -> error ("bad result in scalarhelper: " ++ show dothis)
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
                      Nothing -> handletodos ts

\end{code}
