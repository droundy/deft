{-# LANGUAGE PatternGuards #-}

module Optimize ( optimize,
                  findToDo, findNamedSubexpression,
                  findNamedScalar, findFFTtodo, findFFTinputtodo )
    where

import Expression
import Statement
import qualified Data.Set as Set

numtotake :: Int
numtotake = 20000

subsq :: [a] -> [[a]]
subsq xs = -- map (:[]) xs ++
           rest xs
    where rest (y:ys@(_:_)) = map (:[y]) ys ++ rest ys
          rest _ = []

findNamedScalar :: Type b => Expression b -> Maybe Exprn
findNamedScalar xxxx = mconcat [find volume,
                                find mydV,
                                find mydr,
                                searchExpressionDepthFirst Set.empty helper xxxx]
  where helper e@(Var _ _ _ _ (Just _)) | ES _ <- mkExprn e = Just $ mkExprn e
        helper _ = Nothing
        find :: Expression Scalar -> Maybe Exprn
        find x = if hasexpression x xxxx then Just $ mkExprn x else Nothing
        mydV = substitute volume (scalarVariable volume) dVscalar
        mydr = substitute dVscalar (scalarVariable dVscalar) $ var "dr" "\\Delta r" $ dVscalar ** (1.0/3)

findNamedSubexpression :: Type b => Expression b -> Maybe Exprn
findNamedSubexpression = searchExpressionDepthFirst Set.empty helper
  where helper e@(Var _ _ _ _ (Just _)) = Just $ mkExprn e
        helper _ = Nothing

-- In the following, we refuse to "find" a subexpression that has a
-- "k" in it, according to the hasK function.  This is to avoid
-- certain issues where we are unable to set k to zero correctly,
-- since we have created a variable that happens to be zero when k ==
-- 0, but which is divided by k (or by k^2).  This is a hokey
-- overkill, but I don't have a better idea.
findToDo :: Type b => Set.Set String -> [Exprn] -> Expression b -> Maybe Exprn
findToDo i everything = searchExpression i helper
  where helper (Sum _ ii) | Set.size ii < 2 = Nothing
        helper (Product _ ii) | Set.size ii < 2 = Nothing
        helper (Sum s _) | todo:_ <- filter simplifiable subes = Just todo
          where subes = map (mkExprn . pairs2sum) $ take numtotake $ subsq $
                        take numtotake $ filter (\(_,e) -> countVars [mkExprn e] > 0 &&
                                                           not (hasFFT e) &&
                                                           not (hasK e)) $ sum2pairs s
                simplifiable sube = countVarssube > 1 && countVarssube < 3 && ithelps sube
                  where countVarssube = countVars [sube]
                oldnum = countVars everything
                ithelps e = countAfterRemovalE e everything + 1 < oldnum
        helper (Product p _) | todo:_ <- filter simplifiable subes = Just todo
          where subes = map (mkExprn . pairs2product) $ take numtotake $ subsq $
                        take numtotake $ filter (\(e,_) -> countVars [mkExprn e] > 0 &&
                                                           not (hasFFT e) &&
                                                           not (hasK e)) $ product2pairs p
                simplifiable sube = countVarssube > 1 && countVarssube < 3 && ithelps sube
                  where countVarssube = countVars [sube]
                oldnum = countVars everything
                ithelps e = countAfterRemovalE e everything + 1 < oldnum
        helper _ = Nothing

findFFTtodo :: Type b => Expression b -> Maybe Exprn
findFFTtodo = searchExpression Set.empty helper
  where helper e@(Expression _)
          | EK (Expression (FFT (Var _ _ _ _ Nothing))) <- mkExprn e = Just $ mkExprn e
          | ER (Expression (IFFT (Var _ _ _ _ Nothing))) <- mkExprn e = Just $ mkExprn e
        helper _ = Nothing

findFFTinputtodo :: Type b => Set.Set String -> Expression b -> Maybe Exprn
findFFTinputtodo i = searchExpressionDepthFirst i helper
  where helper e@(Expression _)
          | EK (Expression (FFT e')) <- mkExprn e, not (hasFFT e') = Just $ ER e'
          | ER (Expression (IFFT e')) <- mkExprn e, not (hasFFT e') = Just $ EK e'
          | ES (Expression (Summate e')) <- mkExprn e, not (hasFFT e') = Just $ mkExprn e
        helper _ = Nothing

scalarVariable :: Expression Scalar -> Expression Scalar
scalarVariable (Var _ _ x t _) = Var CannotBeFreed x x t Nothing
scalarVariable _ = error "oopsisse"

optimize :: [Exprn] -> ([Statement], [Exprn])
optimize eee = case optimizeScalars [] 0 eee of
                 (a,_,b) -> (a,b)

-- Evaluate any purely scalar expressions that we can.  These are the
-- simplest computations we can do.  Currently this function also
-- enables all other optimizations, but I'd like to separate them out
-- again later so we can easily experiment with different sets of
-- optimizations.
optimizeScalars :: [Statement] -> Int -> [Exprn] -> ([Statement], Int, [Exprn])
optimizeScalars sts n everything = case mconcat $ map (mapExprn findNamedScalar) everything of
                                     Just (ES s@(Var _ _ _ _ (Just e))) ->
                                       case optimizeHelper Set.empty n [] everything [mkExprn e] of
                                         ([],_,_) -> optimizeScalars (sts++[Initialize (ES v), Assign (ES v) (ES s)]) n
                                                                     (map (mapExprn (mkExprn . substitute s v)) everything)
                                              where v = scalarVariable s
                                         (sts',n',everything') -> optimizeScalars (sts++sts') n' everything'
                                     Nothing -> optimizeHelper Set.empty (n :: Int) sts everything everything
                                     _ -> error "bad result in optimizeScalars"

-- Then we go looking for memory to save or ffts to evaluate...
optimizeHelper :: Set.Set String -> Int -> [Statement] -> [Exprn] -> [Exprn]
               -> ([Statement], Int, [Exprn])
optimizeHelper i n sts everything e = case handleSubstitution sts n everything e todos of
                                        Nothing -> (sts, n, everything)
                                        Just (vs, sts', n', everything', e') -> optimizeHelper vs n' sts' everything' e'
  where todos = if Set.size i == 0
                then [mconcat $ map (mapExprn (findToDo i everything)) e,
                      mconcat $ map (mapExprn findFFTtodo) e,
                      mconcat $ map (mapExprn (findFFTinputtodo i)) e]
                else [mconcat $ map (mapExprn (findToDo i everything)) e,
                      mconcat $ map (mapExprn findFFTtodo) e,
                      mconcat $ map (mapExprn (findFFTinputtodo i)) e,
                      mconcat $ map (mapExprn (findFFTinputtodo Set.empty)) e]

handleSubstitution :: [Statement] -> Int -> [Exprn] -> [Exprn] -> [Maybe Exprn]
                      -> Maybe (Set.Set String, [Statement], Int, [Exprn], [Exprn])
handleSubstitution _ _ _ _ [] = Nothing
handleSubstitution sts n e1 e2 (Nothing:todos) = handleSubstitution sts n e1 e2 todos
handleSubstitution sts n e1 e2 (Just (EK ke):_) = Just (varSet v,
                                                        sts++[Initialize (EK v), Assign (EK v) (EK ke)],
                                                        n+1,
                                                        map (mapExprn (mkExprn . substitute ke v)) e1,
                                                        map (mapExprn (mkExprn . substitute ke v)) e2)
  where v :: Expression KSpace
        v = case ke of
          Var _ xi x t _ -> Var IsTemp xi x t Nothing
          _ -> Var IsTemp ("ktemp" ++ show n++"[i]")
                          ("ktemp"++show n)
                          ("\\tilde{f}_{" ++ show n ++ "}") Nothing
handleSubstitution sts n e1 e2 (Just (ER re):_) = Just (varSet v,
                                                        sts++[Initialize (ER v), Assign (ER v) (ER re)],
                                                        n+1,
                                                        map (mapExprn (mkExprn . substitute re v)) e1,
                                                        map (mapExprn (mkExprn . substitute re v)) e2)
  where v :: Expression RealSpace
        v = case re of
          Var _ xi x t _ -> Var IsTemp xi x t Nothing
          _ -> Var IsTemp ("rtemp" ++ show n++"[i]")
                          ("rtemp"++show n)
                          ("f_{" ++ show n ++ "}") Nothing
handleSubstitution sts n e1 e2 (Just (ES se):_) = Just (varSet v,
                                                        sts++[Initialize (ES v), Assign (ES v) (ES se)],
                                                        n+1,
                                                        map (mapExprn (mkExprn . substitute se v)) e1,
                                                        map (mapExprn (mkExprn . substitute se v)) e2)
  where v :: Expression Scalar
        v = case se of
          Var _ xi x t _ -> Var IsTemp xi x t Nothing
          _ -> Var CannotBeFreed ("s" ++ show n)
                                 ("s"++show n)
                                 ("s_{" ++ show n ++ "}") Nothing

