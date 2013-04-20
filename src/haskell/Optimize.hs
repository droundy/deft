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

findToDo :: Type b => Set.Set String -> [Exprn] -> Expression b -> Maybe Exprn
findToDo i everything = searchExpression i helper
  where helper (Sum _ ii) | Set.size ii < 2 = Nothing
        helper (Product _ ii) | Set.size ii < 2 = Nothing
        helper (Sum s _) | todo:_ <- filter simplifiable subes = Just todo
          where subes = map (mkExprn . pairs2sum) $ take numtotake $ subsq $
                        take numtotake $ filter (\(_,e) -> countVars [mkExprn e] > 0 && not (hasFFT e)) $ sum2pairs s
                simplifiable sube = countVarssube > 1 && countVarssube < 3 && ithelps sube
                  where countVarssube = countVars [sube]
                oldnum = countVars everything
                ithelps e = countAfterRemovalE e everything + 1 < oldnum
        helper (Product p _) | todo:_ <- filter simplifiable subes = Just todo
          where subes = map (mkExprn . pairs2product) $ take numtotake $ subsq $
                        take numtotake $ filter (\(e,_) -> countVars [mkExprn e] > 0 &&
                                                           not (hasFFT e)) $ product2pairs p
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
optimize eee = case scalarhelper [] 0 eee of
            (a,_,b) -> (a,b)
    where -- First, we want to evaluate any purely scalar expressions
          -- that we can, just to simplify things!
          scalarhelper :: [Statement] -> Int -> [Exprn] -> ([Statement], Int, [Exprn])
          scalarhelper sts n everything =
            case mconcat $ map (mapExprn findNamedScalar) everything of
              Just (ES s@(Var _ _ _ _ (Just e))) ->
                case optimizeHelper Set.empty n [] everything [mkExprn e] of
                  ([],_,_) -> scalarhelper (sts++[Initialize (ES v), Assign (ES v) (ES s)]) n
                                           (map (mapExprn (mkExprn . substitute s v)) everything)
                              where v = scalarVariable s
                  (sts',n',everything') -> scalarhelper (sts++sts') n' everything'
              Nothing -> optimizeHelper Set.empty (n :: Int) sts everything everything
              _ -> error "bad result in scalarhelper"

-- Then we go looking for memory to save or ffts to evaluate...
optimizeHelper :: Set.Set String -> Int -> [Statement] -> [Exprn] -> [Exprn]
               -> ([Statement], Int, [Exprn])
optimizeHelper i n sts everything e =
  if Set.size i == 0
  then handletodos [mconcat $ map (mapExprn (findToDo i everything)) e,
                    mconcat $ map (mapExprn findFFTtodo) e,
                    mconcat $ map (mapExprn (findFFTinputtodo i)) e]
  else handletodos [mconcat $ map (mapExprn (findToDo i everything)) e,
                    mconcat $ map (mapExprn findFFTtodo) e,
                    mconcat $ map (mapExprn (findFFTinputtodo i)) e,
                    mconcat $ map (mapExprn (findFFTinputtodo Set.empty)) e]
  where handletodos [] = (sts, n, everything)
        handletodos (todo:ts) =
          case todo of
            Just (EK ke) -> optimizeHelper (varSet v) (n+1)
                            (sts++[Initialize (EK v), Assign (EK v) (EK ke)])
                            (map (mapExprn (mkExprn . substitute ke v)) everything)
                            (map (mapExprn (mkExprn . substitute ke v)) e)
              where v :: Expression KSpace
                    v = case ke of
                      Var _ xi x t _ -> Var IsTemp xi x t Nothing
                      _ -> Var IsTemp ("ktemp" ++ show n++"[i]")
                                      ("ktemp"++show n)
                                      ("\\tilde{f}_{" ++ show n ++ "}") Nothing
            Just (ER re) -> optimizeHelper (varSet v) (n+1)
                            (sts++[Initialize (ER v), Assign (ER v) (ER re)])
                            (map (mapExprn (mkExprn . substitute re v)) everything)
                            (map (mapExprn (mkExprn . substitute re v)) e)
              where v :: Expression RealSpace
                    v = case re of
                      Var _ xi x t _ -> Var IsTemp xi x t Nothing
                      _ -> Var IsTemp ("rtemp" ++ show n++"[i]")
                                      ("rtemp"++show n)
                                      ("f_{" ++ show n ++ "}") Nothing
            Just (ES se) -> optimizeHelper (varSet v) (n+1)
                            (sts++[Initialize (ES v), Assign (ES v) (ES se)])
                            (map (mapExprn (mkExprn . substitute se v)) everything)
                            (map (mapExprn (mkExprn . substitute se v)) e)
              where v :: Expression Scalar
                    v = case se of
                      Var _ xi x t _ -> Var IsTemp xi x t Nothing
                      _ -> Var CannotBeFreed ("s" ++ show n)
                                             ("s"++show n)
                                             ("s_{" ++ show n ++ "}") Nothing
            Nothing -> handletodos ts

