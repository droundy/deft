{-# LANGUAGE PatternGuards #-}
module Latex ( latexEasy, latexSimp ) where

import Expression ( Expression(..), Exprn(..), Type, substitute, latex, k, k_var, cleanvars )
import Statement ( Statement(..), findNamedSubexpression, simp2, freeVectors, reuseVar )

latexEasy :: (Type a) => Expression a -> String
latexEasy e0 = unlines $
              ["\\documentclass{article}",
               "\\usepackage{amsmath}",
               "\\usepackage{color}",
               "\\usepackage{breqn}",
               "\\begin{document}"]
              ++ map latexme (niceExprns e0) ++
              ["\\end{document}"]

latexSimp :: (Type a) => Expression a -> String
latexSimp e = unlines $
              ["\\documentclass{article}",
               "\\usepackage{environ}", -- consider mdframed when it is in debian
               "\\usepackage{amsmath}",
               "\\usepackage{color}",
               "\\usepackage{breqn}",
               "\\begin{document}",
               "\\newbox{\\savetextbox}",
               "\\NewEnviron{bracetext}[1][\\textwidth]{%",
               "  \\begin{lrbox}{\\savetextbox}%",
               "    \\begin{minipage}{#1} \\BODY \\end{minipage}",
               "  \\end{lrbox}%",
               "  \\smallskip%",
               "  \\noindent\\makebox[0pt][r]{$\\left\\{\\rule{0pt}{\\ht\\savetextbox}\\right.$}%",
               "  \\usebox{\\savetextbox}\\par",
               "  \\smallskip%",
               "}",
               ""]
              ++ map latexS (sts) ++
              [latexe e',
               "\\end{document}"]
    where (sts0,e') = simp2 e
          sts = reuseVar $ freeVectors sts0

eqn :: String -> String
-- The following is a safety valve to avoid applying dmath to
-- equations that are so complicated that it simply can't handle them.
-- The number of characters is determined by pure experimentation, and
-- may be quite inaccurate.
eqn e | length e > 10000 =  unlines ["\\begin{equation}",
                                     "  " ++ e,
                                     "\\end{equation}"]
eqn e = unlines ["\\begin{dmath}",
                 "  " ++ e,
                 "\\end{dmath}"]

latexme :: Exprn -> String
latexme (ES e) = latexe e
latexme (EK e) = latexe e
latexme (ER e) = latexe e

latexe :: Type a => Expression a -> String
latexe (Var _ _ _ t (Just e')) = eqn (t ++ " = " ++ latex (simpk e'))
latexe (Var _ _ _ _ Nothing) = ""
latexe e = error ("oops in latexe: " ++ show e)

simpk :: Type a => Expression a -> Expression a
simpk = substitute (k**2) (k_var "k" **2)

niceExprns :: (Type a) => Expression a -> [Exprn]
niceExprns e0@(Var _ _ _ _ _) =
  case findNamedSubexpression e0 of
    Just (ES v@(Var a b c t (Just _))) -> ES v : niceExprns (substitute v (Var a b c t Nothing) e0)
    Just (ER v@(Var a b c t (Just _))) -> ER v : niceExprns (substitute v (Var a b c t Nothing) e0)
    Just (EK v@(Var a b c t (Just _))) -> EK v : niceExprns (substitute v (Var a b c t Nothing) e0)
    _ -> []
niceExprns _ = error "need named input in niceExprns"


latexS :: Statement -> String
latexS (Assign (ER x@(Var a b c t _)) (ER y))
  | tds@(_:_:_) <- niceExprns (Var a b c t (Just y)) = unlines $ ["\\begin{bracetext}"] ++
                                                      map latexme tds ++
                                                      ["\\end{bracetext}"]
  | otherwise  = eqn $ latex x ++ " = " ++ latex (simpk $ cleanvars y)
latexS (Assign (EK x@(Var a b c t _)) (EK y))
  | tds@(_:_:_) <- niceExprns (Var a b c t (Just y)) = unlines $ ["\\begin{bracetext}"] ++
                                                      map latexme tds ++
                                                      ["\\end{bracetext}}"]
  | otherwise  = eqn $ latex x ++ " = " ++ latex (simpk $ cleanvars y)
latexS (Assign (ES x@(Var a b c t _)) (ES y))
  | tds@(_:_:_) <- niceExprns (Var a b c t (Just y)) = unlines $ ["\\begin{bracetext}"] ++
                                                      map latexme tds ++
                                                      ["\\end{bracetext}"]
  | otherwise  = eqn $ latex x ++ " = " ++ latex (simpk $ cleanvars y)
latexS (Initialize (ER e)) = eqn $ "\\text{initialize real space } " ++ latex e
latexS (Initialize (EK e)) = eqn $ "\\text{initialize reciprocal space } " ++ latex e
latexS (Initialize (ES _)) = ""
latexS (Free e) = eqn $ "\\text{free } " ++ latex e
latexS st = error ("bad statement in latexS: "++ show st)
