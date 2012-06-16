module Latex ( latexEasy, latexSimp ) where

import Expression ( Expression(..), Type, substitute, latex, k, k_var, cleanvars )
import Statement ( ToDo(..), Statement(..), findNamedSubexpression, simp2, freeVectors, reuseVar )

latexEasy :: (Type a) => Expression a -> String
latexEasy e0 = unlines $
              ["\\documentclass{article}",
               "\\usepackage{amsmath}",
               "\\usepackage{breqn}",
               "\\begin{document}"]
              ++ map latexme (niceToDos e0) ++
              ["\\end{document}"]

latexSimp :: (Type a) => Expression a -> String
latexSimp e = unlines $
              ["\\documentclass{article}",
               "\\usepackage{amsmath}",
               "\\usepackage{breqn}",
               "\\begin{document}"]
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

latexme :: ToDo -> String
latexme (DoS e) = latexe e
latexme (DoK e) = latexe e
latexme (DoR e) = latexe e
latexme DoNothing = error "oops in latexme"

latexe :: Type a => Expression a -> String
latexe (Var _ _ _ t (Just e')) = eqn (t ++ " = " ++ latex (simpk e'))
latexe (Var _ _ _ _ Nothing) = ""
latexe e = error ("oops in latexe: " ++ show e)

simpk :: Type a => Expression a -> Expression a
simpk = substitute (k**2) (k_var "k" **2)

niceToDos :: (Type a) => Expression a -> [ToDo]
niceToDos e0@(Var _ _ _ _ _) =
  case findNamedSubexpression e0 of
    DoS v@(Var a b c t (Just _)) -> DoS v : niceToDos (substitute v (Var a b c t Nothing) e0)
    DoR v@(Var a b c t (Just _)) -> DoR v : niceToDos (substitute v (Var a b c t Nothing) e0)
    DoK v@(Var a b c t (Just _)) -> DoK v : niceToDos (substitute v (Var a b c t Nothing) e0)
    _ -> []
niceToDos _ = error "need named input in niceToDos"


latexS :: Statement -> String
latexS (AssignR x y) = eqn $ latex x ++ " = " ++ latex (simpk $ cleanvars y)
latexS (AssignK x y) = eqn $ latex x ++ " = " ++ latex (simpk $ cleanvars y)
latexS (AssignS x y) = eqn $ latex x ++ " = " ++ latex (simpk $ cleanvars y)
latexS (InitializeR e) = eqn $ "\\text{initialize real space } " ++ latex e
latexS (InitializeK e) = eqn $ "\\text{initialize reciprocal space } " ++ latex e
latexS (InitializeS _) = ""
latexS (FreeR e) = eqn $ "\\text{free real space } " ++ latex e
latexS (FreeK e) = eqn $ "\\text{free reciprocal space } " ++ latex e
