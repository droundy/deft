module Latex where

import Expression ( Expression(..), Type, substitute, latex, k, k_var )
import Statement ( ToDo(..), findNamedSubexpression )

latexEasy :: (Type a) => Expression a -> String
latexEasy e0 = unlines $
              ["\\documentclass{article}",
               "\\usepackage{amsmath}",
               "\\usepackage{breqn}",
               "\\begin{document}"]
              ++ map latexme (niceToDos e0) ++
              ["\\end{document}"]
    where latexme (DoS e) = latexe e
          latexme (DoK e) = latexe e
          latexme (DoR e) = latexe e
          latexme DoNothing = error "oops in latexme"
          latexe (Var _ _ _ t (Just e')) = unlines ["\\begin{dmath}",
                                                    "  " ++ t ++ " = " ++ latex (substitute (k**2) (k_var "k" **2) e'),
                                                    "\\end{dmath}"]
          latexe _ = error "oops in latexe"

niceToDos :: (Type a) => Expression a -> [ToDo]
niceToDos e0@(Var _ _ _ _ _) =
  case findNamedSubexpression e0 of
    DoS v@(Var a b c t (Just _)) -> DoS v : niceToDos (substitute v (Var a b c t Nothing) e0)
    DoR v@(Var a b c t (Just _)) -> DoR v : niceToDos (substitute v (Var a b c t Nothing) e0)
    DoK v@(Var a b c t (Just _)) -> DoK v : niceToDos (substitute v (Var a b c t Nothing) e0)
    _ -> []
niceToDos _ = error "need named input in niceToDos"
