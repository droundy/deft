import HughesSaft ( saft_fluid, saft_entropy, saft_association, saft_dispersion )
import WhiteBear ( whitebear, tensorwhitebear, whitebear_m2, gSigmaA, gSigmaS, gSigmaA_m2 )
import Latex
import Expression ( (===), grad, joinFFTs )
import System.Environment ( getArgs )
import System.Process ( rawSystem )
import System.Exit ( ExitCode(ExitSuccess) )
import System.FilePath ( dropFileName )

main :: IO ()
main =
  do todo <- getArgs
     let pdf f x = if f `elem` todo
                   then do let texname = reverse (drop 3 $ reverse f) ++ "tex"
                           writeFile texname x
                           ec <- rawSystem "pdflatex" ["-output-directory", dropFileName f, texname]
                           if ec /= ExitSuccess
                              then fail $ f ++ " failed with exit code " ++ show ec
                              else return ()
                   else return ()
     pdf "doc/WhiteBear.pdf" $ latexEasy $ whitebear
     pdf "doc/TensorWhiteBear.pdf" $ latexEasy $ tensorwhitebear
     pdf "doc/GradTensorWhiteBear.pdf" $ latexEasy $ "grad" === (grad "x" tensorwhitebear)
     pdf "doc/WhiteBearMarkII.pdf" $ latexEasy $ whitebear_m2
     pdf "doc/gSigmaS.pdf" $ latexEasy $ gSigmaS
     pdf "doc/gSigmaA.pdf" $ latexEasy $ gSigmaA
     pdf "doc/gSigmaAm2.pdf" $ latexEasy $ gSigmaA_m2
     pdf "doc/Association.pdf" $ latexEasy $ saft_association
     pdf "doc/GradAssociation.pdf" $ latexEasy $ "grad" === (grad "x" saft_association)
     pdf "doc/Dispersion.pdf" $ latexEasy $ saft_dispersion
     pdf "doc/SimpDispersion.pdf" $ latexOptimizedExpression $ saft_dispersion
     pdf "doc/SaftFluid.pdf" $ latexEasy $ saft_fluid
     pdf "doc/EntropySaftFluid.pdf" $ latexEasy $ saft_entropy
     pdf "doc/GradDispersion.pdf" $ latexEasy $ "grad" === (grad "x" saft_dispersion)
     pdf "doc/JoinedGradDispersion.pdf" $ latexEasy $ "grad" === (joinFFTs $ grad "x" saft_dispersion)
     pdf "doc/SimpGradDispersion.pdf" $ latexOptimizedExpression $ "grad" === (joinFFTs $ grad "x" saft_dispersion)
     pdf "doc/GradAssociation.pdf" $ latexEasy $ "grad" === (grad "x" saft_association)
     pdf "doc/JoinedGradAssociation.pdf" $ latexEasy $ "grad" === (joinFFTs $ grad "x" saft_association)
