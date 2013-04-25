import HughesSaft ( saft_fluid, saft_entropy, saft_association, saft_dispersion )
import WhiteBear ( whitebear, whitebear_m2, correlation_A_WB, correlation_S_WB, correlation_A_WB_m2 )
import Quantum
import Latex
import Expression ( (===), grad, joinFFTs, integrate )
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
     pdf "doc/WhiteBearMarkII.pdf" $ latexEasy $ whitebear_m2
     pdf "doc/CorrelationS.pdf" $ latexEasy $ correlation_S_WB
     pdf "doc/CorrelationA.pdf" $ latexEasy $ correlation_A_WB
     pdf "doc/CorrelationAm2.pdf" $ latexEasy $ correlation_A_WB_m2
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
     pdf "doc/JoinedGradHydrogen.pdf" $ latexEasy $ "grad" === (joinFFTs $ grad "x" (integrate (oneElectron hydrogenPotential)))
