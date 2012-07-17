import SomeFunctionals
import WhiteBear ( whitebear, correlation_A_WB, correlation_S_WB )
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
     pdf "doc/CorrelationS.pdf" $ latexEasy $ correlation_S_WB
     pdf "doc/CorrelationA.pdf" $ latexEasy $ correlation_A_WB
     pdf "doc/Association.pdf" $ latexEasy $ saft_association
     pdf "doc/GradAssociation.pdf" $ latexEasy $ "grad" === (grad "x" saft_association)
     pdf "doc/Dispersion.pdf" $ latexEasy $ saft_dispersion
     pdf "doc/SimpDispersion.pdf" $ latexSimp $ saft_dispersion
     pdf "doc/SaftFluid.pdf" $ latexEasy $ saft_fluid
     pdf "doc/GradDispersion.pdf" $ latexEasy $ "grad" === (grad "x" saft_dispersion)
     pdf "doc/JoinedGradDispersion.pdf" $ latexEasy $ "grad" === (joinFFTs $ grad "x" saft_dispersion)
     pdf "doc/SimpGradDispersion.pdf" $ latexSimp $ "grad" === (joinFFTs $ grad "x" saft_dispersion)
     pdf "doc/GradAssociation.pdf" $ latexEasy $ "grad" === (grad "x" saft_association)
     pdf "doc/JoinedGradAssociation.pdf" $ latexEasy $ "grad" === (joinFFTs $ grad "x" saft_association)
     pdf "doc/JoinedGradHydrogen.pdf" $ latexEasy $ "grad" === (joinFFTs $ grad "x" (integrate (oneElectron hydrogenPotential)))
