import SomeFunctionals
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
     pdf "doc/WhiteBear.pdf" $ latexEasy $ "FHS" === whitebear
     pdf "doc/Association.pdf" $ latexEasy $ "Fassoc" === saft_association
     pdf "doc/Dispersion.pdf" $ latexEasy $ "Fdisp" === saft_dispersion
     pdf "doc/SaftFluid.pdf" $ latexEasy $ "Fw" === saft_fluid
     let foo = "Fdisp" === integrate saft_dispersion
     pdf "doc/GradDispersion.pdf" $ latexEasy $ "grad" === (grad "x" foo)
     pdf "doc/JoinedGradDispersion.pdf" $ latexEasy $ "grad" === (joinFFTs $ grad "x" foo)
     pdf "doc/SimpGradDispersion.pdf" $ latexSimp $ "grad" === (joinFFTs $ grad "x" $ integrate saft_dispersion)
     let assoc = "Fassoc" === integrate saft_association
     pdf "doc/GradAssociation.pdf" $ latexEasy $ "grad" === (grad "x" assoc)
     pdf "doc/JoinedGradAssociation.pdf" $ latexEasy $ "grad" === (joinFFTs $ grad "x" assoc)
     pdf "doc/JoinedGradHydrogen.pdf" $ latexEasy $ "grad" === (joinFFTs $ grad "x" (integrate (oneElectron hydrogenPotential)))
