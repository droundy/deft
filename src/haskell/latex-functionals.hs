import SomeFunctionals
import Latex
import System.Environment ( getArgs )
import System.Process ( rawSystem )
import System.FilePath ( dropFileName )

main :: IO ()
main =
  do todo <- getArgs
     let pdf f x = if f `elem` todo
                   then do let texname = reverse (drop 3 $ reverse f) ++ "tex"
                           writeFile texname x
                           rawSystem "pdflatex" ["-output-directory", dropFileName f, texname]
                           return ()
                   else return ()
     pdf "doc/WhiteBear.pdf" $ latexEasy whitebear
     pdf "doc/Association.pdf" $ latexEasy saft_association
     pdf "doc/Dispersion.pdf" $ latexEasy saft_dispersion
