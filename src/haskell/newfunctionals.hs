import NewCode
import SomeFunctionals ( whitebear )
import System.Environment ( getArgs )

a1, a2, a3 :: Expression ThreeVector
a1 = t_var "a1"
a2 = t_var "a2"
a3 = t_var "a3"

main :: IO ()
main =
  do todo <- getArgs
     let gen f x = if f `elem` todo
                   then writeFile f x
                   else return ()
     let gentest f x = if "tests" `elem` todo
                   then writeFile f x
                   else return ()
     gen "src/new/WhiteBearFast.cpp" $
       defineFunctional whitebear ["R"] "HardSpheresNoTensor2"
     gentest "tests/new-generated-haskell/volume_minus_one_sqr.h" $
       defineFunctional ((volume - 1)**2) [] "volume_minus_one_sqr"
         where volume = "V" === a1 `dot` (a2 `cross` a3)
