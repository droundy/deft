import NewCode
import WhiteBear ( whitebear )
import System.Environment ( getArgs )

a1, a2, a3 :: Vector Scalar
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
       defineFunctional whitebear [ES $ s_var "R"] [ER $ r_var "x"] "HardSpheresNoTensor2"
     gentest "tests/new-generated-haskell/WhiteBear.h" $
       defineFunctional whitebear [ES $ s_var "R"] [ER $ r_var "x"] "WhiteBear"
     gentest "tests/new-generated-haskell/integrate_sqr.h" $
       defineFunctional (integrate $ r_var "nn"**2) [] [ER $ r_var "nn"] "integrate_sqr"
     let myvolume = "V" === a1 `dot` (a2 `cross` a3)
     gentest "tests/new-generated-haskell/volume_minus_one_sqr.h" $
       defineFunctional ((myvolume - 1)**2) [] [] "volume_minus_one_sqr"
