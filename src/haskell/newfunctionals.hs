import NewCode
import WhiteBear ( whitebear )
import System.Environment ( getArgs )
import FMT ( n )

a1, a2, a3 :: Vector Scalar
a1 = t_var "a1"
a2 = t_var "a2"
a3 = t_var "a3"

main :: IO ()
main =
  do todo <- getArgs
     let headerAndCpp fn f inp name =
           if ("src/"++fn++".cpp") `elem` todo
           then do writeFile ("src/"++fn++".h") $ createHeader f name
                   writeFile ("src/"++fn++".cpp") $ createCppFile f inp name (fn++".h")
           else return ()
     let gentest f x = if "tests" `elem` todo
                   then writeFile f x
                   else return ()
     let whitebearn = substitute n (r_var "n") whitebear
     headerAndCpp "new/WhiteBearFast" whitebearn [ER $ r_var "n"] "WhiteBear"
     gentest "tests/new-generated-haskell/WhiteBear.h" $
       defineFunctional whitebear [ER $ r_var "x"] "WhiteBear"
     gentest "tests/new-generated-haskell/integrate_sqr.h" $
       defineFunctional (integrate $ r_var "nn"**2) [ER $ r_var "nn"] "integrate_sqr"
     let myvolume = "V" === a1 `dot` (a2 `cross` a3)
     gentest "tests/new-generated-haskell/volume_minus_one_sqr.h" $
       defineFunctional ((myvolume - 1)**2) [] "volume_minus_one_sqr"
