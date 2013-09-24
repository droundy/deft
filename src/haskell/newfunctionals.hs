import NewCode
import WhiteBear ( whitebear )
import WaterSaft ( water_saft, water_saft_by_hand )
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
         homogeneous_whitebear = makeHomogeneous whitebearn
     headerAndCpp "new/WhiteBearFast" whitebearn [ER $ r_var "n"] "WhiteBear"
     headerAndCpp "new/HomogeneousWhiteBearFast" homogeneous_whitebear [ES $ s_var "n"] "HomogeneousWhiteBear"

     let water_saftn = substitute n (r_var "n") water_saft
         water_saft_by_handn = substitute n (r_var "n") water_saft_by_hand
     headerAndCpp "new/HomogeneousWaterSaftFast" (makeHomogeneous water_saftn) [ES $ s_var "n"] "HomogeneousWaterSaft"
     headerAndCpp "new/HomogeneousWaterSaftByHandFast" (makeHomogeneous water_saft_by_handn)
       [ES $ s_var "n"] "HomogeneousWaterSaftByHand"

     -- The following two are really almost exclusively useful for
     -- debugging at the moment, as they do not compute any gradients.
     -- That is what the empty list [] below is.  On the plus side,
     -- they don't take so long to generate the source code this way.
     headerAndCpp "new/WaterSaftFast" water_saftn [] "WaterSaft"
     headerAndCpp "new/WaterSaftByHandFast" water_saft_by_handn [] "WaterSaftByHand"

     gentest "tests/new-generated-haskell/WhiteBear.h" $
       defineFunctional whitebear [ER $ r_var "x"] "WhiteBear"
     gentest "tests/new-generated-haskell/integrate_sqr.h" $
       defineFunctional (integrate $ r_var "nn"**2) [ER $ r_var "nn"] "integrate_sqr"
     let myvolume = "V" === a1 `dot` (a2 `cross` a3)
     gentest "tests/new-generated-haskell/volume_minus_one_sqr.h" $
       defineFunctional ((myvolume - 1)**2) [] "volume_minus_one_sqr"
