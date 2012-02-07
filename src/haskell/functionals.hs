import CodeGen
import SomeFunctionals

main :: IO ()
main = 
  do writeFile "src/HardSpheresNoTensor2Fast.cpp" $ 
       generateHeader whitebear (Just (r_var "R")) "HardSpheresNoTensor2"
     --writeFile "src/HardSphereGas2Fast.cpp" $
     --  generateHeader (whitebear + idealgas + mu*n) ["R", "mu"] "HardSphereGasFast2"
