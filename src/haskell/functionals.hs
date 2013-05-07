import CodeGen
import HughesSaft ( saft_fluid, saft_entropy, yuwu_correlation, mu )
import WaterSaft ( water_saft )
import IdealGas ( idealgas )
import FMT ( n )
import SFMT ( sfmt )
import WhiteBear ( whitebear, correlation_gross, gSigmaA, gSigmaS, 
                   gSigmaA_m2, gSigmaS_m2 )
import Quantum
import System.Environment ( getArgs )

main :: IO ()
main =
  do todo <- getArgs
     let gen f x = if f `elem` todo
                   then writeFile f x
                   else return ()
     let nmu = "nmu" === integrate (n*mu)
     gen "src/SoftFluidFast.cpp" $
       defineFunctional (idealgas + sfmt + nmu) ["R", "V0", "mu"] "SoftFluid"
     gen "src/HardFluidFast.cpp" $
       defineFunctional (idealgas + whitebear + nmu) ["R", "mu"] "HardFluid"
     
     gen "src/HardSpheresNoTensor2Fast.cpp" $
       defineFunctional whitebear ["R"] "HardSpheresNoTensor2"
     gen "src/gSigmaSFast.cpp" $
       generateHeader gSigmaS ["R"] "gSigmaS"
     gen "src/gSigmaAFast.cpp" $
       generateHeader gSigmaA ["R"] "gSigmaA"

     gen "src/gSigmaSm2Fast.cpp" $
       generateHeader gSigmaS_m2 ["R"] "gSigmaSm2"

     gen "src/CorrelationGrossCorrectFast.cpp" $
       generateHeader correlation_gross ["R"] "CorrelationGrossCorrect"

     gen "src/gSigmaAm2Fast.cpp" $
       generateHeader gSigmaA_m2 ["R"] "gSigmaAm2"
     gen "src/YuWuCorrelationFast.cpp" $
       generateHeader yuwu_correlation ["R"] "YuWuCorrelationFast"
     gen "src/SaftFluid2Fast.cpp" $
       defineFunctional saft_fluid ["R", "epsilon_association", "kappa_association",
                                    "epsilon_dispersion", "lambda_dispersion", "length_scaling",
                                    "mu"] "SaftFluid2"
     gen "src/WaterSaftFast.cpp" $
       defineFunctional water_saft ["R", "epsilon_association", "kappa_association",
                                    "epsilon_dispersion", "lambda_dispersion", "length_scaling",
                                    "mu"] "WaterSaft"
     gen "src/EntropySaftFluid2Fast.cpp" $
       defineFunctionalNoGradient saft_entropy
       ["R", "epsilon_association", "kappa_association",
        "epsilon_dispersion", "lambda_dispersion", "length_scaling"] "EntropySaftFluid2"
     let psi = r_var "x"
	 chrisfunc = psi**3/(1+3*psi)*(exp (-psi))*((psi-psi**2)**0.5)
     gen "src/haskell/ChrisFunctional.tex" $ latex $ chrisfunc
     gen "src/haskell/ChrisFunctional.cpp" $ generateHeader chrisfunc ["R"] "ChrisFunctional"
     gen "src/HydrogenFast.cpp" $
       generateHeader (oneElectron hydrogenPotential) [] "Hydrogen"
