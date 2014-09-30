import CodeGen
import HughesSaft ( saft_fluid, saft_entropy, yuwu_correlation, hughes_X, hughes_HB )
import WaterSaft ( water_saft, water_entropy, water_X, mu, water_saft_by_hand )
import IdealGas ( idealgas )
import FMT ( n, n2, n2mxx, n2x )
import SFMT ( sfmt )
import Rosenfeld ( fmt )
import WhiteBear ( whitebear, tensorwhitebear, whitebear_m2, correlation_gross, gSigmaA, gSigmaS,
                   gSigmaA_m2, gSigmaS_m2, gSigmaA_by_hand, gSigmaA_automagic )
import System.Environment ( getArgs )

main :: IO ()
main =
  do todo <- getArgs
     let gen f x = if f `elem` todo
                   then writeFile f x
                   else return ()
     let nmu = "nmu" === integrate (n*mu)
     gen "src/SoftFluidFast.cpp" $
       defineFunctional (idealgas + sfmt + nmu) ["sigma", "epsilon", "mu"] "SoftFluid"
     gen "src/HardFluidFast.cpp" $
       defineFunctional (idealgas + whitebear + nmu) ["R", "mu"] "HardFluid"
     gen "src/HardRosenfeldFluidFast.cpp" $
       defineFunctional (idealgas + fmt + nmu) ["R", "mu"] "HardRosenfeldFluid"
     
     gen "src/HardSpheresNoTensor2Fast.cpp" $
       defineFunctional whitebear ["R"] "HardSpheresNoTensor2"
     gen "src/WhiteBearMarkIIFast.cpp" $
       defineFunctional whitebear_m2 ["R"] "WhiteBearMarkII"
     gen "src/TensorWhiteBearFast.cpp" $
       defineFunctional tensorwhitebear ["R"] "TensorWhiteBear"
     gen "src/n2DensityFast.cpp" $
       defineTransformation n2 ["R"] "n2Density"
     gen "src/TensorDensityXXFast.cpp" $
       defineTransformation n2mxx ["R"] "TensorDensityXX"
     gen "src/VectorDensityXFast.cpp" $
       defineTransformation n2x ["R"] "VectorDensityX"
     gen "src/gSigmaSFast.cpp" $
       defineTransformation gSigmaS ["R"] "gSigmaS"
     gen "src/gSigmaAFast.cpp" $
       defineTransformation gSigmaA ["R"] "gSigmaA"

     gen "src/gSigmaA_by_handFast.cpp" $
       generateHeader gSigmaA_by_hand ["R"] "gSigmaA_by_hand"
     gen "src/gSigmaA_automagicFast.cpp" $
       generateHeader gSigmaA_automagic ["R"] "gSigmaA_automagic"

     gen "src/gSigmaSm2Fast.cpp" $
       defineTransformation gSigmaS_m2 ["R"] "gSigmaSm2"

     gen "src/CorrelationGrossCorrectFast.cpp" $
       defineTransformation correlation_gross ["R"] "CorrelationGrossCorrect"

     gen "src/gSigmaAm2Fast.cpp" $
       defineTransformation gSigmaA_m2 ["R"] "gSigmaAm2"
     gen "src/YuWuCorrelationFast.cpp" $
       defineTransformation yuwu_correlation ["R"] "YuWuCorrelationFast"
     gen "src/SaftFluid2Fast.cpp" $
       defineFunctional saft_fluid ["R", "epsilon_association", "kappa_association",
                                    "epsilon_dispersion", "lambda_dispersion", "length_scaling",
                                    "mu"] "SaftFluid2"
     gen "src/WaterSaftFast.cpp" $
       defineFunctional water_saft ["R", "epsilon_association", "kappa_association",
                                    "epsilon_dispersion", "lambda_dispersion", "length_scaling",
                                    "mu"] "WaterSaft"
     gen "src/WaterSaft_by_handFast.cpp" $
       defineFunctional water_saft_by_hand ["R", "epsilon_association", "kappa_association",
                                            "epsilon_dispersion", "lambda_dispersion", "length_scaling",
                                            "mu"] "WaterSaft_by_hand"
     gen "src/WaterXFast.cpp" $
       defineTransformation water_X ["R", "epsilon_association", "kappa_association",
                                     "epsilon_dispersion", "lambda_dispersion", "length_scaling",
                                     "mu"] "WaterX"
     gen "src/HughesXFast.cpp" $
       defineTransformation hughes_X ["R", "epsilon_association", "kappa_association",
                                      "epsilon_dispersion", "lambda_dispersion", "length_scaling",
                                      "mu"] "HughesX"
     gen "src/HughesHBFast.cpp" $
       defineTransformation hughes_HB ["R", "epsilon_association", "kappa_association",
                                      "epsilon_dispersion", "lambda_dispersion", "length_scaling",
                                      "mu"] "HughesHB"
     gen "src/WaterEntropyFast.cpp" $
       defineFunctional water_entropy ["R", "epsilon_association", "kappa_association",
                                    "epsilon_dispersion", "lambda_dispersion", "length_scaling",
                                    "mu"] "WaterEntropy"
     gen "src/EntropySaftFluid2Fast.cpp" $
       defineFunctionalNoGradient saft_entropy
       ["R", "epsilon_association", "kappa_association",
        "epsilon_dispersion", "lambda_dispersion", "length_scaling"] "EntropySaftFluid2"
