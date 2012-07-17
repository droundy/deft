import CodeGen
import SomeFunctionals
import Quantum
import System.Environment ( getArgs )

main :: IO ()
main =
  do todo <- getArgs
     let gen f x = if f `elem` todo
                   then writeFile f x
                   else return ()
     gen "src/HardSpheresNoTensor2Fast.cpp" $
       defineFunctional whitebear ["R"] "HardSpheresNoTensor2"
     gen "src/Correlation_S_Fast.cpp" $
       generateHeader correlation_S_WB ["R"] "Correlation_S2"
     gen "src/Correlation_A_Fast.cpp" $
       generateHeader correlation_A_WB ["R"] "Correlation_A2"
     gen "src/YuWuCorrelationFast.cpp" $
       generateHeader yuwu_correlation ["R"] "YuWuCorrelationFast"
     gen "src/SaftFluid2Fast.cpp" $
       defineFunctional saft_fluid ["R", "epsilon_association", "kappa_association",
                                    "epsilon_dispersion", "lambda_dispersion", "length_scaling",
                                    "mu"] "SaftFluid2"
     gen "src/Association2Fast.cpp" $
       defineFunctional saft_association
       ["R", "epsilon_association", "kappa_association",
        "epsilon_dispersion", "lambda_dispersion", "length_scaling"] "Association2"
     gen "src/Dispersion2Fast.cpp" $
       defineFunctional saft_dispersion
       ["R", "epsilon_dispersion", "lambda_dispersion", "length_scaling"] "Dispersion2"
     let psi = r_var "x"
	 chrisfunc = psi**3/(1+3*psi)*(exp (-psi))*((psi-psi**2)**0.5)
     gen "src/haskell/ChrisFunctional.tex" $ latex $ chrisfunc
     gen "src/haskell/ChrisFunctional.cpp" $ generateHeader chrisfunc ["R"] "ChrisFunctional"
     gen "papers/contact/formulas.tex" $ unlines $
       map definelatex [("phione", phi1),
                        ("phitwo", phi2),
                        ("phithree", phi2),
                        ("yuwucorrelation", yuwu_correlation),
                        ("dwbdnthree", dwbdn3),
                        ("dwbdntwo", dwbdn2),
                        ("dwbdnone", dwbdn1),
                        ("dwbdnonev", dwbdn1v_over_n2v),
                        ("dwbdntwov", dwbdn2v_over_n2v)] ++
       map definescalar [("whitebear", whitebear)]
     gen "src/HydrogenFast.cpp" $
       generateHeader (oneElectron hydrogenPotential) [] "Hydrogen"
definelatex :: (String, Expression RealSpace) -> String
definelatex (v,e) = "\\newcommand\\" ++ v ++ "{" ++ latex e ++ "}"
definescalar :: (String, Expression Scalar) -> String
definescalar (v,e) = "\\newcommand\\" ++ v ++ "{" ++ latex e ++ "}"
