import CodeGen
import SomeFunctionals

main :: IO ()
main = 
  do writeFile "src/HardSpheresNoTensor2Fast.cpp" $ 
       generateHeader (fmt whitebear) (Just (r_var "R")) "HardSpheresNoTensor2"
     writeFile "src/ContactAtSphereFast.cpp" $
       generateHeader (fmt wb_contact_at_sphere) (Just (r_var "R")) "ContactAtSphere"
     writeFile "papers/contact/formulas.tex" $ unlines $ map definelatex
       [("phione", phi1),
        ("phitwo", phi2),
        ("phithree", phi2),
        ("whitebear", whitebear),
        ("dwbdnthree", dwbdn3),
        ("dwbdntwo", dwbdn2),
        ("dwbdnone", dwbdn1),
        ("dwbdnonev", dwbdn1v_over_n2v),
        ("dwbdntwov", dwbdn2v_over_n2v)]
     --writeFile "src/HardSphereGas2Fast.cpp" $
     --  generateHeader (whitebear + idealgas + mu*n) ["R", "mu"] "HardSphereGasFast2"

definelatex :: (String, Expression RealSpace) -> String
definelatex (v,e) = "\\newcommand\\" ++ v ++ "{" ++ latex e ++ "}"
