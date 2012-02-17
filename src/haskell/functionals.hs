import CodeGen
import SomeFunctionals
import System.Environment ( getArgs )

main :: IO ()
main = 
  do todo <- getArgs
     let gen f x = if f `elem` todo
                   then writeFile f x
                   else return ()
     gen "src/HardSpheresNoTensor2Fast.cpp" $ 
       generateHeader (fmt whitebear) ["R"] "HardSpheresNoTensor2"
     gen "src/ContactAtSphereFast.cpp" $
       generateHeader (fmt wb_contact_at_sphere) ["R"] "ContactAtSphere"
     gen "src/YuWuContactFast.cpp" $
       generateHeader (fmt yuwu_contact) ["R"] "YuWuContact"
     gen "src/SaftFluid2Fast.cpp" $
       generateHeader saft_fluid ["R"] "SaftFluid2"
     gen "papers/contact/formulas.tex" $ unlines $ map definelatex
       [("phione", phi1),
        ("phitwo", phi2),
        ("phithree", phi2),
        ("whitebear", whitebear),
        ("yuwucontact", yuwu_contact),
        ("dwbdnthree", dwbdn3),
        ("dwbdntwo", dwbdn2),
        ("dwbdnone", dwbdn1),
        ("dwbdnonev", dwbdn1v_over_n2v),
        ("dwbdntwov", dwbdn2v_over_n2v)]

definelatex :: (String, Expression RealSpace) -> String
definelatex (v,e) = "\\newcommand\\" ++ v ++ "{" ++ latex e ++ "}"
