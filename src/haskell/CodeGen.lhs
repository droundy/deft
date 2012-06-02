\begin{code}
{-# LANGUAGE GADTs, PatternGuards #-}

module CodeGen (module Statement,
                module Expression,
                generateHeader)
    where

import Statement
import Expression

functionCode :: String -> String -> [(String, String)] -> String -> String
functionCode "" "" [] "" = ""
functionCode "" "" (x:xs) "" = if xs == [] 
                               then fst x ++ " " ++ snd x
                               else fst x ++ " " ++ snd x ++ ", " ++ functionCode "" "" xs ""
functionCode n t a b = t ++ " " ++ n ++ "(" ++ functionCode "" "" a "" ++ ") const {\n" ++ b ++ "\n}\n"

classCode :: Expression RealSpace -> [String] -> String -> String
classCode e arg n = "class " ++ n ++ " : public FunctionalInterface {\npublic:\n" ++ n ++ codeA arg ++ "  {\n\thave_integral = true;\n}\n" ++
                functionCode "I_have_analytic_grad" "bool" [] "\treturn false;" ++
                functionCode "integral" "double" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] 
                    (unlines ["\tdouble output=0;",
                              "\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);",
                              codeStatements codeIntegrate ++ "\t// " ++ show (countFFT codeIntegrate) ++ " Fourier transform used.",
                              "\t// " ++ show (peakMem codeIntegrate),
                              "\treturn output;\n"]) ++
                functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] 
                    (unlines ["\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);",
                              "\tVectorXd output(gd.NxNyNz);",
                              codeStatements codeVTransform  ++ "\t// " ++ show (countFFT codeVTransform) ++ " Fourier transform used.",
                              "\t// " ++ show (peakMem codeVTransform), 
                              "\treturn output;\n"])  ++
                functionCode "transform" "double" [("double", "kT"), ("double", "x")] 
                    (unlines ["\tdouble output = 0;",
                              "\t" ++ codeStatements codeDTransform,
                              "\treturn output;\n"]) ++
                functionCode "append_to_name" "bool" [("", "std::string")] "\treturn false;" ++
                functionCode "derive" "double" [("double", "kT"), ("double", "x")] 
                    (unlines ["\tdouble output = 0;",
                              codeStatements codeDerive,
                              "\treturn output;\n"]) ++
                functionCode "d_by_dT" "double" [("double", ""), ("double", "")] "\tassert(0); // fail\n\treturn 0;\n" ++
                functionCode "derive_homogeneous" "Expression" [("const Expression &", "")] "\tassert(0); // fail\n\treturn Expression(0);\n" ++
                functionCode "grad" "Functional" [("const Functional", "&ingrad"), ("const Functional", "&x"), ("bool", "")] "\treturn ingrad;" ++
                functionCode "grad_T" "Functional" [("const Functional", "&ingradT")] "\treturn ingradT;" ++
                functionCode "grad" "void" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x"), ("const VectorXd", "&ingrad"), ("VectorXd", "*outgrad"), ("VectorXd", "*outpgrad")] ("\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);\n" ++ codeStatements codeGrad ++ "\t// " ++ show (countFFT codeGrad) ++ " Fourier transform used.\n\t// " ++ show (peakMem codeGrad) ++ "\n") ++
                functionCode "printme" "Expression" [("const Expression", "&x")] ("\treturn funexpr(\"" ++ n ++ "()\")(x);") ++
                functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")] "\tFunctionalInterface::print_summary(prefix, energy, name);" ++
                "private:\n"++ codeArgInit arg  ++"}; // End of " ++ n ++ " class\n\t// Total " ++ (show $ (countFFT codeIntegrate + countFFT codeVTransform + countFFT codeGrad)) ++ " Fourier transform used.\n\t// peak memory used: " ++ (show $ maximum $ map peakMem [codeIntegrate, codeVTransform, codeGrad])
    where
      codeIntegrate = reuseVar $ freeVectors (st ++ [AssignS (s_var "output") e'])
          where (st, e') = simp2 (joinFFTs $ integrate e)
      codeVTransform = reuseVar $ freeVectors (st ++ [AssignR (r_var "output") e'])
          where (st, e') = simp2 $ joinFFTs e
      codeDTransform = freeVectors (st ++ [AssignS (s_var "output") e'])
          where (st, e') = simp2 $ makeHomogeneous e
      codeDerive = freeVectors (st ++ [AssignS (s_var "output") e'])
          where (st, e') = simp2 $ derive (s_var "x") 1 $ makeHomogeneous e
      codeGrad = reuseVar $ freeVectors (st ++ [AssignR (r_var "(*outgrad)") (r_var "(*outgrad)" + e')])
          where (st, e') = simp2 (joinFFTs $ cleanvars $ derive (r_var "x") (r_var "ingrad") e)
      codeA [] = "()"
      codeA a = "(" ++ foldl1 (\x y -> x ++ ", " ++ y ) (map (\x -> "double " ++ x ++ "_arg") a) ++ ") : " ++ foldl1 (\x y -> x ++ ", " ++ y) (map (\x -> x ++ "(" ++ x ++ "_arg)") a)
      codeArgInit [] = ""
      codeArgInit a = unlines $ map (\x -> "\tdouble " ++ x ++ ";") a

generateHeader :: Expression RealSpace -> [String] -> String -> String
generateHeader e arg n = "// -*- mode: C++; -*-\n\n#include \"MinimalFunctionals.h\"\n#include \"utilities.h\"\n#include \"handymath.h\"\n\n" ++ 
                     classCode e arg (n ++ "_type") ++
                     "\n\nFunctional " ++ n ++"(" ++ codeA arg ++ ") {\n\treturn Functional(new " ++ n ++ "_type(" ++ codeA' arg ++ "), \"" ++ n ++ "\");\n}\n"
    where codeA [] = ""
          codeA a = foldl1 (\x y -> x ++ ", " ++ y ) (map ("double " ++) a)
          codeA' [] = ""
          codeA' a = foldl1 (\x y -> x ++ ", " ++ y ) a

\end{code}
