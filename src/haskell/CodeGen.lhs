\begin{code}
{-# LANGUAGE PatternGuards #-}

module CodeGen (module Statement,
                module Expression,
                generateHeader, defineFunctional )
    where

import Statement
import Expression
import qualified Data.Set as Set

functionCode :: String -> String -> [(String, String)] -> String -> String
functionCode "" "" [] "" = ""
functionCode "" "" (x:xs) "" = if xs == [] 
                               then fst x ++ " " ++ snd x
                               else fst x ++ " " ++ snd x ++ ", " ++ functionCode "" "" xs ""
functionCode n t a b = t ++ " " ++ n ++ "(" ++ functionCode "" "" a "" ++ ") const {\n" ++ b ++ "}\n"

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
      codeIntegrate = reuseVar $ freeVectors (st ++ [Assign (ES (s_var "output")) (ES e')])
          where (st, e') = simp2 (factorize $ joinFFTs $ integrate e)
      codeVTransform = reuseVar $ freeVectors (st ++ [Assign (ER (r_var "output")) (ER e')])
          where (st, e') = simp2 $ factorize $ joinFFTs e
      codeDTransform = freeVectors (st ++ [Assign (ES (s_var "output")) (ES e')])
          where (st, e') = simp2 $ makeHomogeneous e
      codeDerive = freeVectors (st ++ [Assign (ES (s_var "output")) (ES e')])
          where (st, e') = simp2 $ derive (s_var "x") 1 $ makeHomogeneous e
      codeGrad = reuseVar $ freeVectors (st ++ [Assign (ER (r_var "(*outgrad)"))
                                                       (ER (r_var "(*outgrad)" + e'))])
          where (st, e') = simp2 (factorize $ joinFFTs $ cleanvars $ derive (r_var "x") (r_var "ingrad") e)
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




scalarClass :: Expression Scalar -> [String] -> String -> String
scalarClass e arg n =
  unlines
  ["class " ++ n ++ " : public FunctionalInterface {",
   "public:",
   "" ++ n ++ codeA arg ++ "  {",
   "\thave_integral = true;",
   "}",
   "",
   functionCode "I_have_analytic_grad" "bool" [] "\treturn false;",
   functionCode "integral" "double" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")]
    (unlines ["\tdouble output=0;",
              "\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);",
              codeStatements codeIntegrate ++ "\t// " ++ show (countFFT codeIntegrate) ++ " Fourier transform used.",
              "\t// " ++ show (peakMem codeIntegrate) ++ " temporaries made",
              "\treturn output;", ""]),
   functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")]
    (unlines ["\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);",
              "\tassert(0);"]),
   functionCode "transform" "double" [("double", "kT"), ("double", "x")]
    (unlines ["\tdouble output = 0;",
              "\t" ++ codeStatements codeDTransform,
              "\treturn output;", ""]),
   functionCode "append_to_name" "bool" [("", "std::string")] "\treturn false;",
   functionCode "derive" "double" [("double", "kT"), ("double", "x")]
    (unlines ["\tdouble output = 0;",
              codeStatements codeDerive,
              "\treturn output;", ""]),
   functionCode "d_by_dT" "double" [("double", ""), ("double", "")]
    (unlines ["\tassert(0); // fail", "\treturn 0;", ""]),
   functionCode "derive_homogeneous" "Expression" [("const Expression &", "")]
    (unlines ["\tassert(0); // fail", "\treturn Expression(0);", ""]),
   functionCode "grad" "Functional" [("const Functional", "&ingrad"), ("const Functional", "&x"), ("bool", "")]
    "\treturn ingrad;",
   functionCode "grad_T" "Functional" [("const Functional", "&ingradT")] "\treturn ingradT;",
   functionCode "grad" "void"
     [("const GridDescription", "&gd"),
      ("double", "kT"),
      ("const VectorXd", "&x"),
      ("const VectorXd", "&ingrad"),
      ("VectorXd", "*outgrad"),
      ("VectorXd", "*outpgrad")]
     (unlines ["\tconst double dr = pow(gd.fineLat.volume(), 1.0/3);",
               codeStatements codeGrad ++ "\t// " ++ show (countFFT codeGrad) ++ " Fourier transform used.",
               "\t// " ++ show (peakMem codeGrad),
               ""]),
   functionCode "printme" "Expression" [("const Expression", "&x")] ("\treturn funexpr(\"" ++ n ++ "()\")(x);"),
   functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")]
    (unlines $ ["\tif (name != \"\") printf(\"%s%25s =\", prefix, name.c_str());",
                "\telse printf(\"%s%25s =\", prefix, \"UNKNOWN\");",
                "\tprint_double(\"\", energy);"] ++
               map printEnergy (Set.toList (findNamedScalars e)) ++
               ["\tprintf(\"\\n\");"]),
  "private:",
  ""++ codeArgInit arg ++ codeMutableData (Set.toList $ findNamedScalars e)  ++"}; // End of " ++ n ++ " class",
  "\t// Total " ++ (show $ (countFFT codeIntegrate + countFFT codeGrad)) ++ " Fourier transform used.",
  "\t// peak memory used: " ++ (show $ maximum $ map peakMem [codeIntegrate, codeGrad])
  ]
    where
      printEnergy v = "\tprintf(\"\\n%s%25s =\", prefix, \"" ++ v ++ "\");\n" ++
                      "\tprint_double(\"\", " ++ v ++ ");"
      codeIntegrate = reuseVar $ freeVectors (st ++ [Assign (ES (s_var "output")) (ES e')])
          where (st0, e') = simp2 (factorize $ joinFFTs e)
                st = filter (not . isns) st0
                isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
                isns _ = False
                ns = findNamedScalars e
      codeDTransform = freeVectors (st ++ [Assign (ES (s_var "output")) (ES e')])
          where (st, e') = simp2 $ makeHomogeneous e
      codeDerive = freeVectors (st ++ [Assign (ES (s_var "output")) (ES e')])
          where (st, e') = simp2 $ derive (s_var "x") 1 $ makeHomogeneous e
      codeGrad = reuseVar $ freeVectors (st ++ [Assign (ER (r_var "(*outgrad)"))
                                                       (ER (r_var "(*outgrad)" + e'))])
          where (st, e') = simp2 $ factorize $ joinFFTs $ cleanvars $
                           substitute (s_var "dV" :: Expression RealSpace) 1 $
                           (r_var "ingrad" * derive (r_var "x") 1 e)
      codeA [] = "()"
      codeA a = "(" ++ foldl1 (\x y -> x ++ ", " ++ y ) (map (\x -> "double " ++ x ++ "_arg") a) ++ ") : " ++ foldl1 (\x y -> x ++ ", " ++ y) (map (\x -> x ++ "(" ++ x ++ "_arg)") a)
      codeArgInit a = unlines $ map (\x -> "\tdouble " ++ x ++ ";") a
      codeMutableData a = unlines $ map (\x -> "\tmutable double " ++ x ++ ";") a

defineFunctional :: Expression Scalar -> [String] -> String -> String
defineFunctional e arg n =
  unlines ["// -*- mode: C++; -*-",
           "",
           "#include \"MinimalFunctionals.h\"",
           "#include \"utilities.h\"",
           "#include \"handymath.h\"",
           "",
           "",
           scalarClass e arg (n ++ "_type"),
           "",
           "",
           "Functional " ++ n ++"(" ++ codeA arg ++ ") {",
           "\treturn Functional(new " ++ n ++ "_type(" ++ codeA' arg ++ "), \"" ++ n ++ "\");",
           "}",
           ""]
    where codeA [] = ""
          codeA a = foldl1 (\x y -> x ++ ", " ++ y ) (map ("double " ++) a)
          codeA' [] = ""
          codeA' a = foldl1 (\x y -> x ++ ", " ++ y ) a

findNamedScalars :: Type b => Expression b -> Set.Set String
findNamedScalars e@(Expression _)
    | EK (Expression (FFT e')) <- mkExprn e = findNamedScalars e'
    | ER (Expression (IFFT e')) <- mkExprn e = findNamedScalars e'
    | ES (Expression (Integrate e')) <- mkExprn e = findNamedScalars e'
    | otherwise = Set.empty
findNamedScalars (Sum s _) = Set.unions $ map (findNamedScalars . snd) $ sum2pairs s
findNamedScalars (Product p _) = Set.unions $ map (findNamedScalars . fst) $ product2pairs p
findNamedScalars (Cos e) = findNamedScalars e
findNamedScalars (Sin e) = findNamedScalars e
findNamedScalars (Log e) = findNamedScalars e
findNamedScalars (Exp e) = findNamedScalars e
findNamedScalars (Abs e) = findNamedScalars e
findNamedScalars (Signum e) = findNamedScalars e
findNamedScalars (Var _ _ b _ (Just e))
  | ES _ <- mkExprn e = Set.insert b $ findNamedScalars e
  | otherwise = findNamedScalars e
findNamedScalars (Var _ _ _ _ Nothing) = Set.empty
findNamedScalars (Scalar e) = findNamedScalars e
\end{code}
