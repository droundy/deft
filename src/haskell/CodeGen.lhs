\begin{code}
{-# LANGUAGE PatternGuards #-}

module CodeGen (module Statement,
                module Expression,
                generateHeader, defineFunctional, defineFunctionalNoGradient )
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
                              codeStatements codeIntegrate ++ "\t// " ++ show (countFFT codeIntegrate) ++ " Fourier transform used.",
                              "\t// " ++ show (peakMem codeIntegrate),
                              "\treturn output;\n"]) ++
                functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] 
                    (unlines ["\tVectorXd output(gd.NxNyNz);",
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
                functionCode "grad" "void" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x"), ("const VectorXd", "&ingrad"), ("VectorXd", "*outgrad"), ("VectorXd", "*outpgrad")] (codeStatements codeGrad ++ "\t// " ++ show (countFFT codeGrad) ++ " Fourier transform used.\n\t// " ++ show (peakMem codeGrad) ++ "\n") ++
                functionCode "printme" "Expression" [("const Expression", "&x")] ("\treturn funexpr(\"" ++ n ++ "()\")(x);") ++
                functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")] "\tFunctionalInterface::print_summary(prefix, energy, name);" ++
                "private:\n"++ codeArgInit arg  ++"}; // End of " ++ n ++ " class\n\t// Total " ++ (show $ (countFFT codeIntegrate + countFFT codeVTransform + countFFT codeGrad)) ++ " Fourier transform used.\n\t// peak memory used: " ++ (show $ maximum $ map peakMem [codeIntegrate, codeVTransform, codeGrad])
    where
      defineGrid :: Type a => Expression a -> Expression a
      defineGrid = substitute dVscalar (s_var "gd.dvolume") .
                   substitute dr (s_var "gd.dvolume" ** (1.0/3))
      defineHomogeneousGrid = substitute dVscalar 1 .
                              substitute dr (s_var "gd.dvolume" ** (1.0/3))
      codeIntegrate = reuseVar $ freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = simp2 [ES $ factorize $ joinFFTs $ cleanvars $ defineGrid $ integrate e]
      codeVTransform = reuseVar $ freeVectors (st ++ [Assign (ER (r_var "output")) e'])
          where (st, [e']) = simp2 [mkExprn $ factorize $ joinFFTs $ defineGrid e]
      codeDTransform = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = simp2 [mkExprn $ makeHomogeneous e]
      codeDerive = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = simp2 [ES $ derive (s_var "x") 1 $ makeHomogeneous e]
      codeGrad = reuseVar $ freeVectors (st ++ [Assign (ER (r_var "(*outgrad)"))
                                                       (ER (r_var "(*outgrad)" + e'))])
          where (st, [ER e']) = simp2 [mkExprn $ factorize $ joinFFTs $ cleanvars $
                                       defineHomogeneousGrid $ derive (r_var "x") (r_var "ingrad") e]
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
              codeStatements codeIntegrate ++ "\t// " ++ show (countFFT codeIntegrate) ++ " Fourier transform used.",
              "\t// " ++ show (peakMem codeIntegrate) ++ " temporaries made",
              "\treturn output;", ""]),
   functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")]
    (unlines ["\tassert(0);"]),
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
     (unlines [codeStatements codeGrad ++ "\t// " ++ show (countFFT codeGrad) ++ " Fourier transform used.",
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
      defineGrid = substitute dVscalar (s_var "gd.dvolume") .
                   substitute dr (s_var "gd.dvolume" ** (1.0/3))
      defineHomogeneousGrid = substitute dVscalar 1 .
                              substitute dr (s_var "gd.dvolume" ** (1.0/3))
      printEnergy v = "\tprintf(\"\\n%s%25s =\", prefix, \"" ++ v ++ "\");\n" ++
                      "\tprint_double(\"\", " ++ v ++ ");"
      codeIntegrate = reuseVar $ freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st0, [e']) = simp2 [mkExprn $ factorize $ joinFFTs $ cleanvars $ defineGrid e]
                st = filter (not . isns) st0
                isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
                isns _ = False
                ns = findNamedScalars e
      codeDTransform = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = simp2 [mkExprn $ makeHomogeneous e]
      codeDerive = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = simp2 [ES $ derive (s_var "x") 1 $ makeHomogeneous e]
      codeGrad = reuseVar $ freeVectors (st ++ [Assign (ER (r_var "(*outgrad)"))
                                                       (ER (r_var "(*outgrad)" + e'))])
          where (st, [ER e']) = simp2 [mkExprn $ factorize $ joinFFTs $ cleanvars $
                                       defineHomogeneousGrid $
                                       (r_var "ingrad" * derive (r_var "x") 1 e)]
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





scalarClassNoGradient :: Expression Scalar -> [String] -> String -> String
scalarClassNoGradient e arg n =
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
              codeStatements codeIntegrate ++ "\t// " ++ show (countFFT codeIntegrate) ++ " Fourier transform used.",
              "\t// " ++ show (peakMem codeIntegrate) ++ " temporaries made",
              "\treturn output;", ""]),
   functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")]
    (unlines ["\tassert(0);"]),
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
     (unlines ["\tassert(0); // fail"]),
   functionCode "printme" "Expression" [("const Expression", "&x")] ("\treturn funexpr(\"" ++ n ++ "()\")(x);"),
   functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")]
    (unlines $ ["\tif (name != \"\") printf(\"%s%25s =\", prefix, name.c_str());",
                "\telse printf(\"%s%25s =\", prefix, \"UNKNOWN\");",
                "\tprint_double(\"\", energy);"] ++
               map printEnergy (Set.toList (findNamedScalars e)) ++
               ["\tprintf(\"\\n\");"]),
  "private:",
  ""++ codeArgInit arg ++ codeMutableData (Set.toList $ findNamedScalars e)  ++"}; // End of " ++ n ++ " class",
  "\t// Total " ++ (show $ countFFT codeIntegrate) ++ " Fourier transform used.",
  "\t// peak memory used: " ++ (show $ maximum $ map peakMem [codeIntegrate])
  ]
    where
      defineGrid = substitute dVscalar (s_var "gd.dvolume") .
                   substitute dr (s_var "gd.dvolume" ** (1.0/3))
      printEnergy v = "\tprintf(\"\\n%s%25s =\", prefix, \"" ++ v ++ "\");\n" ++
                      "\tprint_double(\"\", " ++ v ++ ");"
      codeIntegrate = reuseVar $ freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st0, [e']) = simp2 [mkExprn $ factorize $ joinFFTs $ cleanvars $ defineGrid e]
                st = filter (not . isns) st0
                isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
                isns _ = False
                ns = findNamedScalars e
      codeDTransform = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = simp2 [mkExprn $ makeHomogeneous e]
      codeDerive = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = simp2 [ES $ derive (s_var "x") 1 $ makeHomogeneous e]
      codeA [] = "()"
      codeA a = "(" ++ foldl1 (\x y -> x ++ ", " ++ y ) (map (\x -> "double " ++ x ++ "_arg") a) ++ ") : " ++ foldl1 (\x y -> x ++ ", " ++ y) (map (\x -> x ++ "(" ++ x ++ "_arg)") a)
      codeArgInit a = unlines $ map (\x -> "\tdouble " ++ x ++ ";") a
      codeMutableData a = unlines $ map (\x -> "\tmutable double " ++ x ++ ";") a

defineFunctionalNoGradient :: Expression Scalar -> [String] -> String -> String
defineFunctionalNoGradient e arg n =
  unlines ["// -*- mode: C++; -*-",
           "",
           "#include \"MinimalFunctionals.h\"",
           "#include \"utilities.h\"",
           "#include \"handymath.h\"",
           "",
           "",
           scalarClassNoGradient e arg (n ++ "_type"),
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
\end{code}
