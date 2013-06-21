\begin{code}
{-# LANGUAGE PatternGuards #-}

module CodeGen (module Statement,
                module Expression,
                generateHeader, -- generateHeader is deprecated, and will soon be removed!
                defineFunctional, defineFunctionalNoGradient, defineTransformation )
    where

import Statement
import Expression
import qualified Data.Set as Set
import Optimize ( optimize )

functionCode :: String -> String -> [(String, String)] -> String -> String
functionCode "" "" [] "" = ""
functionCode "" "" (x:xs) "" = if xs == [] 
                               then fst x ++ " " ++ snd x
                               else fst x ++ " " ++ snd x ++ ", " ++ functionCode "" "" xs ""
functionCode n t a b = t ++ " " ++ n ++ "(" ++ functionCode "" "" a "" ++ ") const {\n" ++ b ++ "}\n"

classCode :: Expression RealSpace -> [String] -> String -> String
classCode ewithtransforms arg n = "class " ++ n ++ " : public FunctionalInterface {\npublic:\n" ++ n ++ codeA arg ++ "  {\n\thave_integral = true;" ++ definetransforms ++ "\n}\n" ++
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
                functionCode "grad" "Functional" [("const Functional", "&ingrad"), ("const Functional", "&x"), ("bool", "")] "\treturn ingrad;" ++
                functionCode "grad_T" "Functional" [("const Functional", "&ingradT")] "\treturn ingradT;" ++
                functionCode "grad" "void" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x"), ("const VectorXd", "&ingrad"), ("VectorXd", "*outgrad"), ("VectorXd", "*outpgrad")] (codeStatements codeGrad ++ "\t// " ++ show (countFFT codeGrad) ++ " Fourier transform used.\n\t// " ++ show (peakMem codeGrad) ++ "\n") ++
                functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")] "\tFunctionalInterface::print_summary(prefix, energy, name);" ++
                "private:\n"++ codeArgInit arg ++ declaretransforms ++"}; // End of " ++ n ++ " class\n\t// Total " ++ (show $ (countFFT codeIntegrate + countFFT codeVTransform + countFFT codeGrad)) ++ " Fourier transform used.\n\t// peak memory used: " ++ (show $ maximum $ map peakMem [codeIntegrate, codeVTransform, codeGrad])
    where
      defineGrid :: Type a => Expression a -> Expression a
      defineGrid = substitute dVscalar (s_var "gd.dvolume") .
                   substitute dr (s_var "gd.dvolume" ** (1.0/3))
      defineHomogeneousGrid = substitute dVscalar 1 .
                              substitute dr (s_var "gd.dvolume" ** (1.0/3))
      codeIntegrate = reuseVar $ freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = optimize [ES $ factorize $ joinFFTs $ cleanvars $
                                       defineGrid $ integrate e]
      codeVTransform = reuseVar $ freeVectors (st ++ [Assign (ER (r_var "output")) e'])
          where (st, [e']) = optimize [mkExprn $ factorize $ joinFFTs $ defineGrid e]
      codeDTransform = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = optimize [mkExprn $ makeHomogeneous e]
      codeDerive = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = optimize [ES $ derive (s_var "x") 1 $ makeHomogeneous e]
      codeGrad = reuseVar $ freeVectors (st ++ [Assign (ER (r_var "(*outgrad)"))
                                                       (ER (r_var "(*outgrad)" + e'))])
          where (st, [ER e']) = optimize [mkExprn $ factorize $ joinFFTs $ cleanvars $
                                          defineHomogeneousGrid $ derive (r_var "x") (r_var "ingrad") e]
      codeA [] = "()"
      codeA a = "(" ++ foldl1 (\x y -> x ++ ", " ++ y ) (map (\x -> "double " ++ x ++ "_arg") a) ++ ") : " ++ foldl1 (\x y -> x ++ ", " ++ y) (map (\x -> x ++ "(" ++ x ++ "_arg)") a)
      codeArgInit [] = ""
      codeArgInit a = unlines $ map (\x -> "\tdouble " ++ x ++ ";") a
      -- The following is needed to handle spherical fourier
      -- transforms that we store on 1D grids.
      declaretransforms = unlines $ map (\(t,_,_) -> "\tdouble *" ++ t ++ ";") transforms
      chomp str = case reverse str of '\n':rstr -> reverse rstr
                                      _ -> str
      definetransforms = chomp $ concat $ map definet transforms
        where definet (t,s,r) =
                unlines ["\t" ++ t ++ " = new double[" ++ nk ++ "];",
                         "\tfor (int i=0; i<" ++ nk ++"; i++) {",
                         "\t\tconst double k = i*" ++ show (dk s) ++ ";",
                         "\t\t" ++ t ++ "[i] = 0;",
                         "\t\tfor (double r="++halfdr++"; r<" ++ code (rmax s) ++"; r+=" ++ mydr ++ ") {",
                         "\t\t\tconst double rlo = r - " ++ halfdr ++ ";",
                         "\t\t\tconst double rhi = r + " ++ halfdr ++ ";",
                         "\t\t\t" ++ t ++ "[i] += " ++ code r ++ "*sin(k*r)*4*M_PI/3*(rhi*rhi*rhi-rlo*rlo*rlo);",
                         "\t\t}",
                         "I think this bit of code in CodeGen.lhs is broken and is never used.  If it is, this will cause the compiler to fail.",
                         "\t}"]
                  where nk = show (round (kmax s/dk s) :: Int)
                        mydr = code (rresolution s)
                        halfdr = code (rresolution s/2)
      (transforms, e) = mktransforms [1 :: Int ..] (findTransforms ewithtransforms) ewithtransforms
      mktransforms _ [] ee = ([], ee)
      mktransforms (nn:ns) (ft@(Expression (SphericalFourierTransform s r)):ts) ee =
        case mktransforms ns ts ee of
          (xs,ee') -> ((varname, s, r):xs, substitute ft evaluateme ee')
          where varname = "ft" ++ show nn
                evaluateme = setkzero (s_var (varname ++ "[0]")) $
                             s_var (varname ++ "[int(" ++ code k ++ "*" ++ show (1/dk s) ++ "+0.5)]")
      mktransforms _ (z:_) _ = error ("Not a transform: " ++ show z)

generateHeader :: Expression RealSpace -> [String] -> String -> String
generateHeader e arg n = "// -*- mode: C++; -*-\n\n#include \"MinimalFunctionals.h\"\n#include \"utilities.h\"\n#include \"handymath.h\"\n\n" ++
                     classCode e arg (n ++ "_type") ++
                     "\n\nFunctional " ++ n ++"(" ++ codeA arg ++ ") {\n\treturn Functional(new " ++ n ++ "_type(" ++ codeA' arg ++ "), \"" ++ n ++ "\");\n}\n"
    where codeA [] = ""
          codeA a = foldl1 (\x y -> x ++ ", " ++ y ) (map ("double " ++) a)
          codeA' [] = ""
          codeA' a = foldl1 (\x y -> x ++ ", " ++ y ) a




scalarClass :: Expression Scalar -> [String] -> String -> String
scalarClass ewithtransforms arg n =
  unlines
  ["class " ++ n ++ " : public FunctionalInterface {",
   "public:",
   "" ++ n ++ codeA arg ++ "  {",
   "\thave_integral = true;",
   "\t// TODO: code to evaluate Fourier transforms goes here",
   "\toldkT = sqrt(-1.0);  // initialize to NaN so we'll have to define transforms",
   -- insert code here to compute and store the spherical fourier transforms
   initializetransforms,
   "}",
   "~" ++ n ++ "() {",
   deletetransforms,
   "}",
   "",
   functionCode "I_have_analytic_grad" "bool" [] "\treturn false;",
   functionCode "integral" "double" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")]
    (unlines [redefinetransforms,
              "\tdouble output=0;",
              codeStatements codeIntegrate ++ "\t// " ++ show (countFFT codeIntegrate) ++ " Fourier transform used.",
              "\t// " ++ show (peakMem codeIntegrate) ++ " temporaries made",
              "\treturn output;", ""]),
   functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")]
    (unlines ["\tassert(0);"]),
   functionCode "transform" "double" [("double", "kT"), ("double", "x")]
    (unlines [redefinetransforms,
              "\tdouble output = 0;",
              "\t" ++ codeStatements codeDTransform,
              "\treturn output;", ""]),
   functionCode "append_to_name" "bool" [("", "std::string")] "\treturn false;",
   functionCode "derive" "double" [("double", "kT"), ("double", "x")]
    (unlines [redefinetransforms,
              "\tdouble output = 0;",
              codeStatements codeDerive,
              "\treturn output;", ""]),
   functionCode "d_by_dT" "double" [("double", ""), ("double", "")]
    (unlines ["\tassert(0); // fail", "\treturn 0;", ""]),
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
     (unlines [redefinetransforms,
               codeStatements codeGrad ++ "\t// " ++ show (countFFT codeGrad) ++ " Fourier transform used.",
               "\t// " ++ show (peakMem codeGrad),
               ""]),
   functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")]
    (unlines $ ["\tif (name != \"\") printf(\"%s%25s =\", prefix, name.c_str());",
                "\telse printf(\"%s%25s =\", prefix, \"UNKNOWN\");",
                "\tprint_double(\"\", energy);"] ++
               map printEnergy (Set.toList (findNamedScalars e)) ++
               ["\tprintf(\"\\n\");"]),
  "private:",
  ""++ codeArgInit arg ++ codeMutableData (Set.toList $ findNamedScalars e)
  ++ "\t// TODO: add declaration of spherical fourier transform data here\n"
  ++ declaretransforms
  ++"}; // End of " ++ n ++ " class",
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
          where (st0, [e']) = optimize [mkExprn $ factorize $ joinFFTs $ cleanvars $ defineGrid e]
                st = filter (not . isns) st0
                isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
                isns _ = False
                ns = findNamedScalars e
      codeDTransform = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = optimize [mkExprn $ makeHomogeneous e]
      codeDerive = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = optimize [ES $ derive (s_var "x") 1 $ makeHomogeneous e]
      codeGrad = reuseVar $ freeVectors (st ++ [Assign (ER (r_var "(*outgrad)"))
                                                       (ER (r_var "(*outgrad)" + e'))])
          where (st, [ER e']) = optimize [mkExprn $ factorize $ joinFFTs $ cleanvars $
                                          defineHomogeneousGrid $
                                          (r_var "ingrad" * derive (r_var "x") 1 e)]
      codeA [] = "()"
      codeA a = "(" ++ foldl1 (\x y -> x ++ ", " ++ y ) (map (\x -> "double " ++ x ++ "_arg") a) ++ ") : " ++ foldl1 (\x y -> x ++ ", " ++ y) (map (\x -> x ++ "(" ++ x ++ "_arg)") a)
      codeArgInit a = unlines $ map (\x -> "\tdouble " ++ x ++ ";") a
      codeMutableData a = unlines $ map (\x -> "\tmutable double " ++ x ++ ";") a
      -- The following is needed to handle spherical fourier
      -- transforms that we store on 1D grids.
      declaretransforms = chomp $ unlines $ "\tmutable double oldkT;" : map (\(t,_,_) -> "\tmutable double *" ++ t ++ ";") transforms
      chomp str = case reverse str of '\n':rstr -> reverse rstr
                                      _ -> str
      initializetransforms = chomp $ unlines $ map initt transforms
        where initt (t,_,_) = "\t" ++ t ++ " = 0;"
      deletetransforms = chomp $ unlines $ map delt transforms
        where delt (t,_,_) = "\tdelete[] " ++ t ++ ";"
      redefinetransforms = chomp $ unlines $ [
                         "\tif (oldkT != kT) {",
                         "\t\toldkT = kT;",
                         chomp $ concat $ map definet transforms,
                         "\t}"]
        where definet (t,s@(Spherical {}),r) =
                unlines ["\t\tif (!"++t++") " ++ t ++ " = new double[" ++ nk ++ "];",
                         "\t\t// first evaluate k=0 version of " ++ t ++ "...",
                         "\t\t" ++ t ++ "[0] = 0;",
                         "\t\tfor (double r="++halfdr++"; r<" ++ code (rmax s) ++"; r+=" ++ mydr ++ ") {",
                         "\t\t\t\t" ++ t ++ "[0] += " ++ code (r*dvolume) ++ ";",
                         "\t\t}",
                         "\t\t// ... now work on all other k values of " ++ t,
                         "\t\tfor (int i=1; i<" ++ nk ++"; i++) {",
                         "\t\t\tconst double k = i*" ++ show (dk s) ++ ";",
                         "\t\t\t" ++ t ++ "[i] = 0;",
                         "\t\t\tfor (double r="++halfdr++"; r<" ++ code (rmax s) ++"; r+=" ++ mydr ++ ") {",
                         "\t\t\t\t" ++ t ++ "[i] += " ++ code (r*sin(kvar*rvar)/(kvar*rvar)*dvolume)  ++ ";",
                         "\t\t\t}",
                         "\t\t}"]
                  where nk = show (round (kmax s/dk s) :: Int)
                        mydr = code (rresolution s)
                        halfdr = code (rresolution s/2)
                        kvar = s_var "k"
                        rvar = s_var "r"
                        dvolume = "dvolume" === 4*pi/3*(rhi**3 - rlo**3)
                        rhi = "rhi" === rvar + rresolution s/2
                        rlo = "rlo" === rvar - rresolution s/2
              definet (t,s@(VectorS {}),r) =
                unlines ["\t\tif (!"++t++") " ++ t ++ " = new double[" ++ nk ++ "];",
                         "\t\tfor (int i=1; i<" ++ nk ++"; i++) {",
                         "\t\t\tconst double k = i*" ++ show (dk s) ++ ";",
                         "\t\t\t" ++ t ++ "[i] = 0;",
                         "\t\t\tfor (double r="++halfdr++"; r<" ++ code (rmax s) ++"; r+=" ++ mydr ++ ") {",
                         "\t\t\t\t" ++ t ++ "[i] += " ++ code (r*(cos kr - sin kr/kr)*dvolume/kvar**2) ++ ";",
                         "\t\t\t}",
                         "\t\t}"]
                  where nk = show (round (kmax s/dk s) :: Int)
                        mydr = code (rresolution s)
                        halfdr = code (rresolution s/2)
                        kr = kvar*rvar
                        kvar = s_var "k"
                        rvar = s_var "r"
                        dvolume = "dvolume" === 4*pi/3*(rhi**3 - rlo**3)
                        rhi = "rhi" === rvar + rresolution s/2
                        rlo = "rlo" === rvar - rresolution s/2
      (transforms, e) = mktransforms [1 :: Int ..] (findTransforms ewithtransforms) ewithtransforms
      mktransforms _ [] ee = ([], ee)
      mktransforms (nn:ns) (ft@(Expression (SphericalFourierTransform s r)):ts) ee =
        case mktransforms ns ts ee of
          (xs,ee') -> ((varname, s, r):xs, substitute ft evaluateme ee')
          where varname = "ft" ++ show nn
                evaluateme = setkzero (s_var (varname ++ "[0]")) $
                             s_var (varname ++ "[int(" ++ code k ++ "*" ++ show (1/dk s) ++ "+0.5)]")
      mktransforms _ (z:_) _ = error ("Not a transform: " ++ show z)

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
scalarClassNoGradient ewithtransforms arg n =
  unlines
  ["class " ++ n ++ " : public FunctionalInterface {",
   "public:",
   "" ++ n ++ codeA arg ++ "  {",
   "\thave_integral = true;",
   definetransforms,
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
   functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")]
    (unlines $ ["\tif (name != \"\") printf(\"%s%25s =\", prefix, name.c_str());",
                "\telse printf(\"%s%25s =\", prefix, \"UNKNOWN\");",
                "\tprint_double(\"\", energy);"] ++
               map printEnergy (Set.toList (findNamedScalars e)) ++
               ["\tprintf(\"\\n\");"]),
  "private:",
  ""++ codeArgInit arg ++ codeMutableData (Set.toList $ findNamedScalars e)
    ++ "\t// TODO: add declaration of spherical fourier transform data here\n"
    ++ declaretransforms
  ++"}; // End of " ++ n ++ " class",
  "\t// Total " ++ (show $ countFFT codeIntegrate) ++ " Fourier transform used.",
  "\t// peak memory used: " ++ (show $ maximum $ map peakMem [codeIntegrate])
  ]
    where
      defineGrid = substitute dVscalar (s_var "gd.dvolume") .
                   substitute dr (s_var "gd.dvolume" ** (1.0/3))
      printEnergy v = "\tprintf(\"\\n%s%25s =\", prefix, \"" ++ v ++ "\");\n" ++
                      "\tprint_double(\"\", " ++ v ++ ");"
      codeIntegrate = reuseVar $ freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st0, [e']) = optimize [mkExprn $ factorize $ joinFFTs $ cleanvars $ defineGrid e]
                st = filter (not . isns) st0
                isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
                isns _ = False
                ns = findNamedScalars e
      codeDTransform = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = optimize [mkExprn $ makeHomogeneous e]
      codeDerive = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = optimize [ES $ derive (s_var "x") 1 $ makeHomogeneous e]
      codeA [] = "()"
      codeA a = "(" ++ foldl1 (\x y -> x ++ ", " ++ y ) (map (\x -> "double " ++ x ++ "_arg") a) ++ ") : " ++ foldl1 (\x y -> x ++ ", " ++ y) (map (\x -> x ++ "(" ++ x ++ "_arg)") a)
      codeArgInit a = unlines $ map (\x -> "\tdouble " ++ x ++ ";") a
      codeMutableData a = unlines $ map (\x -> "\tmutable double " ++ x ++ ";") a
      -- The following is needed to handle spherical fourier
      -- transforms that we store on 1D grids.
      declaretransforms = unlines $ map (\(t,_,_) -> "\tdouble *" ++ t ++ ";") transforms
      chomp str = case reverse str of '\n':rstr -> reverse rstr
                                      _ -> str
      definetransforms = chomp $ concat $ map definet transforms
        where definet (t,s,r) =
                unlines ["\t" ++ t ++ " = new double[" ++ nk ++ "];",
                         "\tfor (int i=0; i<" ++ nk ++"; i++) {",
                         "\t\tconst double k = i*" ++ show (dk s) ++ ";",
                         "\t\t" ++ t ++ "[i] = 0;",
                         "\t\tfor (double r="++halfdr++"; r<" ++ code (rmax s) ++"; r+=" ++ mydr ++ ") {",
                         "\t\t\tconst double rlo = r - " ++ halfdr ++ ";",
                         "\t\t\tconst double rhi = r + " ++ halfdr ++ ";",
                         "\t\t\t" ++ t ++ "[i] += " ++ code r ++ "*sin(k*r)*4*M_PI/3*(rhi*rhi*rhi-rlo*rlo*rlo);",
                         "\t\t}",
                         "I think this bit of old code in CodeGen.lhs is broken and is never used.  If it is, this will cause the compiler to fail.",
                         "\t}"]
                  where nk = show (round (kmax s/dk s) :: Int)
                        mydr = code (rresolution s)
                        halfdr = code (rresolution s/2)
      (transforms, e) = mktransforms [1 :: Int ..] (findTransforms ewithtransforms) ewithtransforms
      mktransforms _ [] ee = ([], ee)
      mktransforms (nn:ns) (ft@(Expression (SphericalFourierTransform s r)):ts) ee =
        case mktransforms ns ts ee of
          (xs,ee') -> ((varname, s, r):xs, substitute ft evaluateme ee')
          where varname = "ft" ++ show nn
                evaluateme = setkzero (s_var (varname ++ "[0]")) $
                             s_var (varname ++ "[int(" ++ code k ++ "*" ++ show (1/dk s) ++ "+0.5)]")
      mktransforms _ (z:_) _ = error ("Not a transform: " ++ show z)

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




transformationClass :: Expression RealSpace -> [String] -> String -> String
transformationClass ewithtransforms arg n =
  unlines
  ["class " ++ n ++ " : public FunctionalInterface {",
   "public:",
   "" ++ n ++ codeA arg ++ "  {",
   "\thave_integral = true;",
   definetransforms,
   "}",
   "",
   functionCode "I_have_analytic_grad" "bool" [] "\treturn false;",
   functionCode "integral" "double" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")]
    (unlines ["\tassert(0);"]),
   functionCode "transform" "VectorXd" [("const GridDescription", "&gd"), ("double", "kT"), ("const VectorXd", "&x")] 
   (unlines ["\tVectorXd output(gd.NxNyNz);",
             codeStatements codeVTransform  ++ "\t// " ++ show (countFFT codeVTransform) ++ " Fourier transform used.",
             "\t// " ++ show (peakMem codeVTransform),
             "\treturn output;\n"]),
   functionCode "transform" "double" [("double", "kT"), ("double", "x")]
    (unlines ["\tdouble output = 0;",
              "\t" ++ codeStatements codeDTransform,
              "\treturn output;", ""]),
   functionCode "append_to_name" "bool" [("", "std::string")] "\treturn false;",
   functionCode "derive" "double" [("double", "kT"), ("double", "x")]
    (unlines ["\tassert(0);"]),
   functionCode "d_by_dT" "double" [("double", ""), ("double", "")]
    (unlines ["\tassert(0); // fail", "\treturn 0;", ""]),
   functionCode "grad" "Functional" [("const Functional", "&ingrad"), ("const Functional", "&x"), ("bool", "")]
    "\tassert(0);",
   functionCode "grad_T" "Functional" [("const Functional", "&ingradT")] "\tassert(0);",
   functionCode "grad" "void"
     [("const GridDescription", "&gd"),
      ("double", "kT"),
      ("const VectorXd", "&x"),
      ("const VectorXd", "&ingrad"),
      ("VectorXd", "*outgrad"),
      ("VectorXd", "*outpgrad")]
     (unlines ["\tassert(0); // fail"]),
   functionCode "print_summary" "void" [("const char", "*prefix"), ("double", "energy"), ("std::string", "name")]
    (unlines $ ["\tprintf(\"Nothing to say\");"]),
  "private:",
  ""++ codeArgInit arg ++ codeMutableData (Set.toList $ findNamedScalars e)
    ++ "\t// TODO: add declaration of spherical fourier transform data here\n"
    ++ declaretransforms
  ++"}; // End of " ++ n ++ " class"
  ]
    where
      defineGrid = substitute dVscalar (s_var "gd.dvolume") .
                   substitute dr (s_var "gd.dvolume" ** (1.0/3))
      codeVTransform = reuseVar $ freeVectors (st ++ [Assign (ER (r_var "output")) e'])
          where (st, [e']) = optimize [mkExprn $ factorize $ joinFFTs $ defineGrid e]
      codeDTransform = freeVectors (st ++ [Assign (ES (s_var "output")) e'])
          where (st, [e']) = optimize [mkExprn $ makeHomogeneous e]
      codeA [] = "()"
      codeA a = "(" ++ foldl1 (\x y -> x ++ ", " ++ y ) (map (\x -> "double " ++ x ++ "_arg") a) ++ ") : " ++ foldl1 (\x y -> x ++ ", " ++ y) (map (\x -> x ++ "(" ++ x ++ "_arg)") a)
      codeArgInit a = unlines $ map (\x -> "\tdouble " ++ x ++ ";") a
      codeMutableData a = unlines $ map (\x -> "\tmutable double " ++ x ++ ";") a
      -- The following is needed to handle spherical fourier
      -- transforms that we store on 1D grids.
      declaretransforms = unlines $ map (\(t,_,_) -> "\tdouble *" ++ t ++ ";") transforms
      chomp str = case reverse str of '\n':rstr -> reverse rstr
                                      _ -> str
      definetransforms = chomp $ concat $ map definet transforms
        where definet (t,s,r) =
                unlines ["\t" ++ t ++ " = new double[" ++ nk ++ "];",
                         "\tfor (int i=0; i<" ++ nk ++"; i++) {",
                         "\t\tconst double k = i*" ++ show (dk s) ++ ";",
                         "\t\t" ++ t ++ "[i] = 0;",
                         "\t\tfor (double r="++halfdr++"; r<" ++ code (rmax s) ++"; r+=" ++ mydr ++ ") {",
                         "\t\t\tconst double rlo = r - " ++ halfdr ++ ";",
                         "\t\t\tconst double rhi = r + " ++ halfdr ++ ";",
                         "\t\t\t" ++ t ++ "[i] += " ++ code r ++ "*sin(k*r)*4*M_PI/3*(rhi*rhi*rhi-rlo*rlo*rlo);",
                         "\t\t}",
                         "I think this bit of old code in CodeGen.lhs is broken and is never used.  If it is, this will cause the compiler to fail.",
                         "\t}"]
                  where nk = show (round (kmax s/dk s) :: Int)
                        mydr = code (rresolution s)
                        halfdr = code (rresolution s/2)
      (transforms, e) = mktransforms [1 :: Int ..] (findTransforms ewithtransforms) ewithtransforms
      mktransforms _ [] ee = ([], ee)
      mktransforms (nn:ns) (ft@(Expression (SphericalFourierTransform s r)):ts) ee =
        case mktransforms ns ts ee of
          (xs,ee') -> ((varname, s, r):xs, substitute ft evaluateme ee')
          where varname = "ft" ++ show nn
                evaluateme = setkzero (s_var (varname ++ "[0]")) $
                             s_var (varname ++ "[int(" ++ code k ++ "*" ++ show (1/dk s) ++ "+0.5)]")
      mktransforms _ (z:_) _ = error ("Not a transform: " ++ show z)

defineTransformation :: Expression RealSpace -> [String] -> String -> String
defineTransformation e arg n =
  unlines ["// -*- mode: C++; -*-",
           "",
           "#include \"MinimalFunctionals.h\"",
           "#include \"utilities.h\"",
           "#include \"handymath.h\"",
           "",
           "",
           transformationClass e arg (n ++ "_type"),
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
