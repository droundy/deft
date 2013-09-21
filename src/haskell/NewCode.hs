{-# LANGUAGE PatternGuards #-}

module NewCode (module Statement,
                module Expression,
                defineFunctional, createCppFile, createHeader )
    where

import C
import Statement
import Expression
import Optimize ( optimize )
import qualified Data.Set as Set


defineFunctional :: Expression Scalar -> [Exprn] -> String -> String
defineFunctional e variables n = createHeader e n ++ implementation
  where implementation = unlines $ drop 4 $ lines $ createCppFile e variables n ""


functionDeclaration :: String -> String -> [(String, String)] -> String
functionDeclaration "" "" [] = ""
functionDeclaration "" "" (x:xs) =
  if xs == []
  then fst x ++ " " ++ snd x
  else fst x ++ " " ++ snd x ++ ", " ++ functionDeclaration "" "" xs
functionDeclaration n t a =
  "\t" ++ t ++ " " ++ n ++ "(" ++ functionDeclaration "" "" a ++ ") const;"


createHeader :: Expression Scalar -> String -> String
createHeader e = if any is_nonscalar (findOrderedInputs e)
                 then create3dHeader e
                 else create0dHeader e
  where is_nonscalar (ES _) = False
        is_nonscalar _ = True

create3dHeader :: Expression Scalar -> String -> String
create3dHeader e0 n =
  unlines $
  ["// -*- mode: C++; -*-",
   "",
   "#include \"new/NewFunctional.h\"",
   "#include \"utilities.h\"",
   "#include \"handymath.h\"",
   "",
   "",
   "class " ++ n ++ " : public NewFunctional {",
   "public:",
   "\t" ++ n ++ "(int Nx, int Ny, int Nz);",
   "",
   functionDeclaration "true_energy" "double" [],
   functionDeclaration "grad" "Vector" [],
   functionDeclaration "printme" "void" [("const char *", "prefix")]]++
  map setarg (inputs e0) ++
 ["private:",
  ""++ codeMutableData (Set.toList $ findNamedScalars e)  ++"}; // End of " ++ n ++ " class"]
    where
      e = mapExpression' renameVar e0
      renameVar (Var CannotBeFreed a b c x) = Var CannotBeFreed ("var"++a) ("var"++b) c x
      renameVar x = x
      -- setarg creates a method that will get a reference to a given
      -- input argument's value.
      setarg :: Exprn -> String
      setarg ee =
        newcode $ CFunction {
          name = nameE ee,
          returnType = Reference (ctype ee),
          constness = "const",
          args = [],
          contents = "\tint sofar = 0;" : initme (inputs e0)
          }
        where initme (xx@(ES _):_) | xx == ee = ["\treturn data[sofar];"]
              initme (xx:_) | xx == ee = ["\treturn data.slice(sofar," ++ sizeE ee ++ ");"]
              initme (xx@(ES _):rr)
                | nameE xx `elem` ["Nx","Ny","Nz"] =
                  ("\tconst double "++nameE xx++" = data[sofar++];") : initme rr
              initme (xx@(ES _):rr) = ("\tsofar += 1; // " ++ nameE xx) : initme rr
              initme (xx@(ER _):rr) = ("\tsofar += Nx*Ny*Nz; // " ++ nameE xx) : initme rr
              initme _ = error "bug inin setarg initme"
      sizeE (ES _) = "1"
      sizeE (ER _) = "Nx*Ny*Nz"
      sizeE (EK _) = error "no sizeE for EK yet"
      inputs :: Expression Scalar -> [Exprn]
      inputs x = findOrderedInputs x -- Set.toList $ findInputs x -- 
      codeMutableData a = unlines $ map (\x -> "\tmutable double " ++ x ++ ";") a

ctype :: Exprn -> C.Type
ctype (ES _) = Double
ctype _ = Vector

create0dHeader :: Expression Scalar -> String -> String
create0dHeader e0 n =
  unlines $
  ["// -*- mode: C++; -*-",
   "",
   "#include \"new/NewFunctional.h\"",
   "#include \"utilities.h\"",
   "#include \"handymath.h\"",
   "",
   "",
   "class " ++ n ++ " : public NewFunctional {",
   "public:",
   "\t" ++ n ++ "() { data = Vector(" ++ show (length $ findOrderedInputs e) ++ "); };",
   "",
   functionDeclaration "true_energy" "double" [],
   functionDeclaration "grad" "Vector" [],
   functionDeclaration "printme" "void" [("const char *", "prefix")]]++
  map setarg (findOrderedInputs e0) ++
 ["private:",
  ""++ codeMutableData (Set.toList $ findNamedScalars e)  ++"}; // End of " ++ n ++ " class"]
    where
      e = mapExpression' renameVar e0
      renameVar (Var CannotBeFreed a b c x) = Var CannotBeFreed ("var"++a) ("var"++b) c x
      renameVar x = x
      -- setarg creates a method that will get a reference to a given
      -- input argument's value.
      setarg :: Exprn -> String
      setarg ee =
        newcode $ CFunction {
          name = nameE ee,
          returnType = Reference (ctype ee),
          constness = "const",
          args = [],
          contents = "\tint sofar = 0;" : initme (findOrderedInputs e0)
          }
        where initme (xx@(ES _):rr)
                | xx == ee = ["\treturn data[sofar];"]
                | otherwise = ("\tsofar += 1; // " ++ nameE xx) : initme rr
              initme _ = error "bug inin setarg initme"
      codeMutableData a = unlines $ map (\x -> "\tmutable double " ++ x ++ ";") a


createCppFile :: Expression Scalar -> [Exprn] -> String -> String -> String
createCppFile e = if any is_nonscalar (findOrderedInputs e)
                  then create3dCppFile e
                  else create0dCppFile e
  where is_nonscalar (ES _) = False
        is_nonscalar _ = True

create0dCppFile :: Expression Scalar -> [Exprn] -> String -> String -> String
create0dCppFile e variables n headername =
  unlines $
  ["// -*- mode: C++; -*-",
   "",
   "#include \"" ++ headername ++ "\"",
   "",
   "",
   newcode $ CFunction {
     name = n++"::true_energy",
     returnType = Double,
     constness = "const",
     args = [],
     contents = ["\tint sofar = 0;"] ++
                map createInput (findOrderedInputs e) ++
                [newcodeStatements (fst energy),
                 "\treturn " ++ newcode (snd energy) ++ ";\n"]
     },
   newcode $ CFunction {
     name = n++"::grad",
     returnType = Vector,
     constness = "const",
     args = [],
     contents = evalv grade
     },
   newcode $ CFunction {
     name = n++"::printme",
     returnType = Void,
     constness = "const",
     args = [(ConstCharPtr, "prefix")],
     contents = map printEnergy $
                filter (`notElem` ["dV", "dr", "volume"]) $
                (Set.toList (findNamedScalars e))
     }] ++
 ["// End of " ++ n ++ " class",
  "// Total " ++ (show $ (countFFT (fst energy) + countFFT (fst grade))) ++ " Fourier transform used.",
  "// peak memory used: " ++ (show $ maximum $ map peakMem [fst energy, fst grade])
  ]
    where
      createInput ee@(ES _) = "\tdouble " ++ nameE ee ++ " = data[sofar]; sofar += 1;"
      createInput ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      createInputAndGrad ee@(ES _) = "\tdouble " ++ nameE ee ++ " = data[sofar];\n" ++
                                     "\tVector actual_grad_" ++ nameE ee ++ " = data.slice(sofar,1); " ++
                                     "sofar += 1;"
      createInputAndGrad ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      maxlen = 1 + maximum (map length $ "total energy" : Set.toList (findNamedScalars e))
      pad nn s | nn <= length s = s
      pad nn s = ' ' : pad (nn-1) s
      printEnergy v = "\tprintf(\"%s" ++ pad maxlen v ++ " =\", prefix);\n" ++
                      "\tprint_double(\"\", " ++ v ++ ");\n" ++
                      "\tprintf(\"\\n\");"
      energy = codex e
      the_actual_gradients = map (\(ES x) ->
                                   ES $ var ("grad_"++nameE (mkExprn x))
                                            ("grad_"++nameE (mkExprn x)) $derive x 1 e) variables
      grade :: ([Statement], [Exprn])
      grade = if variables == []
              then ([],[])
              else case optimize the_actual_gradients of
                (st0, es) -> let st = filter (not . isns) st0
                                 isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
                                 isns _ = False
                                 ns = findNamedScalars e
                                 _:revst = reverse $ reuseVar $ freeVectors $ st ++ map (\e' -> Assign (justvarname e') e') es
                             in (reverse revst, es)
      justvarname (ES (Var a b c d _)) = ES $ Var a (b++"[0]") c d Nothing
      justvarname _ = error "bad in justvarname"
      evalv :: ([Statement], [Exprn]) -> [String]
      evalv (st,ee) = ["\tVector output(data.get_size());",
                       "\toutput = 0;"]++
                      "\tint sofar = 0;" : map createInputAndGrad (findOrderedInputs e) ++
                      [newcodeStatements (st ++ concatMap assignit ee),
                       "\treturn output;"]
        where assignit eee = [Assign (ES $ Var CannotBeFreed
                                      ("actual_"++nameE eee++"[0]")
                                      ("actual_"++nameE eee++"[0]")
                                      ("actual_"++nameE eee) Nothing) eee]
      codex :: Expression Scalar -> ([Statement], Exprn)
      codex x = (init $ reuseVar $ freeVectors $ st ++ [Assign e' e'], e')
        where (st0, [e']) = optimize [ES $ factorize $ joinFFTs x]
              st = filter (not . isns) st0
              isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
              isns _ = False
              ns = findNamedScalars e


create3dCppFile :: Expression Scalar -> [Exprn] -> String -> String -> String
create3dCppFile e variables n headername =
  unlines $
  ["// -*- mode: C++; -*-",
   "",
   "#include \"" ++ headername ++ "\"",
   "",
   "",
   newcode $ CFunction {
     name = n++"::true_energy",
     returnType = Double,
     constness = "const",
     args = [],
     contents = ["\tint sofar = 0;"] ++
                map createInput (inputs e) ++
                [newcodeStatements (fst energy),
                 "\treturn " ++ newcode (snd energy) ++ ";\n"]
     },
   newcode $ CFunction {
     name = n++"::grad",
     returnType = Vector,
     constness = "const",
     args = [],
     contents = evalv grade
     },
   newcode $ CFunction {
     name = n++"::printme",
     returnType = Void,
     constness = "const",
     args = [(ConstCharPtr, "prefix")],
     contents = map printEnergy $
                filter (`notElem` ["dV", "dr", "volume"]) $
                (Set.toList (findNamedScalars e))
     },
   newcode $ CFunction {
     name = n++"::"++n,
     returnType = None,
     constness = "",
     args = [(Int, "myNx"), (Int, "myNy"), (Int, "myNz")],
     contents =["\tdata = Vector(int(" ++ code (sum $ map actualsize $ inputs e) ++ "));",
                "\tNx() = myNx;",
                "\tNy() = myNy;",
                "\tNz() = myNz; // good"]
     }] ++
 ["// End of " ++ n ++ " class",
  "// Total " ++ (show $ (countFFT (fst energy) + countFFT (fst grade))) ++ " Fourier transform used.",
  "// peak memory used: " ++ (show $ maximum $ map peakMem [fst energy, fst grade])
  ]
    where
      actualsize (ES _) = 1 :: Expression Scalar
      actualsize (ER _) = s_var "myNx" * s_var "myNy" * s_var "myNz"
      actualsize (EK _) = error "need to compute size of EK in actualsize of NewCode"
      createInput ee@(ES _) = "\tdouble " ++ nameE ee ++ " = data[sofar]; sofar += 1;"
      createInput ee@(ER _) = "\tVector " ++ nameE ee ++ " = data.slice(sofar,Nx*Ny*Nz); sofar += Nx*Ny*Nz;"
      createInput ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      createInputAndGrad ee@(ES _) = "\tdouble " ++ nameE ee ++ " = data[sofar];\n" ++
                                     "\tVector grad_" ++ nameE ee ++ " = data.slice(sofar,1); " ++
                                     "sofar += 1;"
      createInputAndGrad ee@(ER _) = "\tVector " ++ nameE ee ++ " = data.slice(sofar,Nx*Ny*Nz);\n"++
                                     "\tVector grad_" ++ nameE ee ++ " = output.slice(sofar,Nx*Ny*Nz); " ++
                                     "sofar += Nx*Ny*Nz;"
      createInputAndGrad ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      inputs :: Expression Scalar -> [Exprn]
      inputs x = findOrderedInputs x -- Set.toList $ findInputs x -- 
      maxlen = 1 + maximum (map length $ "total energy" : Set.toList (findNamedScalars e))
      pad nn s | nn <= length s = s
      pad nn s = ' ' : pad (nn-1) s
      printEnergy v = "\tprintf(\"%s" ++ pad maxlen v ++ " =\", prefix);\n" ++
                      "\tprint_double(\"\", " ++ v ++ ");\n" ++
                      "\tprintf(\"\\n\");"
      energy = codex e
      the_actual_gradients = map (mapExprn (\x -> mkExprn $ var ("grad_"++nameE (mkExprn x))
                                                                ("grad_"++nameE (mkExprn x)) $
                                                  derive x 1 e)) variables
      grade :: ([Statement], [Exprn])
      grade = if variables == []
              then ([],[])
              else case optimize the_actual_gradients of
                (st0, es) -> let st = filter (not . isns) st0
                                 isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
                                 isns _ = False
                                 ns = findNamedScalars e
                                 _:revst = reverse $ reuseVar $ freeVectors $ st ++ map (\e' -> Assign (justvarname e') e') es
                             in (reverse revst, es)
      justvarname (ES (Var a b c d _)) = ES $ Var a (b++"[0]") c d Nothing
      justvarname (ER (Var a b c d _)) = ER $ Var a b c d Nothing
      justvarname (EK (Var a b c d _)) = EK $ Var a b c d Nothing
      justvarname _ = error "bad in justvarname"
      evalv :: ([Statement], [Exprn]) -> [String]
      evalv (st,ee) = ["\tVector output(data.get_size());",
                       "\tfor (int i=0;i<data.get_size();i++) {",
                       "\t\toutput[i] = 0;",
                       "\t}"]++
                      "\tint sofar = 0;" : map createInputAndGrad (inputs e) ++
                      [newcodeStatements (st ++ concatMap assignit ee),
                       "\treturn output;"]
        where assignit eee = [Assign (justvarname eee) eee]
      codex :: Expression Scalar -> ([Statement], Exprn)
      codex x = (init $ reuseVar $ freeVectors $ st ++ [Assign e' e'], e')
        where (st0, [e']) = optimize [ES $ factorize $ joinFFTs x]
              st = filter (not . isns) st0
              isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
              isns _ = False
              ns = findNamedScalars e

