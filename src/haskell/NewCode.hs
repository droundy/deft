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
   "public:"] ++ map (declare . strip_type) (create3dMethods e0 [] n) ++
  map (newcode . setarg) (inputs e0) ++
 ["private:",
  ""++ codeMutableData (Set.toList $ findNamedScalars e)  ++"}; // End of " ++ n ++ " class"]
    where
      strip_type f = f { name = drop (length n+2) $ name f }
      e = mapExpression' renameVar e0
      renameVar (Var CannotBeFreed a b c x) = Var CannotBeFreed ("var"++a) ("var"++b) c x
      renameVar x = x
      -- setarg creates a method that will get a reference to a given
      -- input argument's value.
      setarg :: Exprn -> CFunction
      setarg ee = CFunction {
        name = nameE ee,
        returnType = Reference (ctype ee),
        constness = "const",
        args = [],
        contents = "int sofar = 0;" : initme (inputs e0)
        }
        where initme (xx@(ES _):_) | xx == ee = ["return data[sofar];"]
              initme (xx:_) | xx == ee = ["return data.slice(sofar," ++ sizeE ee ++ ");"]
              initme (xx@(ES _):rr)
                | nameE xx `elem` ["Nx","Ny","Nz"] =
                  ("const double "++nameE xx++" = data[sofar++];") : initme rr
              initme (xx@(ES _):rr) = ("sofar += 1; // " ++ nameE xx) : initme rr
              initme (xx@(ER _):rr) = ("sofar += Nx*Ny*Nz; // " ++ nameE xx) : initme rr
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
   "public:"] ++ map (declare . strip_type) (create0dMethods e0 [] n) ++
  map setarg (findOrderedInputs e0) ++
 ["private:",
  ""++ codeMutableData (Set.toList $ findNamedScalars e)  ++"}; // End of " ++ n ++ " class"]
    where
      strip_type f = f { name = drop (length n+2) $ name f }
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
          contents = "int sofar = 0;" : initme (findOrderedInputs e0)
          }
        where initme (xx@(ES _):rr)
                | xx == ee = ["return data[sofar];"]
                | otherwise = ("sofar += 1; // " ++ nameE xx) : initme rr
              initme _ = error "bug inin setarg initme"
      codeMutableData a = unlines $ map (\x -> "\tmutable double " ++ x ++ ";") a


createCppFile :: Expression Scalar -> [Exprn] -> String -> String -> String
createCppFile e variables n headername = unlines $ ["// -*- mode: C++; -*-",
                                                    "",
                                                    "#include \"" ++ headername ++ "\"",
                                                    "",
                                                    ""] ++ map newcode methods
  where is_nonscalar (ES _) = False
        is_nonscalar _ = True
        methods = if any is_nonscalar (findOrderedInputs e)
                  then create3dMethods e variables n
                  else create0dMethods e variables n

create0dMethods :: Expression Scalar -> [Exprn] -> String -> [CFunction]
create0dMethods e variables n =
  [CFunction {
      name = n++"::true_energy",
      returnType = Double,
      constness = "const",
      args = [],
      contents = ["int sofar = 0;"] ++
                 map createInput (findOrderedInputs e) ++
                 [newcodeStatements (fst energy),
                  "return " ++ newcode (snd energy) ++ ";\n"]
      },
   CFunction {
     name = n++"::grad",
     returnType = Vector,
     constness = "const",
     args = [],
     contents = evalv grade
     },
   CFunction {
     name = n++"::printme",
     returnType = Void,
     constness = "const",
     args = [(ConstCharPtr, "prefix")],
     contents = map printEnergy $
                filter (`notElem` ["dV", "dr", "volume"]) $
                (Set.toList (findNamedScalars e))
     }]
    where
      createInput ee@(ES _) = "double " ++ nameE ee ++ " = data[sofar]; sofar += 1;"
      createInput ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      createInputAndGrad ee@(ES _) = "double " ++ nameE ee ++ " = data[sofar];\n" ++
                                     "Vector actual_grad_" ++ nameE ee ++ " = data.slice(sofar,1); " ++
                                     "sofar += 1;"
      createInputAndGrad ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      maxlen = 1 + maximum (map length $ "total energy" : Set.toList (findNamedScalars e))
      pad nn s | nn <= length s = s
      pad nn s = ' ' : pad (nn-1) s
      printEnergy v = "printf(\"%s" ++ pad maxlen v ++ " =\", prefix);\n" ++
                      "print_double(\"\", " ++ v ++ ");\n" ++
                      "printf(\"\\n\");"
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
      evalv (st,ee) = ["Vector output(data.get_size());",
                       "output = 0;"]++
                      "int sofar = 0;" : map createInputAndGrad (findOrderedInputs e) ++
                      [newcodeStatements (st ++ concatMap assignit ee),
                       "return output;"]
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

create3dMethods :: Expression Scalar -> [Exprn] -> String -> [CFunction]
create3dMethods e variables n =
  [CFunction {
     name = n++"::"++n,
     returnType = None,
     constness = "",
     args = [(Int, "myNx"), (Int, "myNy"), (Int, "myNz")],
     contents =["data = Vector(int(" ++ code (sum $ map actualsize $ inputs e) ++ "));",
                "Nx() = myNx;",
                "Ny() = myNy;",
                "Nz() = myNz;"]
     },
   CFunction {
     name = n++"::"++n,
     returnType = None,
     constness = "",
     args = [(Double, "ax"), (Double, "ay"), (Double, "az"), (Double, "dx")],
     contents =["int myNx = int(ceil(ax/dx));",
                "int myNy = int(ceil(ay/dx));",
                "int myNz = int(ceil(az/dx));",
                "data = Vector(int(" ++ code (sum $ map actualsize $ inputs e) ++ "));",
                "Nx() = myNx;",
                "Ny() = myNy;",
                "Nz() = myNz;",
                "a1() = ax;",
                "a2() = ay;",
                "a3() = az;"]
     },
   CFunction {
      name = n++"::true_energy",
      returnType = Double,
      constness = "const",
      args = [],
      contents = ["int sofar = 0;"] ++
                 map createInput (inputs e) ++
                 [newcodeStatements (fst energy),
                  "return " ++ newcode (snd energy) ++ ";\n"]
      },
   CFunction {
     name = n++"::grad",
     returnType = Vector,
     constness = "const",
     args = [],
     contents = evalv grade
     },
   CFunction {
     name = n++"::printme",
     returnType = Void,
     constness = "const",
     args = [(ConstCharPtr, "prefix")],
     contents = map printEnergy $
                filter (`notElem` ["dV", "dr", "volume"]) $
                (Set.toList (findNamedScalars e))
     }]
    where
      actualsize (ES _) = 1 :: Expression Scalar
      actualsize (ER _) = s_var "myNx" * s_var "myNy" * s_var "myNz"
      actualsize (EK _) = error "need to compute size of EK in actualsize of NewCode"
      createInput ee@(ES _) = "double " ++ nameE ee ++ " = data[sofar]; sofar += 1;"
      createInput ee@(ER _) = "Vector " ++ nameE ee ++ " = data.slice(sofar,Nx*Ny*Nz); sofar += Nx*Ny*Nz;"
      createInput ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      createInputAndGrad ee@(ES _) = "double " ++ nameE ee ++ " = data[sofar];\n" ++
                                     "Vector grad_" ++ nameE ee ++ " = data.slice(sofar,1); " ++
                                     "sofar += 1;"
      createInputAndGrad ee@(ER _) = "Vector " ++ nameE ee ++ " = data.slice(sofar,Nx*Ny*Nz);\n"++
                                     "Vector grad_" ++ nameE ee ++ " = output.slice(sofar,Nx*Ny*Nz); " ++
                                     "sofar += Nx*Ny*Nz;"
      createInputAndGrad ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      inputs :: Expression Scalar -> [Exprn]
      inputs x = findOrderedInputs x -- Set.toList $ findInputs x -- 
      maxlen = 1 + maximum (map length $ "total energy" : Set.toList (findNamedScalars e))
      pad nn s | nn <= length s = s
      pad nn s = ' ' : pad (nn-1) s
      printEnergy v = "printf(\"%s" ++ pad maxlen v ++ " =\", prefix);\n" ++
                      "print_double(\"\", " ++ v ++ ");\n" ++
                      "printf(\"\\n\");"
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
      evalv (st,ee) = ["Vector output(data.get_size());",
                       "for (int i=0;i<data.get_size();i++) {",
                       "\toutput[i] = 0;",
                       "}"]++
                      "int sofar = 0;" : map createInputAndGrad (inputs e) ++
                      [newcodeStatements (st ++ concatMap assignit ee),
                       "return output;"]
        where assignit eee = [Assign (justvarname eee) eee]
      codex :: Expression Scalar -> ([Statement], Exprn)
      codex x = (init $ reuseVar $ freeVectors $ st ++ [Assign e' e'], e')
        where (st0, [e']) = optimize [ES $ factorize $ joinFFTs x]
              st = filter (not . isns) st0
              isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
              isns _ = False
              ns = findNamedScalars e

