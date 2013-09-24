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
      sizeE (EK _) = "Nx*Ny*(int(Nz)/2+1)"
      inputs :: Expression Scalar -> [Exprn]
      inputs x = findOrderedInputs x -- Set.toList $ findInputs x -- 
      codeMutableData a = unlines $ map (\x -> "\tmutable double " ++ x ++ ";") a

ctype :: Exprn -> C.Type
ctype (ES _) = Double
ctype (ER _) = Vector
ctype (EK _) = ComplexVector

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
create0dMethods e variables n = createAnydMethods e variables n

create3dMethods :: Expression Scalar -> [Exprn] -> String -> [CFunction]
create3dMethods e variables n =
  [CFunction {
     name = n++"::"++n,
     returnType = None,
     constness = "",
     args = [(Int, "myNx"), (Int, "myNy"), (Int, "myNz")],
     contents =["data = Vector(int(" ++ code (sum $ map actualsize $ findOrderedInputs e) ++ "));",
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
                "data = Vector(int(" ++ code (sum $ map actualsize $ findOrderedInputs e) ++ "));",
                "Nx() = myNx;",
                "Ny() = myNy;",
                "Nz() = myNz;",
                "a1() = ax;",
                "a2() = ay;",
                "a3() = az;"]
     }] ++ createAnydMethods e variables n
    where
      actualsize (ES _) = 1 :: Expression Scalar
      actualsize (ER _) = s_var "myNx" * s_var "myNy" * s_var "myNz"
      actualsize (EK _) = error "need to compute size of EK in actualsize of NewCode"

createAnydMethods :: Expression Scalar -> [Exprn] -> String -> [CFunction]
createAnydMethods e variables n =
  [CFunction {
      name = n++"::true_energy",
      returnType = Double,
      constness = "const",
      args = [],
      contents = ["int sofar = 0;"] ++
                 map createInput (findOrderedInputs e) ++
                 [newcodeStatements (eval_scalar e)]
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
      contents = concatMap printEnergy $
                 filter (`notElem` ["dV", "dr", "volume"]) $
                 (Set.toList (findNamedScalars e))
      }] ++ map getIntermediate (Set.toList $ findNamed e)
  where
      maxlen = 1 + maximum (map length $ "total energy" : Set.toList (findNamedScalars e))
      pad nn s | nn <= length s = s
      pad nn s = ' ' : pad (nn-1) s
      printEnergy v = ["printf(\"%s" ++ pad maxlen v ++ " =\", prefix);",
                       "print_double(\"\", " ++ v ++ ");",
                       "printf(\"\\n\");"]
      createInput ee@(ES _) = "double " ++ nameE ee ++ " = data[sofar]; sofar += 1;"
      createInput ee@(ER _) = "Vector " ++ nameE ee ++ " = data.slice(sofar,Nx*Ny*Nz); sofar += Nx*Ny*Nz;"
      createInput ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      getIntermediate :: (String, Exprn) -> CFunction
      getIntermediate ("", _) = error "empty string in getIntermediate"
      getIntermediate (rsname, a) = CFunction {
      name = n++"::get_"++rsname,
      returnType = ctype a,
      constness = "const",
      args = [],
      contents = ["int sofar = 0;"] ++
                 map createInput (findOrderedInputs e) ++
                 [newcodeStatements $ eval_named (""++rsname) a]
      }
      createInputAndGrad ee@(ES _) = ["double " ++ nameE ee ++ " = data[sofar];",
                                      "double *grad_" ++ nameE ee ++ " = &data[sofar++];"]
      createInputAndGrad ee@(ER _) = ["Vector " ++ nameE ee ++ " = data.slice(sofar,Nx*Ny*Nz);",
                                      "Vector grad_" ++ nameE ee ++ " = output.slice(sofar,Nx*Ny*Nz);",
                                      "sofar += Nx*Ny*Nz;"]
      createInputAndGrad ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      the_gradients = map (mapExprn (\x -> ("grad_"++nameE (mkExprn x), mkExprn $ derive x 1 e))) variables
      grade :: [Statement]
      grade = if variables == []
              then []
              else case optimize $ map snd the_gradients of
                (st0, es) -> let st = filter (not . isns) st0
                                 isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
                                 isns _ = False
                                 ns = findNamedScalars e
                             in reuseVar $ freeVectors $ st ++ zipWith Assign (map varn the_gradients) es
      varn (vnam, ES _) = ES $ Var CannotBeFreed ('*':vnam) ('*':vnam) vnam Nothing
      varn (vnam, ER _) = ER $ Var CannotBeFreed vnam ("grad_"++vnam++"[i]") vnam Nothing
      varn (vnam, EK _) = EK $ Var CannotBeFreed vnam (vnam++"[i]") vnam Nothing
      evalv :: [Statement] -> [String]
      evalv st = ["Vector output(data.get_size());",
                  "output = 0;"]++
                  "int sofar = 0;" : concatMap createInputAndGrad (findOrderedInputs e) ++
                  [newcodeStatements st,
                   "return output;"]

eval_scalar :: Expression Scalar -> [Statement]
eval_scalar x = reuseVar $ freeVectors $ st ++ [Return e']
  where (st0, [e']) = optimize [ES $ factorize $ joinFFTs x]
        st = filter (not . isns) st0
        isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
        isns _ = False
        ns = findNamedScalars x

eval_named :: String -> Exprn -> [Statement]
eval_named outname (ER x) = reuseVar $ freeVectors $ st ++ [Initialize outvar, Assign outvar e', Return outvar]
  where (st0, [e']) = optimize [ER $ factorize $ joinFFTs x]
        st = filter (not . isns) st0
        isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
        isns _ = False
        ns = findNamedScalars x
        outvar = ER $ Var IsTemp (outname++"[i]") outname outname Nothing
eval_named outname (EK x) = reuseVar $ freeVectors $ st ++ [Initialize outvar, Assign outvar e', Return outvar]
  where (st0, [e']) = optimize [EK $ factorize $ joinFFTs x]
        st = filter (not . isns) st0
        isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
        isns _ = False
        ns = findNamedScalars x
        outvar = EK $ Var IsTemp (outname++"[i]") outname outname Nothing
eval_named _ (ES x) = reuseVar $ freeVectors $ st ++ [Return e']
  where (st0, [e']) = optimize [ES $ factorize $ joinFFTs x]
        st = filter (not . isns) st0
        isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
        isns _ = False
        ns = findNamedScalars x
