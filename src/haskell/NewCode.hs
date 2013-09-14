{-# LANGUAGE PatternGuards #-}

module NewCode (module Statement,
                module Expression,
                defineFunctional, createCppFile, createHeader )
    where

import Statement
import Expression
import Optimize ( optimize )
import qualified Data.Set as Set


defineFunctional :: Expression Scalar -> [Exprn] -> [Exprn] -> String -> String
defineFunctional e arg variables n = createHeader e arg n ++ implementation
  where implementation = unlines $ drop 4 $ lines $ createCppFile e variables n ""


functionCode :: String -> String -> [(String, String)] -> String -> String
functionCode "" "" [] "" = ""
functionCode "" "" (x:xs) "" = if xs == []
                               then fst x ++ " " ++ snd x
                               else fst x ++ " " ++ snd x ++ ", " ++ functionCode "" "" xs ""
functionCode n t a b = t ++ " " ++ n ++ "(" ++ functionCode "" "" a "" ++ ") const {\n" ++ b ++ "}\n"

functionDeclaration :: String -> String -> [(String, String)] -> String
functionDeclaration "" "" [] = ""
functionDeclaration "" "" (x:xs) =
  if xs == []
  then fst x ++ " " ++ snd x
  else fst x ++ " " ++ snd x ++ ", " ++ functionDeclaration "" "" xs
functionDeclaration n t a =
  "\t" ++ t ++ " " ++ n ++ "(" ++ functionDeclaration "" "" a ++ ") const;"


createHeader :: Expression Scalar -> [Exprn] -> String -> String
createHeader e0 arg0 n =
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
   "" ++ n ++ codeA arg ++ "  {",
   "}",
   "",
   functionDeclaration "energy" "double" [("const Vector", "&xxx")],
   functionDeclaration "energy_per_volume" "double" [("const Vector", "&xxx")],
   functionDeclaration "denergy_per_volume_dx" "double" [("const Vector", "&xxx")],
   functionDeclaration "grad" "Vector" [("const Vector", "&xxx")],
   functionDeclaration "printme" "void" [("const char *", "prefix")],
   functionDeclaration "createInput" "Vector" (codeInputArgs (inputs e0))] ++
  map setarg (inputs e0) ++
 ["private:",
  ""++ codeArgInit arg ++ codeMutableData (Set.toList $ findNamedScalars e)  ++"}; // End of " ++ n ++ " class"]
    where
      e = mapExpression' renameVar e0
      arg = map renameExprn arg0
      renameVar (Var CannotBeFreed a b c x) = Var CannotBeFreed ("var"++a) ("var"++b) c x
      renameVar x = x
      renameExprn (ES x) = ES $ mapExpression' renameVar x
      renameExprn (EK x) = EK $ mapExpression' renameVar x
      renameExprn (ER x) = ER $ mapExpression' renameVar x
      -- setarg creates a method that will get a reference to a given
      -- input argument's value.
      setarg :: Exprn -> String
      setarg ee = functionCode (nameE ee) (newreferenceE ee) [("Vector", "xxx")] $
                  unlines ("\tint sofar = 0;" : initme (inputs e0) )
                    where initme (xx@(ES _):_) | xx == ee = ["\treturn xxx[sofar];"]
                          initme (xx:_) | xx == ee = ["\treturn xxx.slice(sofar," ++ sizeE ee ++ ");"]
                          initme (xx@(ES _):rr)
                            | nameE xx `elem` ["Nx","Ny","Nz"] =
                              ("\tconst double "++nameE xx++" = xxx[sofar++];") : initme rr
                          initme (xx@(ES _):rr) = ("\tsofar += 1; // " ++ nameE xx) : initme rr
                          initme (xx@(ER _):rr) = ("\tsofar += Nx*Ny*Nz; // " ++ nameE xx) : initme rr
                          initme _ = error "bug inin setarg initme"
      sizeE (ES _) = "1"
      sizeE (ER _) = "Nx*Ny*Nz"
      sizeE (EK _) = error "no sizeE for EK yet"
      inputs :: Expression Scalar -> [Exprn]
      inputs x = findOrderedInputs x -- Set.toList $ findInputs x -- 
      codeA :: [Exprn] -> String
      codeA [] = "()"
      codeA a = "(" ++ foldl1 (\x y -> x ++ ", " ++ y ) (map (\x -> "double " ++ nameE x ++ "_arg") a) ++ ") : " ++ foldl1 (\x y -> x ++ ", " ++ y) (map (\x -> nameE x ++ "(" ++ nameE x ++ "_arg)") a)
      codeInputArgs = map (\x -> (newdeclareE x, nameE x))
      codeArgInit a = unlines $ map (\x -> "\tdouble " ++ nameE x ++ ";") a
      codeMutableData a = unlines $ map (\x -> "\tmutable double " ++ x ++ ";") a


createCppFile :: Expression Scalar -> [Exprn] -> String -> String -> String
createCppFile e variables n headername =
  unlines $
  ["// -*- mode: C++; -*-",
   "",
   "#include \"" ++ headername ++ "\"",
   "",
   "",
   functionCode (n++"::energy") "double" [("const Vector", "&xxx")]
      (unlines $
       ["\tint sofar = 0;"] ++
       map createInput (inputs e) ++
       [newcodeStatements (fst energy),
       "\treturn " ++ newcode (snd energy) ++ ";\n"]),
   functionCode (n++"::energy_per_volume") "double" [("const Vector", "&xxx")]
      (unlines $
       ["\tint sofar = 0;"] ++
       map createInput (inputs $ makeHomogeneous e) ++
       [newcodeStatements (fst energy_per_volume),
       "\treturn " ++ newcode (snd energy_per_volume) ++ ";\n"]),
   functionCode (n++"::denergy_per_volume_dx") "double" [("const Vector", "&xxx")]
    (newcodeStatements (fst denergy_per_volume_dx) ++
     "\treturn " ++ newcode (snd denergy_per_volume_dx) ++ ";\n"),
   functionCode (n++"::grad") "Vector" [("const Vector", "&xxx")] (evalv grade),
   functionCode (n++"::printme") "void" [("const char *", "prefix")]
      (unlines $ map printEnergy $
       filter (`notElem` ["dV", "dr", "volume"]) $
       (Set.toList (findNamedScalars e))),
   functionCode (n++"::createInput") "Vector" (codeInputArgs (inputs e))
      (unlines $ ["\tconst int newsize = " ++
                  xxx (map getsize $ inputs e),
                  "\tVector out(newsize);",
                  "\tint sofar = 0;"] ++
                 map initializeOut (inputs e) ++
                 ["\treturn out;"])] ++
 ["// End of " ++ n ++ " class",
  "// Total " ++ (show $ (countFFT (fst energy) + countFFT (fst grade))) ++ " Fourier transform used.",
  "// peak memory used: " ++ (show $ maximum $ map peakMem [fst energy, fst grade])
  ]
    where
      initializeOut x@(ES _) = concat ["\tout[sofar] = ",nameE x,";\n\tsofar += 1;"]
      initializeOut x = concat ["\tout.slice(sofar,", getsize x, ") = ",nameE x,";\n\t",
                                "sofar += ",getsize x,";"]
      getsize (ES _) = "1"
      getsize ee = nameE ee ++ ".get_size()"
      createInput ee@(ES _) = "\tdouble " ++ nameE ee ++ " = xxx[sofar]; sofar += 1;"
      createInput ee@(ER _) = "\tVector " ++ nameE ee ++ " = xxx.slice(sofar,Nx*Ny*Nz); sofar += Nx*Ny*Nz;"
      createInput ee = error ("unhandled type in NewCode scalarClass: " ++ show ee)
      createInputAndGrad ee@(ES _) = "\tdouble " ++ nameE ee ++ " = xxx[sofar];\n" ++
                                     "\tVector grad_" ++ nameE ee ++ " = xxx.slice(sofar,1); " ++
                                     "sofar += 1;"
      createInputAndGrad ee@(ER _) = "\tVector " ++ nameE ee ++ " = xxx.slice(sofar,Nx*Ny*Nz);\n"++
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
      energy_per_volume = codex (makeHomogeneous e)
      denergy_per_volume_dx = codex (derive (s_var "xxx" :: Expression Scalar) 1 $ makeHomogeneous e)
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
      evalv :: ([Statement], [Exprn]) -> String
      evalv (st,ee) = unlines (["\tVector output(xxx.get_size());",
                                "\tfor (int i=0;i<xxx.get_size();i++) {",
                                "\t\toutput[i] = 0;",
                                "\t}"]++
                               "\tint sofar = 0;" : map createInputAndGrad (inputs e) ++
                               [newcodeStatements (st ++ concatMap assignit ee),
                                "\treturn output;"])
        where assignit eee = [Assign (justvarname eee) eee]
      codex :: Expression Scalar -> ([Statement], Exprn)
      codex x = (init $ reuseVar $ freeVectors $ st ++ [Assign e' e'], e')
        where (st0, [e']) = optimize [ES $ factorize $ joinFFTs x]
              st = filter (not . isns) st0
              isns (Initialize (ES (Var _ _ s _ Nothing))) = Set.member s ns
              isns _ = False
              ns = findNamedScalars e
      codeInputArgs = map (\x -> (newdeclareE x, nameE x))
      xxx [] = ""
      xxx iii = foldl1 (\x y -> x ++ " + " ++ y) iii ++";"

