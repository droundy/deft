module C ( Type(..), CFunction(..) ) where

import Expression ( Code(..) )

-- Type None means it is a constructor or destructor
data Type = Void | Double | Int | Vector | ConstCharPtr | None | Reference Type

instance Code Type where
  codePrec _ Void = showString "void "
  codePrec _ Double = showString "double "
  codePrec _ Int = showString "int "
  codePrec _ Vector = showString "Vector "
  codePrec _ None = showString ""
  codePrec _ _ = error "codePrec Type incomplete"
  newcodePrec _ Void = showString "void "
  newcodePrec _ Double = showString "double "
  newcodePrec _ Int = showString "int "
  newcodePrec _ Vector = showString "Vector "
  newcodePrec _ None = showString ""
  newcodePrec _ ConstCharPtr = showString "const char *"
  newcodePrec _ (Reference Double) = showString "double &"
  newcodePrec _ (Reference Vector) = showString "Vector "
  newcodePrec _ (Reference x) = error ("Cannot create reference to type " ++ newcode x)
  latexPrec _ Void = showString "\\textrm{void}"
  latexPrec _ Double = showString "\\textrm{double}"
  latexPrec _ Int = showString "\\textrm{int}"
  latexPrec _ Vector = showString "\\textrm{Vector}"
  latexPrec _ None = showString ""
  latexPrec _ _ = error "latexPrec Type incomplete"

data CFunction = CFunction {
  returnType :: Type,
  name :: String,
  constness :: String,
  args :: [(Type, String)],
  contents :: [String] }

instance Code CFunction where
  codePrec _ _ = error "codePrec not implemented for CFunction"
  newcodePrec _ f =
    showString $ unlines $
    [newcode (returnType f) ++ name f ++ "(" ++ showargs (args f) ++ ") " ++ constness f ++ " {"] ++
    map ("\t"++) (contents f) ++
    ["}"]
  latexPrec _ _ = error "latexPrec not implemented for CFunction"

showargs :: [(Type, String)] -> String
showargs [] = ""
showargs [(t,a)] = newcode t ++ a
showargs ((t,a):r) = newcode t ++ a ++ ", " ++ showargs r
