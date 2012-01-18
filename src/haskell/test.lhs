\begin{code}
import CodeGen
import Test.HUnit

{-gradTests :: Test
gradTests = TestList [test_grad "foo" (integrate $ foo*foo) (dV*2*foo),
                      test_grad "foo" (integrate $ foo**2) (dV*2*foo),
                      test_grad "x" (integrate $ x * ifft (5*fft s)) (dV*5*s),
                      test_grad "x" (integrate $ ifft (2 * fft x)) (dV*2)]
  where test_grad v e g = TestCase $ assertEqual (show g) g (grad v e)
        foo = r_var "foo"
        dV = s_var "dV"
        x = r_var "x"
        s = s_var "s"

showTests :: Test
showTests = TestList [show_test "x[i]" x,
                      show_test "x[i]*x[i]" (x**2),
                      show_test "R \"x\"*R \"x\"*(R \"x\"*R \"x\")" (x**4)]
  where show_test str e = TestCase $ assertEqual str str (show e)
        x = r_var "x"

codeTests :: Test
codeTests = TestList [show_test "x[i]" x,
                      show_test "x[i]*x[i]" (x**2),
                      show_test "x[i]*x[i]*(x[i]*x[i])" (x**4)]
  where show_test str e = TestCase $ assertEqual str str (code e)
        x = r_var "x"

showStatements :: Test
showStatements = TestList [ss "for (int i=0; i<gd.NxNyNz; i++) {    \nx = 5.0;\n}\n" 
                               ("x" := (5 :: Expression RealSpace)),
                           ss "for (int i=0; i<gd.NxNyNz; i++) {    \nx[i] = y[i];\n}\n" 
                               ("x" := y),
                           ss "for (int i=0; i<gd.NxNyNz; i++) {    \nx[i] = ifft(exp((-1.0)*(khere[0]*khere[0] + khere[1]*khere[1] + khere[2]*khere[2]))*fft(y[i]));\n}\n" 
                               ("x" := ifft (exp (-ksqr) * fft y)),
                           ss "for (int i=0; i<gd.NxNyNz; i++) {    \nx[i] = ifft(exp((-1.0)*(khere[0]*khere[0] + khere[1]*khere[1] + khere[2]*khere[2]))*fft(y[i]));\n}\n" 
                               ("x" := ifft (exp (-k**2) * fft y)),
                           ss "for (int i=0; i<gd.NxNyNz; i++) {    \nx[i] = y[i];\n}\n\nGrid a(gd);\nfor (int i=0; i<gd.NxNyNz; i++) {    \na[i] := y[i];\n}\n" $ 
                           do "x" := y
                              "a" :?= y]
  where ss str st = TestCase $ assertEqual str str (show st)
        y = r_var "y"
-}
main :: IO ()
main = do writeFile "tests/generated-haskell/nice-sum.h" $ generateHeader (r_var "x" + s_var "kT") Nothing "NiceSum"
          writeFile "tests/generated-haskell/nice-quad.h" $ generateHeader ((r_var "x" + s_var"kT")**2 - r_var "x" + 2*s_var "kT") Nothing "NiceQuad"
          writeFile "tests/generated-haskell/nice-sqrt.h" $ generateHeader (r_var "x"**0.5) Nothing "NiceSqrt"
          writeFile "tests/generated-haskell/nice-sqrtandmore.h" $ generateHeader ((r_var "x"**0.5) - r_var "x" + 2*s_var "kT") Nothing "NiceSqrtandMore"
          writeFile "tests/generated-haskell/nice-log.h" $ generateHeader (log (r_var "x")) Nothing "NiceLog"
          writeFile "tests/generated-haskell/nice-logandsqr.h" $ generateHeader (log (r_var "x") + (r_var "x")**2) Nothing "NiceLogandSqr"
          writeFile "tests/generated-haskell/nice-logandsqrandinverse.h" $ generateHeader (log (r_var "x") + (r_var "x")**2 - (r_var "x")**3 + (r_var "x")**(-1)) Nothing "NiceLogandSqrandInverse"
          writeFile "tests/generated-haskell/nice-logoneminusx.h" $ generateHeader (log (1 - r_var "x")) Nothing "NiceLogOneMinusX"
          let spreading = 6.0
              kdr = k * s_var "dr"
              kR = k * s_var "R"
              nbar = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft (r_var "x"))
          writeFile "tests/generated-haskell/nice-nbar.h" $ generateHeader nbar (Just (r_var "R")) "NiceNbar"
          c <- runTestTT $ TestList []
          if failures c > 0
            then fail $ "Failed " ++ show (failures c) ++ " tests."
            else putStrLn "All tests passed!"
\end{code}
