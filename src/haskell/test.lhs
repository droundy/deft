\begin{code}
import CodeGen
import Test.HUnit
{-
nbarTests :: Test
nbarTests = TestList [nbar_test homolognbar]
  where nbar_test e = TestCase $ assertEqual (setZero Kx $ setZero Ky $ setZero Kz e) (setZero Kz $ setZero Ky $ setZero Kz e)
        spreading = 6.0
        kdr = k * s_var "dr"
        kR = k * s_var "R"
        nbar = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft (r_var "x"))
        dlognbar = derive (R "x") 1 (log nbar)
        homolognbar = setZero (S "dr") dlognbar -- The problem is setZero Kz.
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
          --writeFile "tests/generated-haskell/nice-logoneminusnbar.h" $ generateHeader (log (1-ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft (r_var "x")))) (Just (r_var "R")) "NiceLogOneMinusNbar"
          c <- runTestTT $ TestList []
          if failures c > 0
            then fail $ "Failed " ++ show (failures c) ++ " tests."
            else putStrLn "All tests passed!"
\end{code}
