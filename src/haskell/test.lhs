\begin{code}
import System.Directory ( createDirectoryIfMissing )
import CodeGen
import Test.HUnit

latexTests :: Test
latexTests = TestList [t "x" x,
                       t "0" (0 :: Expression RealSpace),
                       t "\\sqrt{x}" (sqrt x),
                       t "x^2" (x**2),
                       t "\\frac{1}{x^2}" (1/x**2),
                       t "\\frac{y^2}{x^2}" (y**2/x**2),
                       t "x^4" (x**4)]
  where t str e = TestCase $ assertEqual str str (latex e)
        x = r_var "x"
        y = s_var "y" :: Expression RealSpace

codeTests :: Test
codeTests = TestList [t "x[i]" x,
                      t "0" (0 :: Expression RealSpace),
                      t "sqrt(x[i])" (sqrt x),
                      t "x[i]*x[i]" (x**2),
                      t "1/(x[i]*x[i])" (1/x**2),
                      t "y*y/(x[i]*x[i])" (y**2/x**2),
                      t "(x[i]*x[i])*(x[i]*x[i])" (x**4)]
  where t str e = TestCase $ assertEqual str str (code e)
        x = r_var "x"
        y = s_var "y" :: Expression RealSpace

eqTests :: Test
eqTests = TestList [t "x*x == x**2" (x ** 2) (x*x),
                    t "isConstant 0 == Just 0" (Just 0) (isConstant (0 :: Expression RealSpace)),
                    t "isConstant 1 == Just 1" (Just 1) (isConstant (1 :: Expression RealSpace)),
                    t "isConstant 2 == Just 2" (Just 2) (isConstant (2 :: Expression RealSpace)),
                    t "(x*y)*4*(a*b) == 4*x*y*a*b" (4*x*y*a*b) ((x*y)*4*(a*b)),
                    t "-1 + 1 == 0" 0 (-1 + 1 :: Expression RealSpace),
                    t "(x*y)*(a*b) == x*y*a*b" (x*y*a*b) ((x*y)*(a*b)),
                    t "(x+y)+(a+b) == x+y+a+b" (x+y+a+b) ((x+y)+(a+b)),
                    t "(x*y*4)*(a*b) == 4*x*y*a*b" (4*x*y*a*b) ((x*y*4)*(a*b)),
                    t "(x*y*4)*(a*b*4) == 16*x*y*a*b" (16*x*y*a*b) ((x*y*4)*(a*b*4)),
                    t "4*(x*y) == 4*x*y" (4*x*y) (4*(x*y)),
                    t "2*(x + x) == 4*x" (4*x) (2*(x+x)),
                    t "2*(x + y) == 2*x + 2*y" (2*x + 2*y) (2*(x+y)),
                    t "2*(x/2 + y) == x + 2*y" (x + 2*y) (2*(x/2+y)),
                    t "derive x (x+y) == 1" 1 (derive (R "x") 1 (x+y)),
                    t "derive x (x+y) == 1" 1 (derive (R "x") 1 (x+y)),
                    t "derive x (x**2) == 2x" (2*x) (derive (R "x") 1 (x**2)),
                    t "derive x (k cos kx)" 
                       (latex $ kk * cos(kk*x))
                       (latex $ derive (R "x") 1 (sin (kk*x))),
                    t "derive x dV (sin kx)" 
                       (latex $ s_var "dV" * kk * cos(kk*x))
                       (latex $ derive (R "x") (s_var "dV") (sin (kk*x))),
                    t "derive x (sin x cos x) == 4x**3" 
                       (latex $ (cos x)**2 - (sin x)**2) 
                       (latex $ derive (R "x") 1 (sin x * cos x)),
                    t "derive x sin kx"
                       (latex $ (kk*cos (kk*x)))
                       (latex $ derive (R "x") 1 (sin (kk*x))),
                    t "derive x (A*sin kx)"
                       (latex $ (s_var "A"*kk*cos (kk*x)))
                       (latex $ derive (R "x") 1 (s_var "A"*sin (kk*x))),
                    t "derive x (sin kx cos kx)"
                       (latex $ kk*((cos (kk*x))**2 - (sin (kk*x))**2)) 
                       (latex $ derive (R "x") 1 (sin (kk*x) * cos (kk*x))),
                    t "derive x (kx cos kx)" 
                       (latex $ kk*cos(kk*x) - kk**2*x*sin(kk*x)) 
                       (latex $ derive (R "x") 1 (kk*x*cos(kk*x))),
                    t "derive x (sin kx) == k cos kx" 
                       (latex $ kk*cos(kk*x)) 
                       (latex $ derive (R "x") 1 (sin(kk*x))),
                    t "derive x (sin kx - kx cos kx)"
                       (latex $ kk**2*x*sin(kk*x)) 
                       (latex $ derive (R "x") 1 (sin (kk*x) - kk*x*cos(kk*x))),
                    t "derive x (sin kx) - derive x (kx cos kx)" 
                       (latex $ kk**2*x*sin(kk*x)) 
                       (latex $ derive (R "x") 1 (sin (kk*x)) - derive (R "x") 1 (kk*x*cos(kk*x))),
                    t "derive x (log x + x**2)"
                       (latex $ 1/x + 2*x -3*x**2 - 1/x**2)
                       (latex $ derive (R "x") 1 (log x + x**2 - x**3 + 1/x)),
                    t "makeHomogeneous (sin kr)" 0 (makeHomogeneous (sin kr)),
                    t "makeHomogeneous (cos kr)" 1 (makeHomogeneous (cos kr)),
                    t "makeHomogeneous (sin kr/k) == r" (s_var "r") (makeHomogeneous (sin kr/k)),
                    t "makeHomogeneous (1 - cos k) == 0" 0
                       (makeHomogeneous (1 - cos k)),
                    t "makeHomogeneous (r * cos kr - sin kr/k) == 0" 0
                       (makeHomogeneous (r * cos kr - sin kr/k)),
                    t "makeHomogeneous (r * cos kr) - makeHomogeneous (sin kr/k) == 0" 0
                       (makeHomogeneous (r * cos kr) - makeHomogeneous (sin kr/k)),
                    t "derive aaa (1 - sin aaa/aaa)" 
                       (latex $ - (cos aaa - sin aaa/aaa)/aaa)
                       (latex $ derive (R "aaa") 1 (1 - sin aaa/aaa)),
                    t "derive aaa (cos aaa - sin aaa/aaa)" 
                       (latex $ - sin aaa - (cos aaa - sin aaa/aaa)/aaa)
                       (latex $ derive (R "aaa") 1 (cos aaa - sin aaa/aaa)),
                    t "setZero aaa ((1 - sin aaa/aaa)/aaa) == 0" 0
                       (setZero (R "aaa") ((1 - sin aaa/aaa)/aaa)),
                    t "makeHomogeneous (ky*(r * cos kr - sin kr/k)/k**2) == 0" 0
                       (makeHomogeneous (ky*(r * cos kr - sin kr/k)/k**2)),
                    t "makeHomogeneous n2" (4*pi*(s_var "r")**2*s_var "x") (makeHomogeneous n2),
                    t "makeHomogeneous n2z" 0 (makeHomogeneous n2z),
                    t "makeHomogeneous n2y" 0 (makeHomogeneous n2y),
                    t "makeHomogeneous n2x" 0 (makeHomogeneous n2x),
                    t "makeHomogeneous vectorThirdTerm" 
                       ((4*pi*(s_var "r")**2 *s_var "x")**3) 
                       (makeHomogeneous vectorThirdTerm),
                    t "makeHomogeneous passing simplification"
                       0
                       (makeHomogeneous $
                        ((ifft(((-1*sin(kx)/kx + cos(kx)))/kx)*ifft((sin(kx))/kx)*(log(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))*(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))**2 + ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3))))/(ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3))**2*(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))**2))),
                    
                    t "makeHomogeneous derivative of sinc" 0
                       (makeHomogeneous $ (-1*sin(kx)/kx + cos(kx)) / kx),
                    
                    t "makeHomogeneous works unless you multiply by sinc'"
                       0
                       (makeHomogeneous $
                        (fft((ifft(((-1*sin(kx)/kx + cos(kx)))/kx)*ifft((sin(kx))/kx)*(log(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))*(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))**2 + ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3))))/(ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3))**2*(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))**2)))),
                    
                    t "makeHomogeneous with fft"
                       0
                       (makeHomogeneous $
                        (fft((ifft(((-1*sin(kx)/kx + cos(kx)))/kx)*ifft((sin(kx))/kx)*(log(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))*(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))**2 + ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3))))/(ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3))**2*(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))**2))
                         * (-1*sin(kx)/kx + cos(kx)))
                        / kx),
                    
                    t "makeHomogeneous of annoying fft ifft combo"
                       0
                       (makeHomogeneous $
                        (fft(2*(ifft(((-1*sin(kx)/kx + cos(kx)))/kx)*ifft((sin(kx))/kx)*(log(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))*(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))**2 + ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3))))/(ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3))**2*(1 + (-1)*ifft(((-1*kx*cos(kx) + sin(kx)))/(kx**3)))**2))
                         * (-1*sin(kx)/kx + cos(kx)))
                        / kx),
                    
                    t "makeHomogeneous (sin kr / k)" 1 (makeHomogeneous (sin kr/kr)),
                    t "makeHomogeneous (x+y) == x + y" (s_var "x"+s_var "y") (makeHomogeneous (x+y)),
                    t "makeHomogeneous (x**2) == x**2" (s_var "x"**2) (makeHomogeneous (x**2)),
                    t "x+0 == x" (x+0) x]
  where t str e1 e2 = TestCase $ assertEqual str e1 e2
        x = r_var "x"
        y = r_var "y"
        a = r_var "a"
        aaa = r_var "aaa"
        b = r_var "b"
        kk = s_var "k" :: Expression RealSpace
        kr = k*r
        r = s_var "r" :: Expression KSpace
        n2 = ifft ( smear * (4*pi) * r * (sin kr / k) * fft (r_var "x"))
        n2x = ifft ( smear * (4*pi) * i * kx*(r * cos kr - sin kr/k)/k**2 * fft (r_var "x"))
        n2y = ifft ( smear * (4*pi) * i * ky*(r * cos kr - sin kr/k)/k**2 * fft (r_var "x"))
        n2z = ifft ( smear * (4*pi) * i * kz*(r * cos kr - sin kr/k)/k**2 * fft (r_var "x"))
        vectorThirdTerm = n2*(n2**2 - 3*(n2x**2 + n2y**2 + n2z**2))
        smear = exp (-6.0*kdr*kdr)
        kdr = k*s_var "dr"
        i = s_var "complex(0,1)"

fftTests :: Test
fftTests = TestList [t "countFFT x = 0" 0 x,
                     t "countFFT nbar = 1" 2 nbar,
                     t "countFFT nbar = 1" 4 (nbar + log nbar)]
  where t str n e = TestCase $ assertEqual str n (countFFT $ fst $ simp e)
        x = r_var "x"
        spreading = 6.0
        kdr = k * s_var "dr"
        kR = k * s_var "R"
        nbar = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft (r_var "x"))

substitutionTests :: Test
substitutionTests = TestList [t x y (y**2) (x**2),
                              t x y (cos y**2) (cos x**2),
                              t (fft x) (k_var "temp_FFT") nbar_temp nbar,
                              t (k**2) (kk**2) nbarkk nbar,
                              t x y (integrate y) (integrate x),
                              t x y nbary nbar]
  where t a b eresult e = TestCase $ assertEqual (latex a ++ " -> " ++ latex b ++ " on " ++ latex e) (latex eresult) (latex $ substitute a b e)
        x = r_var "x"
        y = r_var "y" :: Expression RealSpace
        kk = k_var "k"
        spreading = 6.0
        dr = s_var "dr" :: Expression KSpace
        kdr = k * dr
        kR = k * s_var "R"
        nbar = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft x)
        nbar_temp = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * k_var "temp_FFT")
        nbarkk = ifft ( exp (-spreading*kk**2*dr**2) * (4*pi) * (sin(kk*s_var "R") - kk*s_var "R" * cos(kk*s_var "R")) / kk**3 * fft x)
        nbary = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft y)

main :: IO ()
main = do createDirectoryIfMissing True "tests/generated-haskell"
          writeFile "tests/generated-haskell/nice-sum.h" $ generateHeader (r_var "x" + s_var "kT") Nothing "NiceSum"
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
              -- nbar = ifft ( (4*pi) * (sin kR - kR * cos kR) / k**3 * fft (r_var "x"))
          writeFile "tests/generated-haskell/nice-nbar.h" $ generateHeader nbar (Just (r_var "R")) "NiceNbar"
          writeFile "tests/generated-haskell/nice-logoneminusnbar.h" $ generateHeader (log (1-ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft (r_var "x")))) (Just (r_var "R")) "NiceLogOneMinusNbar"
          writeFile "tests/generated-haskell/math.tex" $ latexfile [("nbar", nbar)]
          c <- runTestTT $ TestList [eqTests, codeTests, latexTests, fftTests, substitutionTests]
          if failures c > 0 || errors c > 0
            then fail $ "Failed " ++ show (failures c + errors c) ++ " tests."
            else do putStrLn "All tests passed!"
                    putStrLn $ show c

latexfile :: Type a => [(String, Expression a)] -> String
latexfile xs = "\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n\n" ++ unlines (map helper xs) ++ "\n\\end{document}"
  where helper (v,e) = v ++ "\n\\begin{equation}\n" ++ latex (substitute (k**2) (kk**2) e) ++ "\n\\end{equation}\n"
        kk = k_var "k"
\end{code}
