\begin{code}
import System.Directory ( createDirectoryIfMissing )
import CodeGen
import Test.HUnit
import SomeFunctionals

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
                    --t "isConstant 0 == Just 0" (Just 0) (isConstant (0 :: Expression RealSpace)),
                    --t "isConstant 1 == Just 1" (Just 1) (isConstant (1 :: Expression RealSpace)),
                    --t "isConstant 2 == Just 2" (Just 2) (isConstant (2 :: Expression RealSpace)),
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
                    t "factorandsum yx+yx^2" (y*x*(1+x*y)) (factorandsum [y*x,y**2 * x**2]),
                    t "factorandsum yx+2yx^2" (y*x*(1+2*x*y)) (factorandsum [y*x,2 * y**2 * x**2]),
                    t "derive x y(x+yx^2)" (y*(1+2*x*y)) (derive (R "x") 1 (y*(x+y*x**2))),
                    t "derive x (x+y) == 1" 1 (derive (R "x") 1 (x+y)),
                    t "derive x (x+y) == 1" 1 (derive (R "x") 1 (x+y)),
                    t "derive x (x**2) == 2x" (2*x) (derive (R "x") 1 (x**2)),
                    t "derive x (k cos kx)" 
                       (kk * cos(kk*x))
                       (derive (R "x") 1 (sin (kk*x))),
                    t "derive x dV (sin kx)" 
                       (s_var "dV" * kk * cos(kk*x))
                       (derive (R "x") (s_var "dV") (sin (kk*x))),
                    t "derive x (sin x cos x) == 4x**3" 
                       ((cos x)**2 - (sin x)**2) 
                       (derive (R "x") 1 (sin x * cos x)),
                    t "derive x sin kx"
                       ((kk*cos (kk*x)))
                       (derive (R "x") 1 (sin (kk*x))),
                    t "derive x (A*sin kx)"
                       ((s_var "A"*kk*cos (kk*x)))
                       (derive (R "x") 1 (s_var "A"*sin (kk*x))),
                    t "derive x (sin kx cos kx)"
                       (kk*((cos (kk*x))**2 - (sin (kk*x))**2)) 
                       (derive (R "x") 1 (sin (kk*x) * cos (kk*x))),
                    t "derive x (kx cos kx)" 
                       (kk*cos(kk*x) - kk**2*x*sin(kk*x)) 
                       (derive (R "x") 1 (kk*x*cos(kk*x))),
                    t "derive x (sin kx) == k cos kx" 
                       (kk*cos(kk*x)) 
                       (derive (R "x") 1 (sin(kk*x))),
                    t "derive x (sin kx - kx cos kx)"
                       (kk**2*x*sin(kk*x)) 
                       (derive (R "x") 1 (sin (kk*x) - kk*x*cos(kk*x))),
                    t "derive x (sin kx) - derive x (kx cos kx)" 
                       (kk**2*x*sin(kk*x)) 
                       (derive (R "x") 1 (sin (kk*x)) - derive (R "x") 1 (kk*x*cos(kk*x))),
                    t "derive x (log x + x**2)"
                       (1/x + 2*x -3*x**2 - 1/x**2)
                       (derive (R "x") 1 (log x + x**2 - x**3 + 1/x)),
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
                       (- (cos aaa - sin aaa/aaa)/aaa)
                       (derive (R "aaa") 1 (1 - sin aaa/aaa)),
                    t "derive aaa (cos aaa - sin aaa/aaa)" 
                       (- sin aaa - (cos aaa - sin aaa/aaa)/aaa)
                       (derive (R "aaa") 1 (cos aaa - sin aaa/aaa)),
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
  where t :: Type a => String -> Expression a -> Expression a -> Test
        --t str e1 e2 = TestCase $ assertEqual str (latex e1) (latex e2)
        t str e1 e2 = TestCase $ assertEqual str e1 e2
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
                     t "countFFT nbar" 2 nbar,
                     t "countFFT nbar*n2 + nbar" 3 (nbar*n2 + nbar),
                     t "countFFT nbar + log nbar" 2 (nbar + log nbar)]
  where t str nn e = TestCase $ assertEqual str nn (countFFT $ fst $ simp2 e)
        x = r_var "x"
        spreading = 6.0
        kdr = k * s_var "dr"
        kR = k * s_var "R"
        nbar = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft (r_var "x"))
        n2 = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * s_var "R" * (sin kR / k) * fft (r_var "x"))

substitutionTests :: Test
substitutionTests = TestList [t x y (y**2) (x**2),
                              t x y (cos y**2) (cos x**2),
                              t (fft x) (k_var "temp_FFT") nbar_temp nbar,
                              t (k**2) (kk**2) nbarkk nbar,
                              t (r_var "n2x") (xshell x) (xshell x**2) (r_var "n2x" ** 2),
                              t (r_var "n2x") (xshell x) (r_var "n2y"**2 + xshell x**2) (r_var "n2y"**2 + r_var "n2x"**2),
                              t x y (integrate y) (integrate x),
                              t x y nbary nbar]
  where t a b eresult e = TestCase $ assertEqual (latex a ++ " -> " ++ latex b ++ "\non\n" ++ latex e) 
                                                 (eresult) (substitute a b e)
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
          writeFile "tests/generated-haskell/nice-sum.h" $ generateHeader (r_var "x" + s_var "kT") [] "NiceSum"
          writeFile "tests/generated-haskell/nice-quad.h" $ generateHeader ((r_var "x" + s_var"kT")**2 - r_var "x" + 2*s_var "kT") [] "NiceQuad"
          writeFile "tests/generated-haskell/nice-sqrt.h" $ generateHeader (r_var "x"**0.5) [] "NiceSqrt"
          writeFile "tests/generated-haskell/nice-sqrtandmore.h" $ generateHeader ((r_var "x"**0.5) - r_var "x" + 2*s_var "kT") [] "NiceSqrtandMore"
          writeFile "tests/generated-haskell/nice-log.h" $ generateHeader (log (r_var "x")) [] "NiceLog"
          writeFile "tests/generated-haskell/nice-logandsqr.h" $ generateHeader (log (r_var "x") + (r_var "x")**2) [] "NiceLogandSqr"
          writeFile "tests/generated-haskell/nice-logandsqrandinverse.h" $ generateHeader (log (r_var "x") + (r_var "x")**2 - (r_var "x")**3 + (r_var "x")**(-1)) [] "NiceLogandSqrandInverse"
          writeFile "tests/generated-haskell/nice-logoneminusx.h" $ generateHeader (log (1 - r_var "x")) [] "NiceLogOneMinusX"
          let spreading = 6.0
              kdr = k * s_var "dr"
              kR = k * rad
              rad :: Type a => Expression a
              rad = s_var "R"
              nbar = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft (r_var "x"))
              n3 = nbar
              n2 = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * rad * (sin kR/k) * fft (r_var "x"))
          writeFile "tests/generated-haskell/nice-nbar.h" $ generateHeader nbar ["R"] "NiceNbar"
          writeFile "tests/generated-haskell/nice-n2.h" $ generateHeader n2 ["R"] "NiceN2"
          writeFile "tests/generated-haskell/nice-logoneminusnbar.h" $ 
            generateHeader (log (1-nbar)) ["R"] "NiceLogOneMinusNbar"
          writeFile "tests/generated-haskell/nice-phi1.h" $ 
            generateHeader (-n2*log(1-n3)/(4*pi*rad**2)) ["R"] "NicePhi1"
          let smear = exp (-spreading*kdr*kdr)
              i = s_var "complex(0,1)"
              x = r_var "x"
              n2x = ifft ( smear * (4*pi) * i * kx*(rad * cos kR - sin kR/k)/k**2 * fft x)
              n2y = ifft ( smear * (4*pi) * i * ky*(rad * cos kR - sin kR/k)/k**2 * fft x)
              n2z = ifft ( smear * (4*pi) * i * kz*(rad * cos kR - sin kR/k)/k**2 * fft x)
              phi2here = (n2**2 - n2x**2 - n2y**2 - n2z**2)/(1-n3)/(4*pi*rad)
          writeFile "tests/generated-haskell/nice-phi2.h" $ 
            generateHeader phi2here ["R"] "NicePhi2"
          writeFile "tests/generated-haskell/nice-phi3.h" $ 
            generateHeader ((n3 + (1-n3)**2*log(1-n3))/(36*pi* n3**2 * (1-n3)**2)*n2*(n2**2 - 3*(n2x**2+n2y**2+n2z**2))) 
            ["R"] "NicePhi3"
          writeFile "tests/generated-haskell/nice-n2xsqr.h" $
            generateHeader (n2x**2) ["R"] "NiceN2xsqr"
          writeFile "tests/generated-haskell/math.tex" $ latexfile [("n3", n3), ("n2", n2), ("n2x", n2x), 
                                                                    ("grad n2xsqr", derive (R "x") 1 (n2x**2))]
          c <- runTestTT $ TestList [eqTests, codeTests, latexTests, fftTests, substitutionTests]
          if failures c > 0 || errors c > 0
            then fail $ "Failed " ++ show (failures c + errors c) ++ " tests."
            else do putStrLn "All tests passed!"
                    putStrLn $ show c

latexfile :: Type a => [(String, Expression a)] -> String
latexfile xs = "\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n\n" ++ unlines (map helper xs) ++ "\n\\end{document}"
  where helper (v,e) = v ++ "\n\\begin{equation}\n" ++ latex (cleanup e) ++ "\n\\end{equation}\n"
        kk = k_var "k"
        x = r_var "x"
        cleanup = substitute (k**2) (kk**2) . substitute (xshell x) (r_var "n2x")
\end{code}
