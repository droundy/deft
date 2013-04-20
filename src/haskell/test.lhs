\begin{code}
import System.Directory ( createDirectoryIfMissing )
import CodeGen
import Optimize ( optimize, findToDo )
import Latex ( latexOptimizedExpression )
import Test.HUnit
import FMT ( shell )
import SFMT ( n1 )
import System.Environment ( getArgs )
import qualified Data.Set as Set

codeTests :: Test
codeTests = TestList [t "x[i]" x,
                      t "0" (0 :: Expression RealSpace),
                      t "sqrt(x[i])" (sqrt x),
                      t "x[i] + -3.0*(x[i]*x[i]*x[i] + x[i]*x[i])" foo,
                      t "-3.0*x[i]*x[i]*x[i] + -3.0*x[i]*x[i] + x[i]" (cleanvars foo),
                      t "x[i]*x[i]" (x**2),
                      t "1/(x[i]*x[i])" (1/x**2),
                      t "y*y/(x[i]*x[i])" (y**2/x**2),
                      t "(x[i]*x[i])*(x[i]*x[i])" (x**4)]
  where t str e = TestCase $ assertEqual str str (code e)
        x = r_var "x"
        y = s_var "y" :: Expression RealSpace
        foo = x - 3* var "foo" "foo" (x**2 + x**3)

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
                    t "derive x (x+y) == 1" 1 (derive x 1 (x+y)),
                    t "derive x (x+y) == 1" 1 (derive ("xx" === x**2) 1 (("xx" === x**2)+y)),
                    t "derive x (x+y) == 1" 1 (derive x 1 (x+y)),
                    t "derive x (x**2) == 2x" (2*x) (derive x 1 (x**2)),
                    t "derive x (k cos kx)" 
                       (kk * cos(kk*x))
                       (derive x 1 (sin (kk*x))),
                    t "derive x dV (sin kx)" 
                       (s_var "dV" * kk * cos(kk*x))
                       (derive x (s_var "dV") (sin (kk*x))),
                    t "derive x (sin x cos x) == 4x**3" 
                       ((cos x)**2 - (sin x)**2) 
                       (derive x 1 (sin x * cos x)),
                    t "derive x sin kx"
                       ((kk*cos (kk*x)))
                       (derive x 1 (sin (kk*x))),
                    t "derive x (A*sin kx)"
                       ((s_var "A"*kk*cos (kk*x)))
                       (derive x 1 (s_var "A"*sin (kk*x))),
                    t "derive x (kx cos kx)" 
                       (kk*cos(kk*x) - kk**2*x*sin(kk*x)) 
                       (derive x 1 (kk*x*cos(kk*x))),
                    t "derive x (sin kx) == k cos kx" 
                       (kk*cos(kk*x)) 
                       (derive x 1 (sin(kk*x))),
                    t "derive x (sin kx - kx cos kx)"
                       (kk**2*x*sin(kk*x)) 
                       (derive x 1 (sin (kk*x) - kk*x*cos(kk*x))),
                    t "derive x (sin kx) - derive x (kx cos kx)" 
                       (kk**2*x*sin(kk*x)) 
                       (derive x 1 (sin (kk*x)) - derive x 1 (kk*x*cos(kk*x))),
                    t "derive x (log x + x**2)"
                       (1/x + 2*x -3*x**2 - 1/x**2)
                       (derive x 1 (log x + x**2 - x**3 + 1/x)),
                    --t "derive with named subexpression"
                    --   (cleanvars $ derive x 1 (log x + log x * ("xx" === x**2 - x**3) + 1/x))
                    --   (cleanvars $ derive x 1 (log x + log x * (x**2 - x**3) + 1/x)),
                    t "makeHomogeneous (sin kr)" 0 (makeHomogeneous (sin kr)),
                    t "makeHomogeneous (cos kr)" 1 (makeHomogeneous (cos kr)),
                    t "makeHomogeneous (sin kr/k) == r" (s_var "r") (makeHomogeneous (sin kr/k)),
                    t "makeHomogeneous (1 - cos k) == 0" 0
                       (makeHomogeneous (1 - cos k)),
                    t "makeHomogeneous (r * cos kr - sin kr/k) == 0" 0
                       (makeHomogeneous (r * cos kr - sin kr/k)),
                    t "makeHomogeneous (r * cos kr) - makeHomogeneous (sin kr/k) == 0" 0
                       (makeHomogeneous (r * cos kr) - makeHomogeneous (sin kr/k)),
                    t "setZero aaa ((1 - sin aaa/aaa)/aaa) == 0" 0
                       (setZero (ER $ r_var "aaa") ((1 - sin aaa/aaa)/aaa)),
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
        n2 = ifft ( smear * (4*pi) * r * (sin kr / k) * fft x)
        n2x = ifft ( smear * (4*pi) * i * kx*(r * cos kr - sin kr/k)/k**2 * fft x)
        n2y = ifft ( smear * (4*pi) * i * ky*(r * cos kr - sin kr/k)/k**2 * fft x)
        n2z = ifft ( smear * (4*pi) * i * kz*(r * cos kr - sin kr/k)/k**2 * fft x)
        vectorThirdTerm = n2*(n2**2 - 3*(n2x**2 + n2y**2 + n2z**2))
        smear = exp (-6.0*kdr*kdr)
        kdr = k*s_var "dr"
        i = s_var "complex(0,1)"

fftTests :: Test
fftTests = TestList [t "countFFT x = 0" 0 x,
                     t "countFFT nbar" 2 n3,
                     t "countFFT nbar + n2" 2 (nbar + n2),
                     t "countFFT n0raw log n3" 3 (n0raw*log n3),
                     t "countFFT n0 log n3" 3 (n0*log n3),
                     t "countFFT derive n0raw log n3" 6 (gradme $ integrate $ kT*n0raw*log n3),
                     t "countFFT derive n0 log n3" 6 (gradme $ integrate $ kT*n0*log n3),
                     t "countFFT derive assocalike n0raw" 6
                           (gradme $ integrate $ assocalike n0raw),
                     t "countFFT derive assocalike n0" 6
                           (gradme $ integrate $ assocalike n0),
                     t "countFFT n3 + n2a" 2 (n3 + n2a),
                     t "countFFT nbar*n2 + nbar" 3 (nbar*n2 + nbar),
                     t "countFFT nbar + log nbar" 2 (nbar + log nbar)]
  where t str nn e = TestCase $ assertEqual str nn
                     (countFFT $ fst $ optimize [mkExprn $ factorize $ joinFFTs $ cleanvars e])
        x = r_var "x"
        spreading = 6.0
        kdr = k * s_var "dr"
        kR = k * s_var "R"
        nbar = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft x)
        n2 = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * s_var "R" * (sin kR / k) * fft x)
        n3 = "n3" === nbar
        n2a = "n2" === n2
        n0raw = n2 / (4*pi*s_var "R"**2)
        n0 = "n0" === n0raw
        gradme :: Expression Scalar -> Expression RealSpace
        gradme = cleanvars . derive (r_var "x") 1
        kT = s_var "kT"
        assocalike nn = nn*(1-n3)*log(nn*n2a*(1 - n3))

joinFFTtests :: Test
joinFFTtests = TestList [t "joinFFT fft(a)*k + fft(b)*k = fft(a+b)*k" (fft (a+b)*k) (fft a*k + fft b*k),
                         t "joinFFT ifft(fft(a)*k + fft(b)*k) = ifft(fft(a+b)*k)"
                         (ifft (fft (a+b)*k))
                         (ifft (fft a*k + fft b*k)),
                         t "joinFFT integrate(fft(a)*k + fft(b)*k) = integrate(fft(a+b)*k)"
                         (integrate $ ifft (fft (a+b)*k) :: Expression Scalar)
                         (integrate $ ifft (fft a*k + fft b*k)),
                         t "joinFFT fft(a) + fft(b) = fft(a+b)" (fft (a+b)) (fft a + fft b)]
  where t str nn e = TestCase $ assertEqual str nn (joinFFTs e)
        a = r_var "a"
        b = r_var "b"

memTests :: Test
memTests = TestList [t "peakMem x = 0" 0 x,
                     t "peakMem nbar" 1 n3,
                     t "peakMem x1 + x2 and other stuff" 2 (x1 + x2 + cos(x1 + x2) + ifft ( ksqr * fft (x1 + x2 + 5) )),
                     t "peakMem nbar + n2" 1 (nbar + n2),
                     t "peakMem n0raw log n3" 3 (n0raw*log n3), -- was 2
                     t "peakMem n0 log n3" 3 (n0*log n3), --was 2
                     t "peakMem derive n0raw log n3" 4 (gradme $ integrate $ kT*n0raw*log n3),
                     t "peakMem derive n0 log n3" 4 (gradme $ integrate $ kT*n0*log n3),
                     t "peakMem derive assocalike n0raw" 4
                           (gradme $ integrate $ assocalike n0raw),
                     t "peakMem derive assocalike n0" 4
                           (gradme $ integrate $ assocalike n0),
                     t "peakMem n3 + n2a" 1 (n3 + n2a),
                     t "peakMem nbar*n2 + nbar" 3 (nbar*n2 + nbar), --was 2
                     t "peakMem nbar + log nbar" 1 (nbar + log nbar)]
  where t str nn e = TestCase $ assertEqual str nn (peakMem $ reuseVar $ freeVectors $ fst $ optimize
                                                    [mkExprn $ factorize $ joinFFTs $ cleanvars e])
        x = r_var "x"
        spreading = 6.0
        kdr = k * s_var "dr"
        kR = k * s_var "R"
        nbar = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft x)
        n2 = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * s_var "R" * (sin kR / k) * fft x)
        n3 = "n3" === nbar
        n2a = "n2" === n2
        n0raw = n2 / (4*pi*s_var "R"**2)
        n0 = "n0" === n0raw
        gradme :: Expression Scalar -> Expression RealSpace
        gradme = cleanvars . derive (r_var "x") 1
        kT = s_var "kT"
        assocalike nn = nn*(1-n3)*log(nn*n2a*(1 - n3))
        x1 = r_var "x1"
        x2 = r_var "x2"

substitutionTests :: Test
substitutionTests = TestList [t x y (y**2) (x**2),
                              t x y (cos y**2) (cos x**2),
                              t xx xy (xy**2) (xx**2),
                              t (fft x) (k_var "temp_FFT") nbar_temp nbar,
                              t (k**2) (kk**2) nbarkk nbar,
                              t (r_var "n2x") (shell x) (shell x**2) (r_var "n2x" ** 2),
                              t (r_var "n2x") n1 (r_var "n2y"**2 + n1**2) (r_var "n2y"**2 + r_var "n2x"**2),
                              t x y (integrate y :: Expression Scalar) (integrate x),
                              t x y nbary nbar]
  where t a b eresult e = TestCase $ assertEqual (latex a ++ " -> " ++ latex b ++ "\non\n" ++ latex e) 
                                                 (eresult) (substitute a b e)
        x = r_var "x"
        xx = "xx" === x**2
        xy = "xy" === x*y
        y = r_var "y" :: Expression RealSpace
        kk = k_var "k"
        spreading = 6.0
        kdr = k * dr
        kR = k * s_var "R"
        nbar = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft x)
        nbar_temp = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * k_var "temp_FFT")
        nbarkk = ifft ( exp (-spreading*kk**2*dr**2) * (4*pi) * (sin(kk*s_var "R") - kk*s_var "R" * cos(kk*s_var "R")) / kk**3 * fft x)
        nbary = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft y)

hasexpressionTests :: Test
hasexpressionTests = TestList [t x1 (x1 + x2) True,
                               t x4 (x1 + x2) False,
                               t (x1 + x2) x True,
                               t (x1 + x3) (x1 + x2 + x3) True,
                               t (x3 + x1) (x1 + x2 + x3) True,
                               t (x3 + x1) (x1 * 2 + x2 + x3 * 2) True,
                               t (x1 + x2) (x1 * 2 + x2 * 3 + x3 * 2) False,
                               t (x1**2+x2**2)
                                     ((x1**2+x2**2)*rad + 3 + x4*log(1- x3)*(-3*x1**2-3*x2**2+x4**2)/(x3**2*(1-x3)**2)/36/pi)
                                     True,
                               t (x1**2+2*x2**2)
                                     ((x1**2+x2**2)*rad + 3 + x4*log(1- x3)*(-3*x1**2-3*x2**2+x4**2)/(x3**2*(1-x3)**2)/36/pi)
                                     False,
                               t (x1+x2) (x1+x2+cos(2*x1+2*x2)+x3) True,
                               t (x1 + x3) x False,
                               t (x4 + x5) x True,
                               t x2 x True]
  where t a b e = TestCase $ assertEqual ("hasexpression failure " ++ latex a ++ " in " ++ latex b)
                                         e (hasexpression a b)
        x1 = r_var "x1"
        x2 = r_var "x2"
        x3 = r_var "x3"
        x4 = r_var "x4"
        x5 = r_var "x5"
        x = x1 + x2 + x3 * (x4 + x5)
        rad = s_var "R"

findToDoTests :: Test
findToDoTests = TestList [t (Just $ ER $ -3*(x1**2+x2**2))
                            ((x1**2+x2**2)*rad + 3 + x4*log(1- x3)*(-3*x1**2-3*x2**2+x4**2)/(x3**2*(1-x3)**2)/36/pi),
                          --FIXME
                          --t (DoR $ -3*x1-3*x2)
                          --  ((x1+x2)*rad + log(x1+x2)),
                          t Nothing (x1**2/4/pi+x2**2/4/pi + 3 + cos(-3*x1**2-2*x2**2)),
                          t (Just $ ER $ x1+x2) (x1+x2 + x3 + x4 + cos(x1+x2) + rad)]
    where t ee e = TestCase $ assertEqual ("findToDo" ++ latex e) ee
                   (findToDo Set.empty [mkExprn e] e)
          x1 = r_var "rtemp_1"
          x2 = r_var "rtemp_2"
          x3 = r_var "rtemp_3"
          x4 = r_var "rtemp_4"
          rad = r_var "R"

multisubstituteTests :: Test
multisubstituteTests = TestList [t x1 x2 (x2+x3) (x1+x3),
                                 -- t k kv (kv**2) (k**2),
                                 t (k**2) (kv**2) (kv**2) (k**2),
                                 t k kv (kv**3) (k**3),
                                 t (x1+x2) x4 (x4+x3) (x1+x2+x3),
                                 t (x1+x2) x4 (2*x4+x3) (2*x1+2*x2+x3),
                                 t (x1**2+x2**2) z
                                       (z*rad + 3 + x4*log(1- x3)*(-3*z+x4**2)/(x3**2*(1-x3)**2)/36/pi)
                                       ((x1**2+x2**2)*rad + 3 + x4*log(1- x3)*(-3*x1**2-3*x2**2+x4**2)/(x3**2*(1-x3)**2)/36/pi),
                                 t (-x1**2-x2**2) x4 (-x4/4/pi+cos((3*x4+x5)/r_var "R")+x3) (x1**2/4/pi+x2**2/4/pi+cos((-3*x1**2+x5-3*x2**2)/r_var "R")+x3),
                                 t (x1+x2) x4 (x4+cos(x4)+x3) (x1+x2+cos(x1+x2)+x3),
                                 t (x1+x2) x4 (x4+cos(2*x4)+x3) (x1+x2+cos(2*x1+2*x2)+x3),
                                 t (x4+x5) x3 (x1+x2+x3*x3) x,
                                 t (x1+x2) z (z*x3+(z+x4)*x3+x5*(x2+x3)) y,
                                 t (x1+x2) z (3*z + x3 + cos(z*(z+x3))) (3*x1+3*x2 + x3 + cos(z*(x1+x2+x3))),
                                 t (x1*x2) x4 (x4*x3+x2*x3+x1*x3) (x1*x2*x3+x2*x3+x1*x3),
                                 t (sqrt(x1*x2)) x4 (x4**4*x3+x2*x3+x4**2*x3) (x1**2*x2**2*x3+x2*x3+x1*x2*x3),
                                 t (x1*x2) x4 (x4**2*x3+x2*x3+x4*x3) (x1**2*x2**2*x3+x2*x3+x1*x2*x3)]
  where t a b eresult e = TestCase $ assertEqual (latex a ++ " -> " ++ latex b ++ "\non\n" ++ latex e) 
                                                 (eresult) (substitute a b e)
        x1 = r_var "x1"
        x2 = r_var "x2"
        x3 = r_var "x3"
        x4 = r_var "x4"
        x5 = r_var "x5"
        kv = k_var "k"
        z = r_var "z"
        rad = s_var "R"
        x = x1 + x2 + x3 * (x4 + x5)
        y = (x1 + x2) * x3 + x3 * (x1 + x2 + x4) + x5 * (x2 + x3)

main :: IO ()
main = do createDirectoryIfMissing True "tests/generated-haskell"
          args <- getArgs
          let x = r_var "x"
              wf fn s = if "codegen" `elem` args
                        then do putStrLn $ "Creating file " ++ fn
                                writeFile fn s
                        else return ()
          wf "tests/generated-haskell/nice-sum.h" $ generateHeader (r_var "x" + s_var "kT") [] "NiceSum"
          wf "tests/generated-haskell/nice-quad.h" $ generateHeader ((r_var "x" + s_var"kT")**2 - r_var "x" + 2*s_var "kT") [] "NiceQuad"
          wf "tests/generated-haskell/nice-sqrt.h" $ generateHeader (r_var "x"**0.5) [] "NiceSqrt"
          wf "tests/generated-haskell/nice-sqrtandmore.h" $ generateHeader ((r_var "x"**0.5) - r_var "x" + 2*s_var "kT") [] "NiceSqrtandMore"
          wf "tests/generated-haskell/nice-log.h" $ generateHeader (log x) [] "NiceLog"
          wf "tests/generated-haskell/nice-logandsqr.h" $ generateHeader (log x + x**2) [] "NiceLogandSqr"
          wf "tests/generated-haskell/nice-logandsqrandinverse.h" $ generateHeader (log x + x**2 - x**3 + x**(-1)) [] "NiceLogandSqrandInverse"
          wf "tests/generated-haskell/nice-logoneminusx.h" $ generateHeader (log (1 - r_var "x")) [] "NiceLogOneMinusX"
          let spreading = 6.0
              kdr = k * dr
              kR = k * rad
              rad :: Type a => Expression a
              rad = s_var "R"
              nbar = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft x)
              n3 = nbar
              n2 = ifft ( exp (-spreading*kdr*kdr) * (4*pi) * rad * (sin kR/k) * fft x)
          wf "tests/generated-haskell/nice-nbar.h" $ generateHeader nbar ["R"] "NiceNbar"
          wf "tests/generated-haskell/nice-n2.h" $ generateHeader n2 ["R"] "NiceN2"
          wf "tests/generated-haskell/nice-logoneminusnbar.h" $
            generateHeader (log (1-nbar)) ["R"] "NiceLogOneMinusNbar"
          wf "tests/generated-haskell/nice-phi1.h" $
            generateHeader (-n2*log(1-n3)) ["R"] "NicePhi1"
          wf "tests/generated-haskell/nice-phi1.h" $ 
            generateHeader (-n2*log(1-n3)/(4*pi*rad**2)) ["R"] "NicePhi1"
          let smear = exp (-spreading*kdr*kdr)
              n2x = ifft ( smear * (4*pi) * imaginary * kx*(rad * cos kR - sin kR/k)/k**2 * fft x)
              n2y = ifft ( smear * (4*pi) * imaginary * ky*(rad * cos kR - sin kR/k)/k**2 * fft x)
              n2z = ifft ( smear * (4*pi) * imaginary * kz*(rad * cos kR - sin kR/k)/k**2 * fft x)
              phi2here = (n2**2 - n2x**2 - n2y**2 - n2z**2)/(1-n3)/(4*pi*rad)
          wf "tests/generated-haskell/nice-phi2.h" $
            generateHeader phi2here ["R"] "NicePhi2"
          wf "tests/generated-haskell/nice-phi3.h" $
            generateHeader ((n3 + (1-n3)**2*log(1-n3))/(36*pi* n3**2 * (1-n3)**2)*n2*(n2**2 - 3*(n2x**2+n2y**2+n2z**2)))
            ["R"] "NicePhi3"
          wf "tests/generated-haskell/nice-n2xsqr.h" $
            generateHeader (n2x**2) ["R"] "NiceN2xsqr"
          wf "tests/generated-haskell/math.tex" $ latexfile [("n3", n3), ("n2", n2), ("n2x", n2x),
                                                                    ("grad n2xsqr", derive x 1 (n2x**2))]
          wf "tests/generated-haskell/whitebear.tex" $ latexOptimizedExpression $ factorize $ joinFFTs $ derive (r_var "x") (r_var "ingrad") $
            substitute k (k_var "k") $
            (n3 + (1-n3)**2*log(1-n3))/(36*pi* n3**2 * (1-n3)**2)*n2*(n2**2 - 3*(n2x**2+n2y**2+n2z**2))
          if "codegen" `elem` args
            then putStrLn "Not running actual tests, just generating test code...\n"
            else do c <- runTestTT $ TestList [joinFFTtests,
                                               substitutionTests, hasexpressionTests, multisubstituteTests,
                                               findToDoTests, eqTests, codeTests, fftTests, memTests]
                    if failures c > 0 || errors c > 0
                      then fail $ "Failed " ++ show (failures c + errors c) ++ " tests."
                      else do putStrLn "All tests passed!"
                              putStrLn $ show c

latexfile :: Type a => [(String, Expression a)] -> String
latexfile xs = "\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n\n" ++ unlines (map helper xs) ++ "\n\\end{document}"
  where helper (v,e) = v ++ "\n\\begin{equation}\n" ++ latex (cleanup e) ++ "\n\\end{equation}\n"
        kk = k_var "k"
        x = r_var "x"
        cleanup = substitute (k**2) (kk**2) . substitute (shell x) (r_var "n2x")

\end{code}


