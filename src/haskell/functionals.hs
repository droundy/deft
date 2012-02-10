import CodeGen

main :: IO ()
main = 
  do let spreading = 6.0
         kdr = k * s_var "dr"
         rad :: Type a => Expression a
         rad = s_var "R"
         kR = k * rad
         i = s_var "complex(0,1)"
         smear = exp (-spreading*kdr*kdr)
         n3 = ifft ( smear * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft (r_var "x"))
         n2 = ifft ( smear * (4*pi) * rad * (sin kR / k) * fft (r_var "x"))
         n2x = ifft ( smear * (4*pi) * i * kx*(rad * cos kR - sin kR/k)/k**2 * fft (r_var "x"))
         n2y = ifft ( smear * (4*pi) * i * ky*(rad * cos kR - sin kR/k)/k**2 * fft (r_var "x"))
         n2z = ifft ( smear * (4*pi) * i * kz*(rad * cos kR - sin kR/k)/k**2 * fft (r_var "x"))
         vectorThirdTerm = n2*(n2**2 - 3*(n2x**2 + n2y**2 + n2z**2))
         phi1 = (-1/(4*pi*rad**2))*n2*log(1-n3)
         phi2 = (n2**2 - n2x**2 - n2y**2 - n2z**2)/(4*pi*rad*(1-n3))
         phi3 = (n3 + (1-n3)**2*log(1-n3))/(36*pi * n3**2 * (1-n3)**2)*vectorThirdTerm
     writeFile "src/HardSpheresNoTensor2Fast.cpp" $ 
       generateHeader (s_var "kT" * (phi1+phi2+phi3)) (Just (r_var "R")) "HardSpheresNoTensor2"
    -- compile with make src/HardSpheresNoTensor2Fast.o
