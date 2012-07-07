module WhiteBear
       ( whitebear,
         phi1, phi2, phi3,
         xshell, yshell, zshell )
       where

import Expression

kT :: Type a => Expression a
kT = s_tex "kT" "kT"

rad :: Type a => Expression a
rad = s_var "R"

kdr :: Expression KSpace
kdr = k * s_tex "dr" "\\Delta r"

kR :: Expression KSpace
kR = k * rad

smear :: Expression KSpace
smear = exp (-6.0*kdr*kdr)

n, n3, n2, n1, n0, n2x, n2y, n2z, n1v_dot_n2v, sqr_n2v :: Expression RealSpace
n = "n" === r_var "x"
n3 = "n3" === step n
n2 = "n2" === shell n
n1 = "n1" === n2 / (4*pi*rad)
n0 = "n0" === n2 / (4*pi*rad**2)
n2x = "n2x" === xshell n
n2y = "n2y" === yshell n
n2z = "n2z" === zshell n
sqr_n2v = var "n2vsqr" "{\\vec{n}_{1v}\\cdot\\vec{n}_{2v}}" (n2x**2 + n2y**2 + n2z**2)
n1v_dot_n2v = {- var "n1vdn2v" "{\\left|\\vec{n}_{2v}\\right|}" -} (sqr_n2v/(4*pi*rad))

shell, step, xshell, yshell, zshell :: Expression RealSpace -> Expression RealSpace
shell x = ifft ( deltak * fft x)
  where deltak = protect "deltak" "\\delta(k)" $ smear * (4*pi) * rad * sin kR / k
step x = ifft ( stepk * fft x)
  where stepk = protect "step" "\\Theta(k)" $ smear * (4*pi) * (sin kR - kR * cos kR) / k**3
xshell x = ifft ( deltax * fft x)
  where deltax = protect "deltax" "\\delta_x(k)" $
                 smear * (4*pi) * imaginary * kx*(rad * cos kR - sin kR/k)/k**2
yshell x = ifft ( deltay * fft x)
  where deltay = protect "deltay" "\\delta_y(k)" $
                 smear * (4*pi) * imaginary * ky*(rad * cos kR - sin kR/k)/k**2
zshell x = ifft ( deltaz * fft x)
  where deltaz = protect "deltaz" "\\delta_z(k)" $
                 smear * (4*pi) * imaginary * kz*(rad * cos kR - sin kR/k)/k**2

vectorThirdTerm :: Expression RealSpace
vectorThirdTerm = n2*(n2**2 - 3*sqr_n2v)

phi1, phi2, phi3 :: Expression RealSpace
phi1 = var "phi1" "\\Phi_1" $ -n0*log(1-n3)
phi2 = var "phi2" "\\Phi_2" $ (n2*n1 - n1v_dot_n2v)/(1-n3)
phi3 = var "phi3" "\\Phi_3" $ (n3 + (1-n3)**2*log(1-n3))/(36*pi * n3**2 * (1-n3)**2)*vectorThirdTerm


whitebear :: Expression Scalar
whitebear = var "whitebear" "F_{HS}" $ integrate (kT*(phi1+phi2+phi3))
