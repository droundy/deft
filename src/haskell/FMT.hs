{- FMT defines the fundamental measures that are the raw input of
fundamental measure theory.  The FMT free energy functional itself is
defined in WhiteBear (or other modules for other versions thereof). -}

module FMT
       ( n, n3, n2, n2p, n1, n0, rad,
         n2x, n2y, n2z, n1x, n1y, n1z,
         n2px, n2py, n2pz,
         sqr_n2v, n1v_dot_n2v,
         shell, shell_diam, step, xshell, yshell, zshell,
         shellPrime, xshellPrime, yshellPrime, zshellPrime,
         dr )
       where

import Expression

rad :: Type a => Expression a
rad = s_var "R"

kdr :: Expression KSpace
kdr = k * dr

kR :: Expression KSpace
kR = k * rad

smear :: Expression KSpace
smear = exp (-6.0*kdr*kdr)

n, n3, n2, n1, n0, n1x, n1y, n1z, n2x, n2y, n2z, n2p, n2px, n2py, n2pz, n1v_dot_n2v, sqr_n2v :: Expression RealSpace
n = "n" === r_var "x"
n3 = "n3" === step n
n2 = "n2" === shell n
n2p = "n2p" === shellPrime n
n1 = "n1" === shell n / (4*pi*rad)
n0 = "n0" === shell n / (4*pi*rad**2)
n2x = "n2x" === xshell n
n2y = "n2y" === yshell n
n2z = "n2z" === zshell n
n1x = "n1x" === xshell n / (4*pi*rad)
n1y = "n1y" === yshell n / (4*pi*rad)
n1z = "n1z" === zshell n / (4*pi*rad)
n2px = "n2px" === xshellPrime n
n2py = "n2py" === yshellPrime n
n2pz = "n2pz" === zshellPrime n

sqr_n2v = var "n2vsqr" "{\\left|\\vec{n}_{2v}\\right|^2}" (n2x**2 + n2y**2 + n2z**2)
n1v_dot_n2v = var "n1v_dot_n2v" "{\\vec{n}_{1v}\\cdot\\vec{n}_{2v}}" (n2x*n1x + n2y*n1y + n2z*n1z)

shell_diam, shell, step, xshell, yshell, zshell, shellPrime, xshellPrime, yshellPrime, zshellPrime :: Expression RealSpace -> Expression RealSpace
shell x = ifft ( deltak * fft x)
  where deltak = protect "deltak" "\\delta(k)" $ smear * (4*pi) * rad * sin kR / k
shell_diam x = ifft ( deltak * fft x)
  where deltak = protect "deltak2" "\\delta_{2R}(k)" $ smear * (4*pi) * 2*rad * sin (2*kR) / k
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
shellPrime x = ifft ( delta' * fft x)
  where delta' = protect "deltaprime" "\\delta'(k)" $
                 smear * (4*pi) * ( -sin kR/k - rad * cos kR)
xshellPrime x = ifft ( xdelta' * fft x)
  where xdelta' = protect "xdeltaprime" "\\delta_x'(k)" $
                  smear * (4*pi) * imaginary * kx*( rad * sin kR)/k
yshellPrime x = ifft ( ydelta' * fft x)
  where ydelta' = protect "ydeltaprime" "\\delta_y'(k)" $
                  smear * (4*pi) * imaginary * ky*( rad * sin kR)/k
zshellPrime x = ifft ( zdelta' * fft x)
  where zdelta' = protect "zdeltaprime" "\\delta_z'(k)" $
                  smear * (4*pi) * imaginary * kz*( rad * sin kR)/k
