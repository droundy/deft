module SomeFunctionals 
       ( whitebear, idealgas, mu, n )
       where

import CodeGen

kT :: Type a => Expression a
kT = s_var "kT"

rad :: Type a => Expression a
rad = s_var "R"

kdr :: Expression KSpace
kdr = k * s_var "dr"

kR :: Expression KSpace
kR = k * rad

i :: Type a => Expression a
i = s_var "complex(0,1)"

smear :: Expression KSpace
smear = exp (-6.0*kdr*kdr)

n, n3, n2, n2x, n2y, n2z :: Expression RealSpace
n = r_var "x"
n3 = ifft ( smear * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft n)
n2 = ifft ( smear * (4*pi) * rad * (sin kR / k) * fft n)
n2x = ifft ( smear * (4*pi) * i * kx*(rad * cos kR - sin kR/k)/k**2 * fft n)
n2y = ifft ( smear * (4*pi) * i * ky*(rad * cos kR - sin kR/k)/k**2 * fft n)
n2z = ifft ( smear * (4*pi) * i * kz*(rad * cos kR - sin kR/k)/k**2 * fft n)

vectorThirdTerm :: Expression RealSpace
vectorThirdTerm = n2*(n2**2 - 3*(n2x**2 + n2y**2 + n2z**2))

phi1, phi2, phi3 :: Expression RealSpace
phi1 = (-1/(4*pi*rad**2))*n2*log(1-n3)
phi2 = (n2**2 - n2x**2 - n2y**2 - n2z**2)/(4*pi*rad*(1-n3))
phi3 = (n3 + (1-n3)**2*log(1-n3))/(36*pi * n3**2 * (1-n3)**2)*vectorThirdTerm

whitebear :: Expression RealSpace
whitebear = kT * (phi1+phi2+phi3)

nQ :: Expression RealSpace
nQ = (mass*kT/2/pi)**1.5
  where mass = 18.01528*1822.8885 -- uses molecular weight of water

idealgas :: Expression RealSpace
idealgas = kT*n*(log(n/nQ) - 1)

mu :: Type a => Expression a
mu = s_var "mu"
