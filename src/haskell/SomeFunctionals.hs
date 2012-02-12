module SomeFunctionals 
       ( fmt, whitebear, wb_contact_at_sphere, idealgas, mu, n,
         phi1, phi2, phi3,
         xshell, yshell, zshell,
         dwbdn3, dwbdn2, dwbdn1, dwbdn2v_over_n2v, dwbdn1v_over_n2v )
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

n, n3, n2, n1, n0, n2x, n2y, n2z, n1x, n1y, n1z, n1v_dot_n2v, sqr_n2v :: Expression RealSpace
n = r_var "x"
n3 = r_var "n3" -- step n
n2 = r_var "n2" -- shell n
n1 = r_var "n1" -- n2 / (4*pi*rad)
n0 = r_var "n0" -- n2 / (4*pi*rad**2)
n2x = r_var "n2x" -- xshell n
n2y = r_var "n2y" -- yshell n
n2z = r_var "n2z" -- zshell n
n1x = r_var "n1x"
n1y = r_var "n1y"
n1z = r_var "n1z"
n1v_dot_n2v = r_var n1vdotname
sqr_n2v = r_var n2vsqrname

n1vdotname, n2vsqrname :: String
n1vdotname = "{\\vec{n}_{1v}\\cdot\\vec{n}_{2v}}"
n2vsqrname = "{\\left|\\vec{n}_{2v}\\right|}"

fmt :: Expression RealSpace -> Expression RealSpace
fmt = substitute n3 (step n) .
      substitute n2 (shell n) .
      substitute n1 (shell n/(4*pi*rad)) .
      substitute n0 (shell n/(4*pi*rad**2)) .
      substitute n2x (xshell n) .
      substitute n2y (yshell n) .
      substitute n2z (zshell n) .
      substitute n1x (xshell n/(4*pi*rad)) .
      substitute n1y (yshell n/(4*pi*rad)) .
      substitute n1z (zshell n/(4*pi*rad)) .
      substitute n1v_dot_n2v (n1x*n2x + n1y*n2y + n1z*n2z) .
      substitute sqr_n2v (n2x**2 + n2y**2 + n2z**2)

shell, step, xshell, yshell, zshell :: Expression RealSpace -> Expression RealSpace
shell x = ifft ( smear * (4*pi) * rad * (sin kR / k) * fft x)
step x = ifft ( smear * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft x)
xshell x = ifft ( smear * (4*pi) * i * kx*(rad * cos kR - sin kR/k)/k**2 * fft x)
yshell x = ifft ( smear * (4*pi) * i * ky*(rad * cos kR - sin kR/k)/k**2 * fft x)
zshell x = ifft ( smear * (4*pi) * i * kz*(rad * cos kR - sin kR/k)/k**2 * fft x)

vectorThirdTerm :: Expression RealSpace
vectorThirdTerm = n2*(n2**2 - 3*sqr_n2v)

phi1, phi2, phi3 :: Expression RealSpace
phi1 = -n0*log(1-n3)
phi2 = (n2*n1 - n1v_dot_n2v)/(1-n3)
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

dwbdn3, dwbdn2, dwbdn1, dwbdn1v_over_n2v, dwbdn2v_over_n2v :: Expression RealSpace
dwbdn3 =  d phi1 + d phi2 + d phi3
  where d = derive (R "n3") 1

dwbdn2 = d phi1 + d phi2 + d phi3
  where d = derive (R "n2") 1

dwbdn1 = derive (R "n1") 1 (phi1 + phi2 + phi3)

dwbdn1v_over_n2v = derive (R n1vdotname) 1 (phi1 + phi2 + phi3)

dwbdn2v_over_n2v = 2*derive (R  n2vsqrname) 1 (phi1 + phi2 + phi3)

--dwbdn2v, dwbdn1v :: Expression RealSpace
--dwbdn2v = 1/(4*pi*rad)/(1-n3) - 6*n2*(n3 + (1-n3)**2*log(1-n3)/(36*pi*n3**2*(1-n3)**2))

--dwbdn1v = 1/(1-n3)

wb_contact_at_sphere :: Expression RealSpace
wb_contact_at_sphere =
  1/(4*4*pi*rad**2)*(shell (dwbdn3 + 2*dwbdn2/rad + dwbdn1/(4*pi*rad**2)) +
                     xshell (n2x * vecpart) + yshell (n2y * vecpart) + zshell (n2z * vecpart))
    where vecpart = 0*dwbdn1v_over_n2v / (4*pi*rad**2) + 2/rad*dwbdn2v_over_n2v
  -- + xshell (n2x*(dwbdn1v/(4*pi*rad**2) + dwbdn2v*2/rad))
  -- + yshell (n2y*(dwbdn1v/(4*pi*rad**2) + dwbdn2v*2/rad))
  -- + zshell (n2z*(dwbdn1v/(4*pi*rad**2) + dwbdn2v*2/rad))
