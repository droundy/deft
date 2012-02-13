module SomeFunctionals 
       ( fmt, whitebear, wb_contact_at_sphere, idealgas, mu, n,
         phi1, phi2, phi3,
         xshell, yshell, zshell,
         yuwu_zeta, yuwu_contact,
         saft_dispersion, saft_association,
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

yuwu_zeta :: Expression RealSpace
yuwu_zeta = (n2**2 - n2x**2 - n2y**2 - n2z**2)/n2**2

yuwu_contact :: Expression RealSpace
yuwu_contact = n0*yuwu_zeta*ghs
  where ghs = invdiff*(1 + 0.5*(invdiff*zeta2)*yuwu_zeta*(3 + invdiff*zeta2))
        zeta2 = rad*n2/3
        invdiff = 1/(1-n3)

lambda_dispersion, epsilon_dispersion :: Type a => Expression a
lambda_dispersion = s_var "lambda_dispersion"
epsilon_dispersion = s_var "epsilon_dispersion"

length_scaling :: Type a => Expression a
length_scaling = s_var "length_scaling"

eta_for_dispersion, eta_effective :: Expression RealSpace
eta_for_dispersion = r_var "{\\eta_d}" -- ifft (exp (-k**2/2/(length_scaling*lambda_dispersion*rad)**2) * fft n)

-- The following equation is equation 36 in Gil-Villegas 1997 paper.
eta_effective = (c1 + c2*eta_for_dispersion + c3*eta_effective**2)*eta_for_dispersion
  where c1 = 2.25855 - 1.50349*lambda_dispersion + 0.249434*lambda_dispersion**2
        c2 = -0.669270 + 1.40049*lambda_dispersion - 0.827739*lambda_dispersion**2
        c3 = 10.1576 - 15.0427*lambda_dispersion + 5.30827*lambda_dispersion**2

saft_dispersion, saft_association, xsaft, a1, a2 :: Expression RealSpace
saft_dispersion = subst $ n*(a1 + a2/kT)
  where subst = substitute eta_for_dispersion (ifft (exp (-k**2/2/(length_scaling*lambda_dispersion*rad)**2) * fft n))

a1 = -4*(lambda_dispersion**3-1)*epsilon_dispersion* eta_for_dispersion *(1 - eta_effective/2)/(1-eta_effective)**3

a2 = 0.5*khs*eta_for_dispersion* derive (R "{\\eta_d}") 1 a1
     where khs = (1 - eta_for_dispersion)**4/(1 + 4*(eta_for_dispersion + eta_for_dispersion**2))

saft_association = 4*kT*n0*yuwu_zeta*(log xsaft - xsaft/2 + 1/2)

kappa_association, epsilon_association :: Type a => Expression a
kappa_association = s_var "kappa_association"
epsilon_association = s_var "epsilon_association"

xsaft = (sqrt(1 + 8*yuwu_contact*kappa_association*boltz) - 1) / (4*n0*yuwu_contact*kappa_association*boltz)
  where boltz = exp(epsilon_association/kT)-1
