module SomeFunctionals 
       ( whitebear, wb_contact_at_sphere, idealgas, mu, n,
         phi1, phi2, phi3,
         of_effective_potential,
         xshell, yshell, zshell,
         yuwu_zeta, yuwu_correlation,
         eta_for_dispersion, lambda_dispersion, a1, a2, eta_effective,
         saft_dispersion, saft_association, saft_fluid,
         dwbdn3, dwbdn2, dwbdn1, dwbdn2v_over_n2v, dwbdn1v_over_n2v )
       where

import CodeGen

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
shell x = ifft ( smear * (4*pi) * rad * (sin kR / k) * fft x)
step x = ifft ( smear * (4*pi) * (sin kR - kR * cos kR) / k**3 * fft x)
xshell x = ifft ( smear * (4*pi) * imaginary * kx*(rad * cos kR - sin kR/k)/k**2 * fft x)
yshell x = ifft ( smear * (4*pi) * imaginary * ky*(rad * cos kR - sin kR/k)/k**2 * fft x)
zshell x = ifft ( smear * (4*pi) * imaginary * kz*(rad * cos kR - sin kR/k)/k**2 * fft x)

vectorThirdTerm :: Expression RealSpace
vectorThirdTerm = n2*(n2**2 - 3*sqr_n2v)

phi1, phi2, phi3 :: Expression RealSpace
phi1 = var "phi1" "\\Phi_1" $ -n0*log(1-n3)
phi2 = var "phi2" "\\Phi_2" $ (n2*n1 - n1v_dot_n2v)/(1-n3)
phi3 = var "phi3" "\\Phi_3" $ (n3 + (1-n3)**2*log(1-n3))/(36*pi * n3**2 * (1-n3)**2)*vectorThirdTerm


whitebear :: Expression Scalar
whitebear = var "whitebear" "F_{HS}" $ integrate (kT*(phi1+phi2+phi3))

nQ :: Expression RealSpace
nQ = (mass*kT/2/pi)**1.5
  where mass = var "mH2O" "m_{H_2O}" (18.01528 * gpermol) -- uses molecular weight of water
        gpermol = var "gpermol" "\\frac{\\textrm{g}}{\\textrm{mol}}" 1822.8885

idealgas :: Expression Scalar
idealgas = integrate (kT*n*(log(n/nQ) - 1))

mu :: Type a => Expression a
mu = s_tex "mu" "\\mu"

dwbdn3, dwbdn2, dwbdn1, dwbdn1v_over_n2v, dwbdn2v_over_n2v :: Expression RealSpace
dwbdn3 =  d phi1 + d phi2 + d phi3
  where d = derive n3 1

dwbdn2 = d phi1 + d phi2 + d phi3
  where d = derive n2 1

dwbdn1 = derive n1 1 (phi1 + phi2 + phi3)

dwbdn1v_over_n2v = derive n1v_dot_n2v 1 (phi1 + phi2 + phi3)

dwbdn2v_over_n2v = 2*derive sqr_n2v 1 (phi1 + phi2 + phi3)

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
yuwu_zeta = var "zeta_yuwu" "{\\zeta}" $ (n2**2 - sqr_n2v)/n2**2

yuwu_correlation :: Expression RealSpace
yuwu_correlation = var "ghsyuwu" "g_{HS}^{\\textit{YuWu}}" ghs
    where ghs = 1/(1-n3) + rad/2*n2*yuwu_zeta/(1-n3)**2 +
                rad**2/18*n2**2*yuwu_zeta/(1-n3)**3

lambda_dispersion, epsilon_dispersion :: Type a => Expression a
lambda_dispersion = s_tex "lambda_dispersion" "\\lambda_d"
epsilon_dispersion = s_tex "epsilon_dispersion" "\\epsilon_d"

eta_for_dispersion, eta_effective :: Expression RealSpace
--eta_for_dispersion = var "eta_d" "{\\eta_d}" $ (4*pi*rad**3/3)*n
eta_for_dispersion = var "eta_d" "{\\eta_d}" $
                     ((4*pi*rad**3/3)*ifft (exp (-0.5*k**2*(2*length_scaling*lambda_dispersion*rad)**2) * fft n))
                       where length_scaling = s_tex "length_scaling" "l_s"

-- The following equation is equation 36 in Gil-Villegas 1997 paper.
eta_effective = var "eta_eff" "\\eta_{\\textit{eff}}" $
                (c1 + c2*eta_for_dispersion + c3*eta_for_dispersion**2)*eta_for_dispersion
  where c1 = "c1" === 2.25855 - 1.50349*lambda_dispersion + 0.249434*lambda_dispersion**2
        c2 = "c2" === -0.669270 + 1.40049*lambda_dispersion - 0.827739*lambda_dispersion**2
        c3 = "c3" === 10.1576 - 15.0427*lambda_dispersion + 5.30827*lambda_dispersion**2

{-
deta_effective_by_deta_for_dispersion :: Expression RealSpace
deta_effective_by_deta_for_dispersion = (c1 + 2*c2*eta_for_dispersion + 3*c3*eta_for_dispersion**2)
  where c1 = "c1" === 2.25855 - 1.50349*lambda_dispersion + 0.249434*lambda_dispersion**2
        c2 = "c2" === -0.669270 + 1.40049*lambda_dispersion - 0.827739*lambda_dispersion**2
        c3 = "c3" === 10.1576 - 15.0427*lambda_dispersion + 5.30827*lambda_dispersion**2

deta_effective_by_dlambda :: Expression RealSpace
deta_effective_by_dlambda = (c1 + 2*c2*eta_for_dispersion + 3*c3*eta_for_dispersion**2)*eta_for_dispersion
  where c1 = -1.50349 + 2*0.249434*lambda_dispersion
        c2 = 1.40049 - 2*0.827739
        c3 = -15.0427 + 2*5.30827*lambda_dispersion
-}

saft_dispersion, saft_association :: Expression Scalar
saft_dispersion = "Fdisp" === (("a1integrated" === integrate (n*a1)) + ("a2integrated" === integrate (n*a2/kT)))

a1 = "a1" === -4*(lambda_dispersion**3-1)*epsilon_dispersion*eta_for_dispersion*ghs
  where ghs = var "ghs" "g_{HS}" $ (1 - eta_effective/2)/(1-eta_effective)**3

da1_dlambda, da1_detad :: Expression RealSpace
da1_dlambda = substitute (r_var "etad") eta_for_dispersion $ derive (lambda_dispersion :: Expression RealSpace) (1 :: Expression RealSpace) $ substitute eta_for_dispersion (r_var "etad") a1

da1_detad = derive eta_for_dispersion 1 a1

a2 = "a2" === 0.5*khs*epsilon_dispersion*eta_for_dispersion*da1_detad
     where khs = var "KHS" "{\\kappa_{HS}}" $ (1 - eta_for_dispersion)**4/(1 + 4*(eta_for_dispersion + eta_for_dispersion**2))

saft_association = "Fassoc" === integrate (4*kT*n0*yuwu_zeta*(log xsaft - xsaft/2 + 1/2))

kappa_association, epsilon_association :: Type a => Expression a
kappa_association = s_tex "kappa_association" "\\kappa_a"
epsilon_association = s_tex "epsilon_association" "\\epsilon_a"

xsaft, deltasaft, a1, a2 :: Expression RealSpace
xsaft = "X" === (sqrt(1 + 8*n0*yuwu_zeta*deltasaft) - 1) / (4*n0*yuwu_zeta*deltasaft)

deltasaft = var "deltasaft" "{\\Delta}" $ gSW*kappa_association*boltz
  where boltz = "boltz" === exp(epsilon_association/kT)-1

gSW :: Expression RealSpace
gSW = "gSW" ===
      yuwu_correlation
      + (1/4/kT)*(da1_detad - lambda_dispersion/(3*eta_for_dispersion)*da1_dlambda)

saft_fluid :: Expression Scalar
saft_fluid = "FSAFT" === (idealgas + whitebear + saft_association + saft_dispersion + integrate (n*mu))

of_effective_potential :: Expression Scalar -> Expression Scalar
of_effective_potential = substitute (r_var "veff") (r_var "x") .
                         substitute (r_var "x") (exp (- r_var "veff" / kT))
