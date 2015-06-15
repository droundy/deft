module IdealGas
       ( idealgas, idealgas_of_veff, n, kT, of_effective_potential )
       where

import Expression

nQ :: Expression RealSpace
nQ = (mass*kT/2/pi)**1.5
  where mass = var "mH2O" "m_{H_2O}" (18.01528 * gpermol) -- uses molecular weight of water
        gpermol = var "gpermol" "\\frac{\\textrm{g}}{\\textrm{mol}}" 1822.8885

idealgas :: Expression Scalar
idealgas = "Fideal" === integrate (kT*n*(log(n/nQ) - 1))

of_effective_potential :: Expression Scalar -> Expression Scalar
of_effective_potential = substitute (r_var "veff") (r_var "x") .
                         substitute (r_var "x") (exp (- r_var "veff" / kT))

idealgas_of_veff :: Expression Scalar
idealgas_of_veff = "ideal" === integrate (-n_of_veff*(veff  + kT*(1 + log(nQ))))
  where veff = r_var "Veff"
        n_of_veff = exp(-veff/kT)

kT :: Type a => Expression a
kT = s_tex "kT" "kT"

n :: Expression RealSpace
n = "n" === r_var "x"
