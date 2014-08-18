module SW_liquid
       ( mu, sw_dispersion )
       where

import Expression
import IdealGas ( kT, idealgas )
import FMT ( rad )

sig :: Expression RealSpace

beta :: Type a => Expression a
beta = s_var 1/kT

ksig :: Expression KSpace
ksig = k * sig

ksiglam :: Expression KSpace
ksiglam = k * sig *lambda

sw_dispersion :: Expression Scalar
sw_dispersion = var "SW_dispersion" "F_{\\text{SW}}"

mu :: Type a => Expression a
mu = s_tex "mu" "\\mu"

xi1 :: Expression RealSpace -> Expression RealSpace
xi1 x = ifft (d1k * fft x)
  where d1k = var "xi1k" "\\tilde{\\xi_1}(k)" $
              4*pi*((2 - lambda*ksig**2*(lambda - 1))*cos ksiglam 
                    - ksig*(sin ksig + (1-2*lambda)*sin ksiglam)
                    - cos ksig)/(k*k*k*ksig)
        
d2 :: Expression RealSpace -> Expression RealSpace
d2 x = ifft (d2k * fft x)
  where d2k = var "b2k" "\\tilde{b_2}(k)" $
              (6*sin(ksig) + (ksig**2*(1 - 4*lambda + 3*lambda**2) - 6)*sin(ksiglam) - 2*ksig*cos(ksig) - ksig*(4 + lambda*(ksig**2*(lambda - 1)**2 - 6))*cos(ksiglam))/(k*k*ksig**2)

d3 :: Expression RealSpace -> Expression RealSpace
d3 x = ifft (d3k * fft x)
  where d3k = var "b3k" "\\tilde{b_3}(k)" $
              (ksig*(6*sin(ksig) + (18 -25*lambda + ksig**2*(lambda - 1)**2*(4*lambda - 1))*sin(ksiglam)) + 24*cos(ksig) + (6*ksig**2*(lambda - 1)*(2*lambda - 1) - lambda*ksig**4*(lambda - 1)**3)*cos(ksiglam))/(k*k*ksig**3)

d4 :: Expression RealSpace -> Expression Realspace
d4 x = ifft (d4k * fft x)
  where d4k = var "b4k" "\\tilde{b_4}(k)" $
              ((ksig**2*(lambda - 1)*(36 - 60*lambda + ksig**2*(lambda - 1)**2*(5*lambda - 1)) + 120)*sin(ksiglam) - 120*sin(ksig) + ksig*(24*cos(ksig) - (24*(5*lambda - 4) - 4*ksig**2*(lambda - 1)**2*(5*lambda - 2) + lambda*ksig**4*(lambda - 1)**4))*cos(ksiglam))/(k*k*ksig**4)
