{- SFMT defines the fundamental measures that are the raw input of
soft fundamental measure theory. -}

module SFMT
       ( sfmt, sfmt_fluid_n, sfmt_fluid_Veff, homogeneous_sfmt_fluid,
         phi1, phi2, phi3,
         n0, n1, n2, n3, n2v, n1v, sqr_n2v, n1v_dot_n2v )
       where

import Expression
import IdealGas ( kT, idealgas )

phi1, phi2, phi3 :: Expression Scalar
phi1 = var "phi1" "\\Phi_1" $ integrate $ -kT*n0*log(1-n3)
phi2 = var "phi2" "\\Phi_2" $ integrate $ kT*(n2*n1 - n1v_dot_n2v)/(1-n3)
phi3 = var "phi3" "\\Phi_3" $ integrate $ kT*(n2**3/3 - sqr_n2v*n2)/(8*pi*(1-n3)**2)

sfmt :: Expression Scalar
sfmt = var "sfmt" "F_{\\text{soft}}" $ (phi1 + phi2 + phi3)

sfmt_fluid_n :: Expression Scalar
sfmt_fluid_n = substitute n (r_var "n") $
               sfmt + idealgas + integrate (n * (r_var "Vext" - s_var "mu"))

sfmt_fluid_Veff :: Expression Scalar
sfmt_fluid_Veff = substitute n ("n" === exp(- r_var "Veff"/kT)) $
                  sfmt + idealgas + integrate (n * (r_var "Vext" - s_var "mu"))

homogeneous_sfmt_fluid :: Expression Scalar
homogeneous_sfmt_fluid = makeHomogeneous $ substitute n (r_var "n") $
                         sfmt + idealgas - integrate (n * s_var "mu")

sigma :: Type a => Expression a
sigma = s_var "sigma"

betaeps :: Type a => Expression a
betaeps = s_var "epsilon"/kT

xi :: Type a => Expression a
xi = Scalar $ var "Xi" "\\Xi" (alpha/(sqrt(pi)*(sqrt(betaeps * log 2 ) + 6*log 2)))

alpha :: Type a => Expression a
alpha = Scalar $ var "alpha" "\\alpha" (sigma*(2/(1+6*sqrt(log 2 / betaeps)))**(1.0/6.0))

n, n3, n2, n1, n0, n1v_dot_n2v, sqr_n2v :: Expression RealSpace
n = "n" === r_var "x"
n3 = "n3" === w3 n
n2 = "n2" === w2 n
n1 = "n1" === w1 n
n0 = "n0" === 2/alpha*(2*n1 - n2/(2*pi*alpha))

n2v, n1v :: Vector RealSpace
n2v = "n2v" `nameVector` w2v n
n1v = "n1v" `nameVector` w1v n

sqr_n2v = var "n2vsqr" "{\\left|\\vec{n}_{2v}\\right|^2}" (n2v `dot` n2v)
n1v_dot_n2v = var "n1v_dot_n2v" "{\\vec{n}_{1v}\\cdot\\vec{n}_{2v}}" (n2v `dot` n1v)


w3 :: Expression RealSpace -> Expression RealSpace
w3 x = ifft ( w3k * fft x)
  where w3k = var "w3k" "\\tilde{w_3}(k)" $
              4*pi/(k**2)*exp(-(xi*k/2)**2)*((1+xi**2*k**2/2)*sin kalphao2/k - alpha/2*cos kalphao2)
        kalphao2 = k*alpha/2

w1 :: Expression RealSpace -> Expression RealSpace
w1 x = ifft ( w1k * fft x)
  where w1k = var "w1k" "\\tilde{w_1}(k)" $
              1/k*exp(-(xi*k/2)**2)*sin(k*alpha/2)

w2 :: Expression RealSpace -> Expression RealSpace
w2 x = ifft ( w2k * fft x)
  where w2k = var "w2k" "\\tilde{w_2}(k)" $
              2*pi*exp(-(xi*k/2)**2)*(xi**2*cos(k*alpha/2) + alpha*sin(k*alpha/2)/k)

w2v :: Expression RealSpace -> Vector RealSpace
w2v x = vector_convolve w2vk x
  where w2vk = kvec *. (imaginary * pi*exp(-(xi*k/2)**2)*((alpha**2-xi**4*k**2)*cos(k*alpha/2)-
                                               2*alpha*(xi**2*k+1/k)*sin(k*alpha/2))/k**2)

w1v :: Expression RealSpace -> Vector RealSpace
w1v x = vector_convolve w1vk x
  where w1vk = kvec *. (imaginary * exp(-(xi*k/2)**2)*(alpha/2*cos(k*alpha/2)-(xi**2*k/2+1/k)*sin(k*alpha/2))/k**2)
