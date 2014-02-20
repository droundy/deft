{- SFMT defines the fundamental measures that are the raw input of
soft fundamental measure theory. -}

module SFMT
       ( sfmt, n0, n1, n2, n3, kR, n2v, n1v, sqr_n2v, n1v_dot_n2v )
       where

import Expression
import WhiteBear ( kT )
import FMT ( rad )

phi1, phi2, phi3 :: Expression Scalar
phi1 = var "phi1" "\\Phi_1" $ integrate $ -kT*n0*log(1-n3)
phi2 = var "phi2" "\\Phi_2" $ integrate $ kT*(n2*n1 - n1v_dot_n2v)/(1-n3)
phi3 = var "phi3" "\\Phi_3" $ integrate $ kT*(n2**3/3 - sqr_n2v*n2)/(8*pi*(1-n3)**2)

sfmt :: Expression Scalar
sfmt = var "sfmt" "F_{\\text{soft}}" $ (phi1 + phi2 + phi3)

betaV0 :: Type a => Expression a
betaV0 = s_var "V0"/kT

kR :: Expression KSpace
kR = k * rad

a :: Type a => Expression a
a = "a" === 2*rad/sqrt(betaV0*pi*log 2)

sigma :: Type a => Expression a
sigma = var "sigma" "\\sigma" (2*rad*(1 - sqrt(log 2/betaV0)))

n, n3, n2, n1, n0, n1v_dot_n2v, sqr_n2v :: Expression RealSpace
n = "n" === r_var "x"
n3 = "n3" === w3 n
n2 = "n2" === w2 n
n1 = "n1" === w1 n
n0 = "n0" === 2/sigma*(2*n1 - n2/(2*pi*sigma))

n2v, n1v :: Vector RealSpace
n2v = "n2v" `nameVector` w2v n
n1v = "n1v" `nameVector` w1v n

sqr_n2v = var "n2vsqr" "{\\left|\\vec{n}_{2v}\\right|^2}" (n2v `dot` n2v)
n1v_dot_n2v = var "n1v_dot_n2v" "{\\vec{n}_{1v}\\cdot\\vec{n}_{2v}}" (n2v `dot` n1v)


w3 :: Expression RealSpace -> Expression RealSpace
w3 x = ifft ( w3k * fft x)
  where w3k = var "w3k" "\\tilde{w_3}(k)" $
              4*pi/(k**2)*exp(-(a*k/2)**2)*((1+a**2*k**2/2)*sin ksigmao2/k - sigma/2*cos ksigmao2)
        ksigmao2 = k*sigma/2

w1 :: Expression RealSpace -> Expression RealSpace
w1 x = ifft ( w1k * fft x)
  where w1k = var "w1k" "\\tilde{w_1}(k)" $
              1/k*exp(-(a*k/2)**2)*sin(k*sigma/2)

w2 :: Expression RealSpace -> Expression RealSpace
w2 x = ifft ( w2k * fft x)
  where w2k = var "w2k" "\\tilde{w_2}(k)" $
              2*pi*exp(-(a*k/2)**2)*(a**2*cos(k*sigma/2) + sigma*sin(k*sigma/2)/k)

w2v :: Expression RealSpace -> Vector RealSpace
w2v x = vifft ( w2vk *. fft x)
  where w2vk = kvec *. (-2*pi*exp(-(a*k/2)**2)*((a**4*k**2+sigma**2)/2*cos(k*sigma/2)+
                                               (a**2*k+1/k)*sigma*sin(k*sigma/2))/k**2)

w1v :: Expression RealSpace -> Vector RealSpace
w1v x = vifft ( w1vk *. fft x)
  where w1vk = kvec *. (exp(-(a*k/2)**2)*(sigma/2*cos(k*sigma/2)-(a**2*k/2+1/k)*sin(k*sigma/2))/k**2)
