{- SFMT defines the fundamental measures that are the raw input of
soft fundamental measure theory. -}

module SFMT
       ( sfmt, n0, n1, n2, n3, kR, n2v, n1v, sqr_n2v, n1v_dot_n2v, w1v )
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
betaV0 = s_var "V0"/s_var "kT"

kR :: Expression KSpace
kR = k * rad

n, n3, n2, n1, n0, n1v_dot_n2v, sqr_n2v :: Expression RealSpace
n = "n" === r_var "x"
n3 = "n3" === w3 n
n2 = "n2" === w2 n
n1 = "n1" === w1 n
n0 = "n0" === w0 n
{-
--These are temporary names for the soft weighted densities
soft_n, soft_n3, soft_n2, soft_n1, soft_n0 :: Expression RealSpace
soft_n "soft_n" === r_var "x"
soft_n3 "soft_n3" === soft_w3 soft_n
soft_n2 "soft_n2" === soft_w2 soft_n
soft_n1 "soft_n1" === soft_w1 soft_n
soft_n0 "soft_n0" === soft_w0 soft_n
-}
n2v, n1v :: Vector RealSpace
n2v = "n2v" `nameVector` w2v n
n1v = "n1v" `nameVector` w1v n

sqr_n2v = var "n2vsqr" "{\\left|\\vec{n}_{2v}\\right|^2}" (n2v `dot` n2v)
n1v_dot_n2v = var "n1v_dot_n2v" "{\\vec{n}_{1v}\\cdot\\vec{n}_{2v}}" (n2v `dot` n1v)

gamma, b :: Type a => Expression a
b = 2*gamma/(sqrt(pi*gamma)-1)/rad**2
gamma = "gamma" === 2*((sqrt(pi*betaV0)+sqrt(pi*betaV0-16*sqrt(betaV0)))/8)**2

mydr :: Expression Scalar
mydr = 0.001

mydk :: Double
mydk = 0.1

mykmax :: Double
mykmax = 1000

mys :: Symmetry
mys = Spherical { dk = mydk, kmax = mykmax, rresolution = mydr*rad, rmax = rad }

myvs :: Symmetry
myvs = VectorS { dk = mydk, kmax = mykmax, rresolution = mydr*rad, rmax = rad }

r :: Expression Scalar
r = s_var "r"

w3 :: Expression RealSpace -> Expression RealSpace
w3 x = ifft ( w3k * fft x)
  where w3k = transform mys $ -1/(sqrt(pi*gammax)-1)*(
            1 - exp(-gammax*(1-r/rad)**2) - sqrt(pi*gammax)*erf(sqrt gammax*(1-r/rad))
          )
        gammax = 360.38

w0 :: Expression RealSpace -> Expression RealSpace
w0 x = ifft ( w0k * fft x)
  where w0k = transform mys $ b*exp(-gamma*(1-r/rad)**2)*heaviside(rad - r)/(4*pi*r)

w1 :: Expression RealSpace -> Expression RealSpace
w1 x = ifft ( w1k * fft x)
  where w1k = transform mys $ b*exp(-gamma*(1-r/rad)**2)*heaviside(rad - r)/(4*pi)

w2 :: Expression RealSpace -> Expression RealSpace
w2 x = ifft ( w2k * fft x)
  where w2k = transform mys $ b*r*exp(-gamma*(1-r/rad)**2)*heaviside(rad - r)

w2v :: Expression RealSpace -> Vector RealSpace
w2v x = vifft ( w2vk *. fft x)
  where w2vk = kvec *. transform myvs (b*r*exp(-gamma*(1-r/rad)**2)*heaviside(rad - r)/r)

w1v :: Expression RealSpace -> Vector RealSpace
w1v x = vifft ( w1vk *. fft x)
  where w1vk = kvec *. transform myvs (b*r*exp(-gamma*(1-r/rad)**2)*heaviside(rad - r)/(4*pi*r**2))
{-
soft_w2 ::Expression RealSpace -> Expression RealSpace
soft_w2 x = ifft ( stepk * fft x)
  where stepk = protect "step" "\\Theta(k)" $ smear * (4*pi) * (sin kR - kR * cos kR) / k**3

soft_w1 ::Expression RealSpace -> Expression RealSpace
soft_w1 x = ifft ( stepk * fft x)
  where stepk = protect "step" "\\Theta(k)" $ smear * (4*pi) * (sin kR - kR * cos kR) / k**3

soft_w0 ::Expression RealSpace -> Expression RealSpace
soft_w0 x = ifft ( stepk * fft x)
  where stepk = protect "step" "\\Theta(k)" $ smear * (4*pi) * (sin kR - kR * cos kR) / k**3
-}
