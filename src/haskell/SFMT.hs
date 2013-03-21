{- SFMT defines the fundamental measures that are the raw input of
soft fundamental measure theory. -}

module SFMT
       ( n0, n1, n2, n3, kR, n2v, n1v, sqr_n2v, n1v_dot_n2v, w1v )
       where

import Expression

rad :: Type a => Expression a
rad = s_var "R"

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

n2v, n1v :: Vector RealSpace
n2v = "n2v" `nameVector` w2v n
n1v = "n1v" `nameVector` w1v n

sqr_n2v = var "n2vsqr" "{\\left|\\vec{n}_{2v}\\right|^2}" (n2v `dot` n2v)
n1v_dot_n2v = var "n1v_dot_n2v" "{\\vec{n}_{1v}\\cdot\\vec{n}_{2v}}" (n2v `dot` n1v)

gamma, b :: Type a => Expression a
b = 2*gamma/(sqrt(pi*gamma)-1)/rad**2
gamma = "gamma" === 2*((sqrt(pi*betaV0)+sqrt(pi*betaV0-16*sqrt(betaV0)))/8)**2

w3 :: Expression RealSpace -> Expression RealSpace
w3 x = ifft ( w3k * fft x)
  where w3r = "w3r" === b*rad**2/(2*gamma)*(
          1 - exp(-gamma*(1-rmag/rad)**2) - sqrt(pi*gamma)*erf(sqrt gamma*(1-rmag/rad))
          )*heaviside(rad - rmag)
        w3k = "w3k" === fft w3r

w0 :: Expression RealSpace -> Expression RealSpace
w0 x = ifft ( w0k * fft x)
  where w0r = "w0r" === b*exp(-gamma*(1-rmag/rad)**2)*heaviside(rad - rmag)/(4*pi*rmag)
        w0k = "w0k" === fft w0r

w1 :: Expression RealSpace -> Expression RealSpace
w1 x = ifft ( w1k * fft x)
  where w1r = "w1r" === b*exp(-gamma*(1-rmag/rad)**2)*heaviside(rad - rmag)/(4*pi)
        w1k = "w1k" === fft w1r

w2 :: Expression RealSpace -> Expression RealSpace
w2 x = ifft ( w2k * fft x)
  where w2r = "w2r" === b*rmag*exp(-gamma*(1-rmag/rad)**2)*heaviside(rad - rmag)
        w2k = "w2k" === fft w2r

w2v :: Expression RealSpace -> Vector RealSpace
w2v x = vifft ( w2vk *. fft x)
  where w2vr = "w2vr" `nameVector`
               (b*rmag*exp(-gamma*(1-rmag/rad)**2)*heaviside(rad - rmag)/rmag).*rvec
        w2vk = "w2vk" `nameVector` vfft w2vr

w1v :: Expression RealSpace -> Vector RealSpace
w1v x = vifft ( w1vk *. fft x)
  where w1vr = "w1vr" `nameVector`
               (b*rmag*exp(-gamma*(1-rmag/rad)**2)*heaviside(rad - rmag)/(4*pi*rmag**2)).*rvec
        w1vk = "w1vk" `nameVector` vfft w1vr
