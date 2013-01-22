{- SFMT defines the fundamental measures that are the raw input of
soft fundamental measure theory. -}

module SFMT
       ( n0, n1, n2, n3 )
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
n3 = "n3" === w1 n --wrong, will change later
n2 = "n2" === w1 n --wrong, will change later
n1 = "n1" === w1 n
n0 = "n0" === w1 n --wrong, will change later

n2v, n1v :: Vector RealSpace
n2v = "n2v" `nameVector` w1v n --wrong, will change later
n1v = "n1v" `nameVector` w1v n

sqr_n2v = var "n2vsqr" "{\\left|\\vec{n}_{2v}\\right|^2}" (n2v `dot` n2v)
n1v_dot_n2v = var "n1v_dot_n2v" "{\\vec{n}_{1v}\\cdot\\vec{n}_{2v}}" (n2v `dot` n1v)

w1 :: Expression RealSpace -> Expression RealSpace
w1 x = ifft ( w1k * fft x)
  where w1k = protect "w1k" "w_1(k)" $
              (2*betaV0/rad)**(1/2)/((2*betaV0/rad)**2 + k**2) *
                     ((2*betaV0/rad/k*sin(k*rad) -cos(k*rad) + 1)*exp(-betaV0/2))

w1v :: Expression RealSpace -> Vector RealSpace
w1v x = vifft $ deltav *. fft x
  where deltav = vprotect "w1vk" (\i -> "\\mathbf{w}_{1V}(k)\\cdot \\mathbf{\\hat{"++i++"}}") $
                 (2*betaV0/rad)**(1/2)/((2*betaV0/rad)**2 + k**2) *
                 ((2*betaV0/rad/k*sin(k*rad) -cos(k*rad) + 1)*erf(-betaV0/2)) .* kvec
