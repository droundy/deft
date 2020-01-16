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
-- The following was Rosenfelds early vector version of the functional
-- phi3 = var "phi3" "\\Phi_3" $ integrate $ kT*(n2**3/3 - sqr_n2v*n2)/(8*pi*(1-n3)**2)

-- This is the fixed version, which comes from dimensional crossover
-- phi3 = var "phi3" "\\Phi_3" $ integrate $ kT*n2**3*(1 - sqr_n2v/n2**2)**3/(24*pi*(1-n3)**2) + n2xx

-- This is the tensor version
phi3 = var "phi3" "\\Phi_3" $ integrate $ "phi3_density" ===
  kT*(n2**3 - 3*sqr_n2v*n2 + 9/2*(tensor_vector - trace_tensor_cubed))/(24*pi*(1-n3)**2)

sfmt :: Expression Scalar
sfmt = var "sfmt" "F_{\\text{soft}}" $ (phi1 + phi2 + phi3)

sfmt_fluid_n :: Expression Scalar
sfmt_fluid_n = substitute n (r_var "n") $
               sfmt + idealgas + ("external" === integrate (n * (r_var "Vext" - s_var "mu")))

sfmt_fluid_Veff :: Expression Scalar
sfmt_fluid_Veff = substitute n ("n" === exp(- r_var "Veff"/kT)) $
                  sfmt + idealgas + ("external" === integrate (n * (r_var "Vext" - s_var "mu")))

homogeneous_sfmt_fluid :: Expression Scalar
homogeneous_sfmt_fluid = makeHomogeneous $ substitute n (r_var "n") $
                         sfmt + idealgas - integrate (n * s_var "mu")

sigma :: Type a => Expression a
sigma = s_var "sigma"

betaeps :: Type a => Expression a
betaeps = s_var "epsilon"/kT

xi :: Type a => Expression a
xi = s_tex "Xi" "\\Xi" -- scalar $ var "Xi" "\\Xi" (alpha/(6*sqrt(pi)*(sqrt(betaeps * log 2 ) + log 2)))

alpha :: Type a => Expression a
alpha = scalar $ var "alpha" "\\alpha" (sigma*(2/(1+sqrt(log 2 / betaeps)))**(1.0/6.0))

n, n3, n2, n1, n0, n1v_dot_n2v, sqr_n2v :: Expression RealSpace
n = "n" === r_var "x"
n3 = "n3" === w3 n
n2 = "n2" === w2 n
n1 = "n1" === w1 n
n0 = "n0" === w0 n

n2xx = "n2xx" === ifft ( w2transverse * (1 - 3*(kx**2)/k**2) * fft n)
n2yy = "n2yy" === ifft ( w2transverse * (1 - 3*(ky**2)/k**2) * fft n)
n2zz = "n2zz" === ifft ( w2transverse * (1 - 3*(kz**2)/k**2) * fft n)
n2xy = "n2xy" === ifft ( w2transverse * (3*kx*ky/k**2) * fft n)
n2yz = "n2yz" === ifft ( w2transverse * (3*ky*kz/k**2) * fft n)
n2zx = "n2zx" === ifft ( w2transverse * (3*kz*kx/k**2) * fft n)

tensor_vector = var "tensor_vector" "n_{2v}\\cdot n_{2m}\\cdot n_{2v}"
    (n2xx*n2x*n2x + n2yy*n2y*n2y + n2zz*n2z*n2z +
      2*n2xy*n2x*n2y + 2*n2yz*n2y*n2z + 2*n2zx*n2z*n2x)
  where n2x = xhat `dot` n2v
        n2y = yhat `dot` n2v
        n2z = zhat `dot` n2v

trace_tensor_cubed = var "trace_tensor_cubed" "\\mathrm{Tr}(n_{2m}^3)"
     (n2xx*n2xx*n2xx +   -- xx xx xx
       n2yy*n2yy*n2yy +   -- yy yy yy
       n2zz*n2zz*n2zz +   -- zz zz zz
       3*n2xx*n2xy*n2xy + -- xx xy yx three places to put xx
       3*n2xx*n2zx*n2zx + -- xx xz zx three places to put xx
       3*n2xy*n2yy*n2xy + -- xy yy yx three places to put yy
       3*n2yy*n2yz*n2yz + -- yy yz zy three places to put yy
       3*n2zx*n2zz*n2zx + -- xz zz zx three places to put zz
       3*n2yz*n2zz*n2yz + -- yz zz zy three places to put zz
       3*n2xy*n2yz*n2zx + -- xy yz zx three cyclic permutations
       3*n2zx*n2yz*n2xy)  -- xz zy yx == the above, but these are the non-cyclic permutations

n2v, n1v :: Vector RealSpace
n2v = "n2v" `nameVector` w2v n
n1v = "n1v" `nameVector` w1v n

sqr_n2v = var "n2vsqr" "{\\left|\\vec{n}_{2v}\\right|^2}" (n2v `dot` n2v)
n1v_dot_n2v = var "n1v_dot_n2v" "{\\vec{n}_{1v}\\cdot\\vec{n}_{2v}}" (n2v `dot` n1v)


w3 :: Expression RealSpace -> Expression RealSpace
w3 x = ifft ( w3k * fft x)
  where w3k = var "w3k" "\\tilde{w_3}(k)" $
              4*pi/(k**2)*exp(-((xi/sqrt 2)*k/2)**2)*((1+(xi/sqrt 2)**2*k**2/2)*sin kalphao2/k - alpha/2*cos kalphao2)
        kalphao2 = k*alpha/2

w1 :: Expression RealSpace -> Expression RealSpace
w1 x = ifft ( w1k * fft x)
  where w1k = var "w1k" "\\tilde{w_1}(k)" $
              1/k*exp(-((xi/sqrt 2)*k/2)**2)*sin(k*alpha/2)

w0 :: Expression RealSpace -> Expression RealSpace
w0 x = ifft ( w0k * fft x)
  where w0k = var "w0k" "\\tilde{w_0}(k)" $
              (sqrt (pi/2)/(k*xi)) *exp(-0.5*alpha**2/xi**2)*(real_part_complex_erf(k*xi/2/sqrt 2 + imaginary*alpha/sqrt 2/xi))

w2 :: Expression RealSpace -> Expression RealSpace
w2 x = ifft ( w2k * fft x)
  where w2k = var "w2k" "\\tilde{w_2}(k)" $
              2*pi*exp(-((xi/sqrt 2)*k/2)**2)*((xi/sqrt 2)**2*cos(k*alpha/2) + alpha*sin(k*alpha/2)/k)

w2v :: Expression RealSpace -> Vector RealSpace
w2v x = vector_convolve w2vk x
  where w2vk = kvec *. (imaginary * pi*exp(-((xi/sqrt 2)*k/2)**2)*((alpha**2-(xi/sqrt 2)**4*k**2)*cos(k*alpha/2)-
                                               2*alpha*((xi/sqrt 2)**2*k+1/k)*sin(k*alpha/2))/k**2)

w1v :: Expression RealSpace -> Vector RealSpace
w1v x = vector_convolve w1vk x
  where w1vk = kvec *. (imaginary * exp(-((xi/sqrt 2)*k/2)**2)*(alpha/2*cos(k*alpha/2)-((xi/sqrt 2)**2*k/2+1/k)*sin(k*alpha/2))/k**2)

w2transverse :: Expression KSpace
w2transverse = -- var "w2transversek" "\\tilde{w_{2\\perp}}(k)" $
               --(4*pi*exp(-(xi*k)**2/8)*(((xi*k)/(3*sqrt 2) - 1)*cos(k*alpha/2) + k**2*xi*alpha/(12*sqrt 2)*sin(k*alpha/2))/k**2
               --      +
               --     (2*pi)**(3/2)/((k**3)*xi)*exp(-alpha**2/(2*(xi**2)))*
               --     2*real_part_complex_erf(k*xi/2**1.5 + imaginary*alpha/(xi*sqrt 2))
               (-(pi*xi**2)/3 - 4*pi/k**2)*exp(-(xi*k)**2/8)*cos(k*alpha/2) - 2*pi*alpha/(3*k)*exp(-(xi*k)**2/8)*sin(k*alpha/2)
                     +
                    (2*pi)**(3/2)/((k**3)*xi)*exp(-alpha**2/(2*(xi**2)))*
                    2*real_part_complex_erf(k*xi/2/sqrt 2 + imaginary*alpha/(xi*sqrt 2))
