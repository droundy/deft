{- FMT defines the fundamental measures that are the raw input of
fundamental measure theory.  The FMT free energy functional itself is
defined in WhiteBear (or other modules for other versions thereof). -}

module FMT
       ( n, n3, n2, n2p, n1, n0, rad,
         n2v, n1v, n2vp,
         sqr_n2v, n1v_dot_n2v,
         shell, shell_diam, step, vshell, vshelldot,
         shellPrime, vshellPrime, vshellPrimedot,
         dr )
       where

import Expression

rad :: Type a => Expression a
rad = s_var "R"

kdr :: Expression KSpace
kdr = k * dr

kR :: Expression KSpace
kR = k * rad

smear :: Expression KSpace
smear = exp (-6.0*kdr*kdr)

n, n3, n2, n1, n0, n2p, n1v_dot_n2v, sqr_n2v :: Expression RealSpace
n = "n" === r_var "x"
n3 = "n3" === step n
n2 = "n2" === shell n
n2p = "n2p" === shellPrime n
n1 = "n1" === shell n / (4*pi*rad)
n0 = "n0" === shell n / (4*pi*rad**2)

n2v, n1v, n2vp :: Vector RealSpace
n2v = "n2v" `nameVector` vshell n
n1v = "n1v" `nameVector` vshell n /. (4*pi*rad)
n2vp = "n2vp" `nameVector` vshellPrime n

sqr_n2v = var "n2vsqr" "{\\left|\\vec{n}_{2v}\\right|^2}" (n2v `dot` n2v)
n1v_dot_n2v = var "n1v_dot_n2v" "{\\vec{n}_{1v}\\cdot\\vec{n}_{2v}}" (n2v `dot` n1v)

shell_diam, shell, step, shellPrime :: Expression RealSpace -> Expression RealSpace
shell x = ifft ( deltak * fft x)
  where deltak = protect "deltak" "\\delta(k)" $ smear * (4*pi) * rad * sin kR / k
shell_diam x = ifft ( deltak * fft x)
  where deltak = protect "deltak2" "\\delta_{2R}(k)" $ smear * (4*pi) * 2*rad * sin (2*kR) / k
step x = ifft ( stepk * fft x)
  where stepk = protect "step" "\\Theta(k)" $ smear * (4*pi) * (sin kR - kR * cos kR) / k**3

vshell, vshellPrime :: Expression RealSpace -> Vector RealSpace
vshell n = vifft $ deltav *. fft n
  where deltav = vprotect "delta" "\\delta(k)" $
                 smear * (4*pi) * imaginary * (rad * cos kR - sin kR/k)/k**2 .* kvec

shellPrime x = ifft ( delta' * fft x)
  where delta' = protect "deltaprime" "\\delta'(k)" $
                 smear * (4*pi) * ( -sin kR/k - rad * cos kR)
vshellPrime x = vifft ( delta' *. fft x)
  where delta' = vprotect "deltaprime" "\\delta'(k)" $
                 smear * (4*pi) * imaginary * ( rad * sin kR)/k .* kvec

vshellPrimedot, vshelldot :: Vector RealSpace -> Expression RealSpace
vshelldot n = ifft $ deltav `dot` vfft n
  where deltav = vprotect "delta" "\\delta(k)" $
                 smear * (4*pi) * imaginary * (rad * cos kR - sin kR/k)/k**2 .* kvec
vshellPrimedot x = ifft ( delta' `dot` vfft x )
  where delta' = vprotect "deltaprime" "\\delta'(k)" $
                 smear * (4*pi) * imaginary * ( rad * sin kR)/k .* kvec
