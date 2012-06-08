module Quantum 
       ( oneElectron, hydrogenPotential )
       where

import CodeGen

hydrogenPotential :: Expression RealSpace
hydrogenPotential = ifft (setkzero 0 $ -4*pi/k**2)

oneElectron :: Expression RealSpace -> Expression RealSpace
oneElectron potential = (ifft (fft psi * k**2/2) + potential*psi)* psi/integrate (psi**2)
      where psi = r_var "x"