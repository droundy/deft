module Quantum 
       ( oneElectron )
       where

import CodeGen

oneElectron :: Expression RealSpace
oneElectron = ifft (fft psi * k**2/2) * psi/integrate (psi**2)
      where psi = r_var "x"