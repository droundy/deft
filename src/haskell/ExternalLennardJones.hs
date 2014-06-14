-- This module extends a functional to have an external potential
-- describing interaction with a Lennard-Jones atom.  It also redefines
-- the density in terms of an effective potential

module ExternalLennardJones
       ( lennard_jones, lennard_jones_density, lennard_jones_potential,
         lennard_jones_hughes, lennard_jones_water_saft,
         lennard_jones_water_saft_by_hand,
         ljepsilon, ljsigma )
       where

import FMT ( n )
import WhiteBear ( kT )
import WaterSaft ( water_saft, water_saft_by_hand )
import HughesSaft ( saft_fluid )
import Expression

ljepsilon, ljsigma :: Type a => Expression a
ljepsilon = s_tex "ljepsilon" "\\epsilon_{LJ}"
ljsigma = s_tex "ljsigma" "\\sigma_{LJ}"

lennard_jones :: Expression Scalar -> Expression Scalar
lennard_jones f = substitute n ljn f + integrate (ljn*v)
  where v = lennard_jones_potential
        ljn = lennard_jones_density

lennard_jones_potential :: Expression RealSpace
lennard_jones_potential = ljepsilon*( (ljsigma/rmag)**12 - (ljsigma/rmag)**6 )

lennard_jones_density :: Expression RealSpace
lennard_jones_density = "n" === exp(-(r_var "Veff")/kT)
-- lennard_jones_density = exp(-(r_var "Veff" + lennard_jones_potential)/kT)

lennard_jones_hughes :: Expression Scalar
lennard_jones_hughes = lennard_jones saft_fluid

lennard_jones_water_saft :: Expression Scalar
lennard_jones_water_saft = lennard_jones water_saft

lennard_jones_water_saft_by_hand :: Expression Scalar
lennard_jones_water_saft_by_hand = lennard_jones water_saft_by_hand
