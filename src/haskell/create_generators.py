#!/usr/bin/python2

for name, module, hsfunctional, inputs in [
    # The following are just for testing purposes
    ("ExternalPotentialTest", "ExternalPotentialTest", "external_potential", '[(ER $ r_var "n", ER 1)]'),
    ("Quadratic", "Quadratic", "quadratic", '[(ER $ r_var "x", ER 1)]'),
    ("QuadraticN0", "QuadraticN0", "quadratic_n0", '[(ER $ r_var "x", ER 1)]'),
    ("QuadraticGaussian", "QuadraticGaussian", "quadratic_gaussian", '[(ER $ r_var "x", ER 1)]'),
    ("LogN0", "LogN0", "log_n0", '[(ER $ r_var "x", ER 1)]'),
    ("Phi1", "WhiteBear", "kTphi1", '[(ER $ r_var "x", ER 1)]'),
    ("Phi2", "WhiteBear", "kTphi2", '[(ER $ r_var "x", ER 1)]'),
    ("Phi3", "WhiteBear", "kTphi3", '[(ER $ r_var "x", ER 1)]'),
    ("SPhi1", "SFMT", "phi1", '[(ER $ r_var "x", ER 1)]'),
    ("SPhi2", "SFMT", "phi2", '[(ER $ r_var "x", ER 1)]'),
    ("SPhi3", "SFMT", "phi3", '[(ER $ r_var "x", ER 1)]'),
    # The rest are "real" functionals of sorts
    ("HomogeneousWhiteBear", "WhiteBear", "homogeneous_whitebear", '[]'),
    ("HomogeneousWhiteBearFluid", "WhiteBear", "homogeneous_whitebear_fluid", '[]'),
    ("WhiteBear", "WhiteBear", "whitebear_n", '[(ER $ r_var "n", ER 1)]'),
    ("WhiteBearFluid", "WhiteBear", "whitebear_fluid_n", '[(ER $ r_var "n", ER 1)]'),
    ("WhiteBearFluidVeff", "WhiteBear", "whitebear_fluid_Veff",
           '[(ER $ r_var "Veff", ER (exp(-r_var "Veff"/s_var "kT")))]'),
    #("SW_dispersion", "SW_liquid", "sw_dispersion", '[(ER $ r_var "n", ER 1)]'),
    ("SFMTFluid", "SFMT", "sfmt_fluid_n", '[(ER $ r_var "n", ER 1)]'),
    ("SFMTFluidVeff", "SFMT", "sfmt_fluid_Veff",
           '[(ER $ r_var "Veff", ER (exp(-r_var "Veff"/s_var "kT")))]'),
    ("HomogeneousSFMTFluid", "SFMT", "homogeneous_sfmt_fluid", '[]'),
    ("SW_liquid", "SW_liquid", "sw_liquid_n", '[(ER $ r_var "n", ER 1)]'),
    ("SW_liquidVeff", "SW_liquid", "sw_liquid_Veff",
           '[(ER $ r_var "Veff", ER (exp(-r_var "Veff"/s_var "kT")))]'),
    ("HomogeneousSW_liquid", "SW_liquid", "homogeneous_sw_liquid", '[]'),
    ("WaterSaft", "WaterSaft", "water_saft_n", '[]'), # no gradients:  for debugging!
    ("WaterSaftByHand", "WaterSaft", "water_saft_by_hand_n", '[]'), # no gradients:  for debugging!
    ("HomogeneousWaterSaft", "WaterSaft", "homogeneous_water_saft_n", '[]'),
    ("HomogeneousWaterSaftByHand", "WaterSaft", "homogeneous_water_saft_by_hand_n", '[]')]:
    # I'm sloppy and just recreate the generate_%s.hs files every time
    f = open('generate_%s.hs' % name, "w")
    f.write("""import NewCode
import %s ( %s )

main :: IO ()
main = createHeaderAndCppFiles %s %s "%s"
""" % (module, hsfunctional, hsfunctional, inputs, name))
    f.close()
