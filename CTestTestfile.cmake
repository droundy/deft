# CMake generated Testfile for 
# Source directory: /mnt/home/droundy/src/deft
# Build directory: /mnt/home/droundy/src/deft
# 
# This file includes the relevent testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(run-pair-monte-carlo "pair-monte-carlo" "10" "1000000" "0.01" "/tmp/foo" "/tmp/dafoo" "periodxy" "20" "wallz" "20" "flatdiv")
ADD_TEST(run-monte-carlo "monte-carlo" "10" "100000" "0.01" "/tmp/test.out")
ADD_TEST(surface-tension "tests/surface-tension.test")
ADD_TEST(functional-arithmetic "tests/functional-arithmetic.test")
ADD_TEST(saft "tests/saft.test")
ADD_TEST(eos "tests/eos.test")
ADD_TEST(eps "tests/eps.test")
ADD_TEST(fftinverse "tests/fftinverse.test")
ADD_TEST(ideal-gas "tests/ideal-gas.test")
ADD_TEST(precision "tests/precision.test")
ADD_TEST(print-iter "tests/print-iter.test")
ADD_TEST(convolve-finite-difference "tests/convolve-finite-difference.test")
ADD_TEST(new-hard-spheres "tests/new-hard-spheres.test")
ADD_TEST(memory "tests/memory.test")
ADD_TEST(functional-of-double "tests/functional-of-double.test")
ADD_TEST(new-code "tests/new-code.test")
ADD_TEST(convolve "tests/convolve.test")
ADD_TEST(generated-code "tests/generated-code.test")
