#!/bin/bash

set -ev

for i in eta_effective gHSetaeff a1_sw a2_sw da1dlam gsw delta x assoc; do
    echo working on $i
    echo diff --side-by-side comparisons/${i}_deft.dat comparisons/${i}_vrpack.dat
    diff --side-by-side comparisons/${i}_deft.dat comparisons/${i}_vrpack.dat
done

echo this test passes
