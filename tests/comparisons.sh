#!/bin/bash

set -ev

for i in eta_effective gHSetaeff a1_sw a2_sw; do
    echo working on $i
    echo diff --side-by-side comparisons/${i}_vrpack.dat comparisons/${i}_deft.dat
    diff --side-by-side comparisons/${i}_deft.dat comparisons/${i}_vrpack.dat
done

echo this test passes
