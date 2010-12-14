#!/bin/bash

set -ev

diff -u comparisons/eta_effective*.dat

tests/print-iter.test > tests/print-iter.sh.log.log

diff -u tests/print-iter.sh.correct tests/print-iter.sh.log.log

echo this test passes
