#!/bin/bash

set -ev

cd tests/generated-haskell

for i in *.tex; do
    echo working on $i
    pdflatex $i
done

echo this test passes
