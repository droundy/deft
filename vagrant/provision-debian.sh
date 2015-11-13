#!/usr/bin/env bash

set -ev

apt-get update
apt-get install -y git python3 ruby-sass
apt-get install -y graphviz
apt-get install -y build-essential fftw3-dev libpopt-dev
# apt-get install -y scons
apt-get install -y haskell-platform
apt-get install -y gnuplot-nox
apt-get install -y texlive-latex-base texlive-publishers texlive-science
apt-get install -y feynmf
apt-get install -y python-matplotlib python-scipy python-markdown
apt-get install -y python3-matplotlib python3-scipy python3-markdown
apt-get install -y python-sympy
apt-get install -y inkscape # for the documentation

git clone git://github.com/droundy/fac
cd fac
sh build-linux.sh
./fac
cp fac /usr/local/bin/

echo We finished setting up fac.
