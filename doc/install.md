# Installing deft

[TOC]

## Getting deft

You can get deft with:

    git clone git://github.com/droundy/deft.git

On the OSU Physics cluster, you can get the code with:

    git clone /home/droundy/git/deft.git

Or if you have an account on the Roundy computers, you can get the
code onto your own computer with:

    git clone username@knightley.physics.oregonstate.edu:/home/droundy/git/deft.git

where `username` is your username on the cluster.

## Building deft

You need to start by getting fac by following the instructions at the
[fac website](https://physics.oregonstate.edu/~roundyd/fac/building.html).
You should probably put fac in your path using something like

    sudo cp fac /usr/local/bin/

Once you have fac available, you can build deft by typing `fac`
anywhere in the deft directory.

## Building on an Ubuntu or Debian machine

You need ...

    fftw3-dev
    libpopt-dev # for command-line argument processing
    g++
    haskell-platform
    texlive-full
    python-matplotlib
    python-scipy python3-scipy
    python-sympy python3-sympy
    python-markdown python3-markdown # for the documentation
    inkscape # for the documentation
