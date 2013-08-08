# Installing deft

Getting deft
============

You can get deft with:

    git clone git://github.com/droundy/deft.git

On the OSU Physics cluster, you can get the code with:

    git clone /home/droundy/git/deft.git

Building deft
=============

To build, just run:

    scons

Building on Ubuntu Precise Pangolin
===================================

You need ...

    haskell-mode # optional, but handy
    fftw3-dev
    g++
    gfortran
    ghc
    libghc-hunit-dev
    texlive-full
    gnuplot
    python-matplotlib
    python-scipy
    python-markdown # for the documentation

Running the test suite
======================

You can run the tests with

    scons check

To write a new test
-------------------

Create a new file such as

    tests/my-new-feature.cpp

Then add that test to `SConstruct`

    for test in Split(""" memory saft eos  ... my-new-feature  """):
        env.BuildTest(test, all_sources)
    
At this point, `scons check` should do what you want.  There are
several such `for` loops, so that you can choose one that lists the
sources required, in order to minimize the build time for your test.


Using the git pre-commit hook
=============================

The autogen.sh script will add a pre-commit hook that causes git to
build the code and papers and run the test suite prior to actually
committing your changes when you run `git commit`.  Unfortunately,
this is a multi-hour process, so I've instituted a few short-cuts.
Ideally, you'd run the full tests once for each day that you make
changes (presumably on the last change you commit for the day) to make
sure you haven't broken anything---usually by forgetting to add a new
file.

To only compile the code, you can run

    TEST=compile git commit # possibly with -m

This will probably take somewhere between 20 and 40 minutes, but won't
catch any mistakes in the papers, which is a bit of a bummer, since
the most common mistakes are the failure to add a `.dat` file.  To
catch these mistakes and still save a bit of time, you can run:

    TEST=build git commit # possibly with -m

And finally, for really simple changes (or if you have many changes to
commit), you can use

    TEST=none git commit # possibly with -m

to skip all the tests.  This is probably the most common option you
use, but if you've got just one commit to make, please consider
letting the entire test suite run, and you can push the next day!
