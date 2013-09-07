# Installing deft

Getting deft
------------

You can get deft with:

    git clone git://github.com/droundy/deft.git

On the OSU Physics cluster, you can get the code with:

    git clone /home/droundy/git/deft.git

Building deft
-------------

To build, just run:

    scons

Building on an Ubuntu or Debian machine
---------------------------------------

You need ...

    haskell-mode # optional, but handy
    emacs-goodies-el # optional, but handy for markdown-mode
    scons
    fftw3-dev
    libpopt-dev # for command-line argument processing
    g++
    haskell-platform
    texlive-full
    gnuplot
    python-matplotlib
    python-scipy
    python-markdown # for the documentation

Running the test suite
----------------------

You can run the tests with

    scons check

### To write a new test

Create a new file such as

    tests/my-new-feature.cpp

Then add that test to `SConstruct`

    for test in Split(""" memory saft eos  ... my-new-feature  """):
        env.BuildTest(test, all_sources)
    
At this point, `scons check` should do what you want.  There are
several such `for` loops, so that you can choose one that lists the
sources required, in order to minimize the build time for your test.


Using the git pre-commit hook
-----------------------------

`scons` will add a pre-commit hook that causes git to build the code
and papers and run the test suite prior to actually committing your
changes when you run `git commit`.  Thanks to `scons` caching, this
should not take long (a few minutes), but when it took much longer, I
instituted a short-cut to bypass the process.  If you have a few quick
(and safe) commits, you could skip the tests on all but the final one.

To skip the tests, you simply run:

    TEST-none git commit # possibly with -m
