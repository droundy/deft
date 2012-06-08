#!/bin/sh

set -ev

# Here we will set up git hooks to do nice things...
if test -d git && test -d .git/hooks; then
    for i in `ls git | grep -v '~'`; do
        echo Setting up $i hook...
        ln -sf ../../git/$i .git/hooks/$i
    done
else
    echo We do not seem to be in a git repository.
fi

# automake wants a file called README.
cp README.md README

aclocal
autoheader
automake --add-missing
autoconf

WARNINGEXCEPTIONS=' -Wno-unused-variable -Wno-unused-parameter -Wno-return-type '

export CXXFLAGS="-ansi -W -Wall $WARNINGEXCEPTIONS -Werror -pipe -O2 -DNDEBUG"
export CXXFLAGS="-ansi -W -Wall $WARNINGEXCEPTIONS -Werror -pipe -O2"
if env | grep CCACHE_; then
    echo Using ccache to speed up compilation.
    CXX='ccache g++' ./configure
else
    echo 'Consider installing ccache to speed up compilation!'
    ./configure
fi

set +v
echo If you are not a deft developer, you should now run ./configure
echo to avoid overly-pedantic compiler errors!
