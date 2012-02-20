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

#CXXFLAGS='-ansi -pipe -W -g -Wall -O2 -Werror' ./configure
#CXXFLAGS='-ansi -pipe -W -Wall -O2 -Werror' ./configure
export CXXFLAGS='-ansi -pipe -W -Wall -O2 -ffast-math -Werror'
export CXXFLAGS='-ansi -pipe -W -Wall -O2 -Werror'
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
