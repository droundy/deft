set -ev

# automake wants a file called README.
cp README.md README

aclocal
autoheader
automake --add-missing || true #workaround for buggy old automake
autoconf

#CXXFLAGS='-ansi -pipe -W -g -Wall -O2 -Werror' ./configure
CXXFLAGS='-ansi -pipe -W -Wall -O2 -Werror' ./configure

set +v
echo If you are not a deft developer, you should now run ./configure
echo to avoid overly-pedantic compiler errors!
