set -ev

# automake wants a file called README.
cp README.mdwn README

aclocal
autoheader
automake --add-missing || true #workaround for buggy old automake
autoconf
