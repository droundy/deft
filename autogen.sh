set -ev

aclocal
autoheader
automake --add-missing || true #workaround for buggy old automake
autoconf
