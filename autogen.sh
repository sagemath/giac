#! /bin/sh
aclocal -I /usr/local/share/aclocal
autoheader
autoconf
automake
./configure $*
