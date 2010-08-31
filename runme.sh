#!/usr/bin/env tcsh

git submodule init
git submodule update
aclocal
libtoolize
automake --add-missing
automake
autoconf

./configure --prefix=/opt/gigide/3.0 --with-cxx=g++ --with-libdirs="-L/usr/lib" --with-cc=gcc --with-f77=gfortran --with-cxx-flags="-O2 -g" --with-intder=$HOME/Programs/Intder/intder --with-libs='-lblas -llapack' --with-anharm=$HOME/Programs/Anharm/anharm

