#!/usr/bin/env tcsh

git submodule init
git submodule update
aclocal
glibtoolize
automake --add-missing
automake
autoconf

./configure --prefix=/opt/gigide/3.0 --with-cxx=g++ --with-libdirs="-L/usr/lib" --with-cc=gcc --with-f77=gfortran --with-cxx-flags="-O0 -g" --with-libs='-lblas -llapack' 

