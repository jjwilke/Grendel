#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT(gigide, 1.0, Jeremiah Wilke jjwilke@uga.edu)
#sanity check to ensure that source files exist
AC_CONFIG_SRCDIR([src/main.cc])
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([config.h])

LT_INIT()

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_LIBTOOL

AC_LANG([C++])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/smartptr/src/Makefile
                 lib/Makefile
                 test.py])

AC_ARG_WITH(include,
[  --with-include          Specifies include directories (-Idir1 -Idir2).],
CPPFLAGS=$withval
echo Using extra include directories: $withval
)

AC_ARG_WITH(libs,
[  --with-libs             Specifies libraries (-llib1 -llib2).],
LIBS=$withval
echo Using extra libraries: $withval
)

AC_ARG_WITH(libdirs,
[  --with-libdirs          Specifies library directories (-Ldir1 -Ldir2).],
LDFLAGS=$withval
echo Using extra library directories: $withval
)

AC_ARG_WITH(cxx,
[  --with-cxx              Gives the name of the C++ compiler to use.],
CXX=$withval
)

AC_ARG_WITH(cxx-flags,
[  --with-cxx-flags        Gives the compiler flags for C++],
CXXFLAGS=$withval
)

AC_ARG_WITH(f77,
[  --with-f77              Gives the name of the FORTRAN 77 compiler to use.],
F77=$withval
)

AC_ARG_WITH(cc,
[  --with-cc               Gives the name of the C compiler to use.],
CC=$withval
)

havemkl=no
AC_ARG_WITH(mkl,
[  --with-mkl             Give an mkl library for blas calls],
echo Using mkl for blas and lapack
AC_DEFINE(HAVE_MKL, 1, [Whether to use mkl])
CPPFLAGS="$CPPFLAGS -I$withval/include"
LIBS="-lmkl_lapack -lguide -lpthread $LIBS"
mkl=$withval
havemkl=yes,
AC_DEFINE(HAVE_MKL, 0, [Whether to use mkl])
)

AC_ARG_WITH(boost,
[ --with-boost          Gives the root folder for boost],
CXXFLAGS="$CXXFLAGS -I/$withval",
echo Using standard include folder for boost
)

if test X$havemkl = Xyes ; then
    AC_ARG_WITH(mkl-arch,
        [  --with-mkl-arch            Give an mkl architecture (e.g. em64t, 32).  This defines lib folder],
        LDFLAGS="$LDFLAGS -L$mkl/lib/$withval"
        LIBS="$LIBS -lmkl_$withval",
        LDFLAGS="$LDFLAGS -L$mkl/lib"
        LIBS="$LIBS -lmkl"
      )
fi
if test X$havemkl = Xyes ; then
LIBS="-lguide -lpthread $LIBS"
fi

#check for python
AC_CHECK_PROGS([PYTHON], python, [:])
if test "$PYTHON" = :; then
    AC_MSG_ERROR([Package requires at least Python v2.4])
fi
PYVER=esyscmd(""./getpyver"")
PYDIST="python$PYVER"
echo checking python version... found $PYDIST
PYBASE=esyscmd(""./getpyinclude"")
PYDEFINC=$PYBASE/$PYDIST
PYDEFLIBDIR=esyscmd(""./getpylib"")
echo found python include dir $PYDEFINC

AC_ARG_WITH(python-include,
[  --with-python-include   Gives the path to the python include directory.] ,
PYINCLUDE=$withval,
PYINCLUDE=$PYDEFINC
)

AC_ARG_WITH(python-libdir,
[  --with-python-libdir    Gives the path to the python include directory.] ,
PYLIBDIR=$withval,
PYLIBDIR=$PYDEFLIBDIR
)

AC_ARG_WITH(intder,
[  --with-intder   Gives the path to intder.] ,
intderpath=$withval
hasintder=1,
intderpath="none"
hasintder=0
)
AC_DEFINE_UNQUOTED(INTDER_PATH,"${intderpath}", [The location of intder])
AC_DEFINE_UNQUOTED(HAS_INTDER,${hasintder}, [Whether intder exists])
AC_SUBST(HAS_INTDER, $hasintder)

AC_ARG_WITH(anharm,
[  --with-anharm   Gives the path to anharm.] ,
anharmpath=$withval
hasanharm=1,
anharmpath="none"
hasanharm=0
)
AC_DEFINE_UNQUOTED(ANHARM_PATH,"${anharmpath}", [The location of anharm])
AC_DEFINE_UNQUOTED(HAS_ANHARM,${hasanharm}, [Whether anharm exists])
AC_SUBST(HAS_ANHARM, $hasanharm)

LIBS="$LIBS -l$PYDIST"
CPPFLAGS="$CPPFLAGS -I$PYINCLUDE"
LDFLAGS="$LDFLAGS -L$PYLIBDIR"

#PKG_CHECK_MODULES([DEPS], libxml++-2.6 >= 2.14.0)
#AC_SUBST(DEPS_CFLAGS)
#AC_SUBST(DEPS_LIBS)

gigdatadir="$prefix/lib"
AC_DEFINE_UNQUOTED(GIGDATADIR,"${gigdatadir}",[The directory into which all data files will be installed])

AC_DEFINE([PROGRAM],"GIGID-E","The name of the program")

AC_CHECK_HEADER([boost/intrusive_ptr.hpp],
,
AC_MSG_ERROR([Package requires boost])
)

AC_SEARCH_LIBS(DGEMM,,
echo Found blas library,
AC_MSG_ERROR([Need blas library to continue])
)


AC_OUTPUT

