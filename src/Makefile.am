bin_PROGRAMS=grendel
grendel_SOURCES=coordinates.cc molecule.cc main.cc fileio.cc qff.cc keyword.cc utilities.cc hessian.cc taylor.cc fit.cc input.cc geometry.cc symmetry.cc exception.cc gigmatrix.cc displacement.cc runtime.cc optimize.cc derivative.cc units.cc 
INCLUDES = $(DEPS_CFLAGS)
LDADD = $(DEPS_LIBS)
grendel_LDADD = smartptr/src/libsmartptr.la  smartptr/src/libs/regex/src/libregexp.la smartptr/src/libs/smart_ptr/src/libptr.la
grendel_LDFLAGS = -static
