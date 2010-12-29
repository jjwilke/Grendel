PROGRAM = anharm
FSRC = anharm.F anhasy.F  anhlin.F  anhsym.F \
       anhdia.F  anhsph.F anhlib.F

FC = gfortran
FFLAGS = -g 

LIBS = 

INCLUDES = -I.

BINOBJ = $(FSRC:%.F=%.o)

$(PROGRAM): $(BINOBJ) $(LIBS)
	$(FC) $(FFLAGS) $^ -o $(PROGRAM)

%.o: %.F
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

clean:
	/bin/rm -f $(PROGRAM) *.o core *~
