
#ifndef _gigide_src_defines_h_
#define _gigide_src_defines_h_

#undef heisenbug
#define heisenbug cout << "Heisenbug: " << __FILE__ << " " << __LINE__ << endl

#define delcheck cout << "Destructor: " << __FILE__ << " " << __LINE__ << endl
#define concheck cout << "Constructor: " << __FILE__ << " " << __LINE__ << endl

#define debug_print cout

#define MAX_LINE_SIZE 300

#define DEFAULT_DISPLACEMENT_SIZE 0.01

#define PRINT_DETAIL 0

#define HARTREE_TO_AJ 4.35974581

#define XYZ_DIM 3

//the conversion from the hartree per bohr units to mks units
#define AU_TO_MKS 2.37492753340298e+28

//the conversion from hertz to wavenumber
#define HERTZ_TO_WN 1.0/2.99792458E10

#define BOHR_TO_ANGSTROM 0.52917720859

#define ANGSTROM_TO_BOHR 1.88972613289

#define RADIAN_TO_DEGREE 57.2957795130823

#define PI 3.14159265358979324

#define DEFAULT_INPUT_FILE "intdif.inp"

#define DEFAULT_INTDER_INPUT_FILE "intder.inp"

#define DEFAULT_INTDER_OUTPUT_FILE "intder.out"

#define DEFAULT_ANHARM_INPUT_FILE "input.dat"

#define SHIFT_DAMP 0.50

#define FIT_POSIX_FILE "fit.posix"

#define DISPITER_POSIX_FILE "dispiter.posix"

#define sigxystr "SIGXY"
#define sigxzstr "SIGXZ"
#define sigyzstr "SIGYZ"
#define c2xstr "C2X"
#define c2ystr "C2Y"
#define c2zstr "C2Z"
#define invopstr "I"
#define eopstr "E"

#define GGF1 138
#define GGF2 36

#define bohr "bohr"
#define angstrom "angstrom"

#endif
