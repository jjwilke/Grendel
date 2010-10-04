#include <pyregexp.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <config.h>
#include <stdio.h>
#include <src/printstream.h>
#include <src/utilities.h>
#include <src/fileio.h>
#include <src/qff.h>
#include <src/coordinates.h>
#include <src/keyword.h>
#include <src/defines.h>
#include <src/input.h>
#include <src/symmetry.h>
#include <src/exception.h>
#include <src/displacement.h>
#include <src/derivative.h>
#include <src/deftypes.h>
#include <src/gigmatrix.h>

#include <src/smartptr/src/set.h>
#include <src/smartptr/src/xmlarchive.h>

namespace gigide {

class GigideRuntime {

    private:
        static ArchivePtr archive_;

        static MoleculePtr mol_;

        static std::vector<SimpleInternalCoordinatePtr> simples_;

        static std::vector<InternalCoordinatePtr> coords_;

        static ForceFieldPtr qff_;

        static GigideInputFilePtr input_;

        static KeywordSetPtr keymap_;
    
    public:
        
        static void run_intder();
        
        static void run_anharm();

        static void xml_commit();

        static void init_input();

        static void init_keymap();

        static void init_qff();

        static void init_molecule();

        static void init_coordinates();

        static void init_symmetry();

        static void run_dispcart();

        static void run_calc();

        static void run_xmlprint();

        static void run_checkfit();

        static void run_bvectest();

        static void run_displace();

        static void run_generate();

        static ConstMoleculePtr molecule();

};

}

