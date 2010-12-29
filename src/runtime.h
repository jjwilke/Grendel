#ifndef gigide_runtime_h
#define gigide_runtime_h

#include "gigide.hpp"

namespace gigide {

class GigideRuntime {

    private:
        static XMLArchivePtr archive_;

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

#endif
