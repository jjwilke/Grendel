#ifndef gigide_units_h
#define gigide_units_h


class Units {

    public:
        typedef enum {aJ, Hartree} energy_t;

    private:
        static void init();

        static double getConversion(energy_t src, energy_t dst);

        static energy_t eunits_;

        static bool initdone_;

    public:

        static double convert(double e, energy_t src, energy_t dst);

        static double standardize(double e, energy_t units);

        static double convert(double e, energy_t units);

};


#endif

