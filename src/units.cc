#include <src/units.h>
#include <src/defines.h>
#include <src/exception.h>
#include <src/keyword.h>

using namespace gigide;

void
Units::init()
{
    std::string eunitstr = GigideKeyword::getEnergyUnits();
    if (eunitstr == "aj")
        eunits_ = aJ;
    else
        eunits_ = Hartree;
}

double
Units::getConversion(energy_t src, energy_t dst)
{
    if (src == dst)
        return 1.0;

    if (src == aJ and dst == Hartree)
        return 1.0/HARTREE_TO_AJ;
    else if (src == Hartree and dst == aJ)
        return HARTREE_TO_AJ;
    else
        except("Invalid unit conversion");
}

double
Units::convert(double e, energy_t src, energy_t dst)
{
    double conversion = getConversion(src, dst);
    return e * conversion;
}

double
Units::standardize(double e, energy_t units)
{
    return convert(e, units, eunits_);
}

double
Units::convert(double e, energy_t units)
{
    return convert(e, eunits_, units);
}


Units::energy_t Units::eunits_;
bool Units::initdone_ = false;



