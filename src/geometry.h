
#ifndef gigide_geometry_h
#define gigide_geometry_h

#include <math.h>
#include <src/defines.h>
#include <src/molecule.h>
#include <src/deftypes.h>
#include <vector>

namespace gigide {

class XYZPoint : public smartptr::Countable {

    protected:
        ConstVectorPtr xyz_;

    public:
        XYZPoint();

        XYZPoint(ConstVectorPtr xyz);

        XYZPoint(const ConstAtomPtr& atom);

        ConstVectorPtr getXYZ() const;

};

class Midpoint : public XYZPoint {
    
    public:
        Midpoint(
            const Set<ConstAtomPtr>& atoms
        );

};

class Axis : public smartptr::Countable {

    private:
        ConstVectorPtr vec_;
    
    public:
        Axis(const ConstXYZPointPtr& pt1, const ConstXYZPointPtr& pt2);

        Axis(ConstVectorPtr vec);

        ConstVectorPtr getVector() const;
        
};

}

#endif
