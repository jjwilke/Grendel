
#include <src/geometry.h>
#include <src/exception.h>

using namespace gigide;
using namespace std;
using namespace smartptr;

XYZPoint::XYZPoint(
    ConstVectorPtr xyz
) : xyz_(xyz)
{
}

XYZPoint::XYZPoint(
)
{
}

XYZPoint::XYZPoint(
    const ConstAtomPtr& atom
) : xyz_(atom->getXYZ())
{
}

ConstVectorPtr
XYZPoint::getXYZ() const
{
    return xyz_;
}

Midpoint::Midpoint(
    const Set<ConstAtomPtr>& atoms
)
{
    if (!atoms.size())
    {
        except("Cannot compute midpoint with no atoms");
    }
    
    VectorPtr xyz = atoms[0]->getXYZ().copy();
    for (int i=1; i < atoms.size(); ++i)
        xyz.accumulate(atoms[i]->getXYZ());
    xyz.scale(1.0/atoms.size());

    xyz_ = xyz;
}

Axis::Axis(
    const ConstXYZPointPtr& pt1,
    const ConstXYZPointPtr& pt2
)
{
    vec_ = pt2->getXYZ() - pt1->getXYZ();
}

Axis::Axis(
    ConstVectorPtr vec
) : vec_(vec)
{
}

ConstVectorPtr
Axis::getVector() const
{
    return vec_;
}
