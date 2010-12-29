#ifndef gigide_coordinates_hpp
#define gigide_coordinates_hpp

namespace gigide {

class InternalCoordinate;
class SimpleInternalCoordinate;
class SymmetryInternalCoordinate;
class LinX;
class LinY;
class Lin1;
class BondAngle;
class BondLength;
class Torsion;
class OutOfPlaneBend;
class PeriodicCoordinate;

typedef boost::intrusive_ptr<InternalCoordinate>          InternalCoordinatePtr;
typedef boost::intrusive_ptr<SimpleInternalCoordinate>    SimpleInternalCoordinatePtr;
typedef boost::intrusive_ptr<SymmetryInternalCoordinate>  SymmetryInternalCoordinatePtr;
typedef boost::intrusive_ptr<LinX>                        LinXPtr;
typedef boost::intrusive_ptr<LinY>                        LinYPtr;
typedef boost::intrusive_ptr<Lin1>                        Lin1Ptr;
typedef boost::intrusive_ptr<BondAngle>                   BondAnglePtr;
typedef boost::intrusive_ptr<BondLength>                  BondLengthPtr;
typedef boost::intrusive_ptr<Torsion>                     TorsionPtr;
typedef boost::intrusive_ptr<OutOfPlaneBend>              OutOfPlaneBendPtr;
typedef boost::intrusive_ptr<PeriodicCoordinate>          PeriodicCoordinatePtr;

}

#endif

