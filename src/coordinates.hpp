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
class CoordinateSubspace;

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

typedef boost::intrusive_ptr<const InternalCoordinate> ConstInternalCoordinatePtr;
typedef boost::intrusive_ptr<const SimpleInternalCoordinate> ConstSimpleInternalCoordinatePtr;
typedef boost::intrusive_ptr<const SymmetryInternalCoordinate> ConstSymmetryInternalCoordinatePtr;
typedef boost::intrusive_ptr<const LinX> ConstLinXPtr;
typedef boost::intrusive_ptr<const LinY> ConstLinYPtr;
typedef boost::intrusive_ptr<const Lin1> ConstLin1Ptr;
typedef boost::intrusive_ptr<const BondAngle> ConstBondAnglePtr;
typedef boost::intrusive_ptr<const BondLength>    ConstBondLengthPtr;
typedef boost::intrusive_ptr<const Torsion> ConstTorsionPtr;
typedef boost::intrusive_ptr<const OutOfPlaneBend>   ConstOutOfPlaneBendPtr;
typedef boost::intrusive_ptr<const PeriodicCoordinate>    ConstPeriodicCoordinatePtr;

typedef boost::intrusive_ptr<CoordinateSubspace>   CoordinateSubspacePtr;
typedef boost::intrusive_ptr<const CoordinateSubspace>   ConstCoordinateSubspacePtr;

}

#endif

