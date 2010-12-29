#ifndef gigide_symmetry_hpp
#define gigide_symmetry_hpp

namespace gigide {

class SymmetryOperation;
class ImproperRotation;
class Rotation;
class Reflection;
class Inversion;
class IdentityElement;
class PointGroupClass;
class PointGroup;

typedef boost::intrusive_ptr<SymmetryOperation> SymmetryOperationPtr;
typedef boost::intrusive_ptr<ImproperRotation> ImproperRotationPtr;
typedef boost::intrusive_ptr<Rotation> RotationPtr;
typedef boost::intrusive_ptr<Reflection> ReflectionPtr;
typedef boost::intrusive_ptr<Inversion> InversionPtr;
typedef boost::intrusive_ptr<IdentityElement> IdentityElementPtr;
typedef boost::intrusive_ptr<PointGroupClass> PointGroupClassPtr;
typedef boost::intrusive_ptr<PointGroup> PointGroupPtr;

typedef boost::intrusive_ptr<const SymmetryOperation> ConstSymmetryOperationPtr;
typedef boost::intrusive_ptr<const ImproperRotation> ConstImproperRotationPtr;
typedef boost::intrusive_ptr<const Rotation> ConstRotationPtr;
typedef boost::intrusive_ptr<const Reflection> ConstReflectionPtr;
typedef boost::intrusive_ptr<const Inversion> ConstInversionPtr;
typedef boost::intrusive_ptr<const IdentityElement> ConstIdentityElementPtr;
typedef boost::intrusive_ptr<const PointGroupClass> ConstPointGroupClassPtr;
typedef boost::intrusive_ptr<const PointGroup> ConstPointGroupPtr;

}

#endif

