#ifndef gigide_displacement_hpp
#define gigide_displacement_hpp

namespace gigide {

class DisplacementMapping;
class Displacement;
class DisplacementIterator;

typedef boost::intrusive_ptr<DisplacementMapping>    DisplacementMappingPtr;
typedef boost::intrusive_ptr<Displacement> DisplacementPtr;
typedef boost::intrusive_ptr<DisplacementIterator>    DisplacementIteratorPtr;

typedef boost::intrusive_ptr<const DisplacementMapping>    ConstDisplacementMappingPtr;
typedef boost::intrusive_ptr<const Displacement> ConstDisplacementPtr;
typedef boost::intrusive_ptr<const DisplacementIterator>    ConstDisplacementIteratorPtr;

}

#endif

