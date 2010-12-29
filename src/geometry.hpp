#ifndef gigide_geometry_hpp
#define gigide_geometry_hpp

namespace gigide {

class XYZPoint;
class Axis;
class Midpoint;

typedef boost::intrusive_ptr<XYZPoint> XYZPointPtr;
typedef boost::intrusive_ptr<Axis> AxisPtr;
typedef boost::intrusive_ptr<Midpoint> MidpointPtr;

typedef boost::intrusive_ptr<const XYZPoint> ConstXYZPointPtr;
typedef boost::intrusive_ptr<const Axis> ConstAxisPtr;
typedef boost::intrusive_ptr<const Midpoint> ConstMidpointPtr;

}

#endif

