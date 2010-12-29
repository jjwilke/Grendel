#ifndef gigide_fit_hpp
#define gigide_fit_hpp

namespace gigide {

class Fit;
class FitFactory;
class FitPoint;

typedef boost::intrusive_ptr<Fit> FitPtr;
typedef boost::intrusive_ptr<FitFactory> FitFactoryPtr;
typedef boost::intrusive_ptr<FitPoint> FitPointPtr;

typedef boost::intrusive_ptr<const Fit> ConstFitPtr;
typedef boost::intrusive_ptr<const FitFactory> ConstFitFactoryPtr;
typedef boost::intrusive_ptr<const FitPoint> ConstFitPointPtr;

}

#endif

