#ifndef gigide_taylor_hpp
#define gigide_taylor_hpp

namespace gigide {

class TaylorTerm;
class TaylorSeriesEnergy;

typedef boost::intrusive_ptr<TaylorTerm>   TaylorTermPtr;
typedef boost::intrusive_ptr<TaylorSeriesEnergy>   TaylorSeriesEnergyPtr;

typedef boost::intrusive_ptr<const TaylorTerm>   ConstTaylorTermPtr;

}

#endif

