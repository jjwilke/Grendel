#ifndef gigide_derivative_hpp
#define gigide_derivative_hpp

namespace gigide {

class Derivative;
class DerivativeIterator;

typedef boost::intrusive_ptr<Derivative>         DerivativePtr;
typedef boost::intrusive_ptr<DerivativeIterator> DerivativeIteratorPtr;

typedef boost::intrusive_ptr<const Derivative> ConstDerivativePtr;
typedef boost::intrusive_ptr<const DerivativeIterator>   ConstDerivativeIteratorPtr;

}

#endif

