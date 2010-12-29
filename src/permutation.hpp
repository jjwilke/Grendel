#ifndef gigide_permutation_hpp
#define gigide_permutation_hpp

#include "permutation.h"

namespace gigide {

typedef boost::intrusive_ptr<PermutationGenerator<int> >   PermutationGeneratorIntPtr;
typedef boost::intrusive_ptr<PermutationGenerator<double> >   PermutationGeneratorDoublePtr;

typedef boost::intrusive_ptr<const PermutationGenerator<int> >   ConstPermutationGeneratorIntPtr;
typedef boost::intrusive_ptr<const PermutationGenerator<double> >   ConstPermutationGeneratorDoublePtr;

}

#endif
