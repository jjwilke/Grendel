#ifndef gigide_molecule_hpp
#define gigide_molecule_hpp

namespace gigide {

class Molecule;
class Atom;

typedef boost::intrusive_ptr<Molecule>   MoleculePtr;
typedef boost::intrusive_ptr<const Molecule>   ConstMoleculePtr;
typedef boost::intrusive_ptr<Atom>   AtomPtr;

typedef boost::intrusive_ptr<const Molecule>   ConstMoleculePtr;
typedef boost::intrusive_ptr<const Atom>   ConstAtomPtr;

}

#endif

