
#include <src/molecule.h>
#include <src/coordinates.h>
#include <src/deftypes.h>

namespace gigide {

/**
    Predicts the force constant between two simple internal coordinates 

    @class EmpiricalHessianTerm 
*/
class EmpiricalHessianTerm : public smartptr::Countable {

    protected:
        /**
            The first internal coordinate in the force constant        
        */
        ConstSimpleInternalCoordinatePtr c1_;

        /**
            The second internal coordinate in the force constant        
        */
        ConstSimpleInternalCoordinatePtr c2_;

        /**
            Evaluates whether two atoms are the same in a given internal coordinate.
            Two indices are given relative with internal coordinate numbering.  For example,
            suppose you have two bond angles in methane: a1 = H2-C1-H3, a2 = H3-C1-H4.
            a1[0] = H2, a1[1] = C1, a1[2] = H3, a2[0] = H3, a2[1] = C1, a2[2] = H4 in the
            inernal coordinate numbering scheme.  The method returns the following: <br>
            matches(0,0) = false <br>
            matches(1,1) = true <br>
            matches(2,0) = true <br>
            matches(2,2) = false <br>

            @param id1 Atom index in coordinate 1
            @param id2 Atom index in coordinate 2
            @param at Returns the atom if a match is found
        */
        bool matches(
            int id1, 
            int id2,
            ConstAtomPtr& at
        ) const;

        /**
            Contructor for empirical hessian value between two coordinates.
            This is only called from derived class constructor.

            @param coord1
            @param coord2
        */
        EmpiricalHessianTerm(
            const ConstSimpleInternalCoordinatePtr& coord1,
            const ConstSimpleInternalCoordinatePtr& coord2
        );

        /**
            Default constructor necessary for terms 
            that do not depend on internal coordinates.
        */
        EmpiricalHessianTerm();

    public:
        /**
            Retrieve the value of the term.
            @return force constant value
        */
        virtual double getValue() const = 0;

        /**
            Given two internal coordinates, configure an empirical
            hessian class to compute the force constant. For example,
            if two bonds are passed to getTerm, this would
            return a BondBond_EmpiricalHessian.

            @param coord1
            @param coord2
            @return The term to compute empirical hessian force constant
        */
        static EmpiricalHessianTermPtr getTerm(
            const ConstSimpleInternalCoordinatePtr& coord1,
            const ConstSimpleInternalCoordinatePtr& coord2
        );
};

/**
    Encapsulates a hessian term that always returns a constant, usually 0.
    
    @class ConstantTerm
*/
class ConstantTerm : public EmpiricalHessianTerm {

    private:
        friend class EmpiricalHessianTerm;

        double val_;

        ConstantTerm(double val);

    public:
        double getValue() const;
};

/**
    Encapsulates a force constant between two bonds

    @class BondBond_EmpiricalHessian
*/
class BondBond_EmpiricalHessian : public EmpiricalHessianTerm {

    private:
        friend class EmpiricalHessianTerm;

        /**
            Return the diagonal force constant for a bond stretch
            @param at1 Atom 1 in the bond
            @param at2 Atom 2 in the bond
            @return force constant
        */
        double xxtype(const ConstAtomPtr& at1, const ConstAtomPtr& at2) const;

        /**
            Return the force constant for two bond stretches that
            share one common atom, for example in water the
            force constant between bond stretches R(O1-H2) and R(O1-H3) 
            would be xytpe
            @param at1 The common atom
            @param at2 The other atom in bond 1
            @param at3 The other atom in bond 2
            @return force constant
        */
        double xytype(const ConstAtomPtr& at1, const ConstAtomPtr& at2, const ConstAtomPtr& at3) const;

        /**
            @param bond1
            @param bond2
        */
        BondBond_EmpiricalHessian(
            const ConstSimpleInternalCoordinatePtr& coord1,
            const ConstSimpleInternalCoordinatePtr& coord2
        );
    
    public:

        /**
        */
        double getValue() const;
};

/**
    @class BendBend_EmpiricalHessian
*/
class BendBend_EmpiricalHessian : public EmpiricalHessianTerm {

    private:
        friend class EmpiricalHessianTerm;

        /**
            Compute diagonal force constant for equivalent bends
            for bond angle at1-at2-at3.

            @param at1 Terminal atom in bend
            @param at2 Middle atom in bend
            @param at3 Terminal atom in bend
        */
        double xxxtype(
            const ConstAtomPtr& at1,
            const ConstAtomPtr& at2,
            const ConstAtomPtr& at3
        ) const;
    
        /**
            The constructor takes simple internal coordinate pointers
            in order 
            @param bend1
            @param bend2
        */
        BendBend_EmpiricalHessian(
            const ConstSimpleInternalCoordinatePtr& bend1,
            const ConstSimpleInternalCoordinatePtr& bend2
        );

    public:


        /**
        */
        double getValue() const;
};

/**
    @class BondBend_EmpiricalHessian
*/
class BondBend_EmpiricalHessian : public EmpiricalHessianTerm {

    private:
        friend class EmpiricalHessianTerm;

        /**
        */
        double xxytype(
            const ConstAtomPtr& at1,
            const ConstAtomPtr& at2,
            const ConstAtomPtr& at3
        ) const;

        /**
        */
        BondBend_EmpiricalHessian(
            const ConstSimpleInternalCoordinatePtr& coord1,
            const ConstSimpleInternalCoordinatePtr& coord2
        );
    
    public:

        /**
        */
        double getValue() const;
};

/**
    @class TorsTors_EmpiricalHessian
*/
class TorsTors_EmpiricalHessian : public EmpiricalHessianTerm {

    private:
        friend class EmpiricalHessianTerm;

        /**
            Compute diagonal force constant for a torsion
            @param at1
            @param at2
            @param at3
            @param at4
        */
        double xxxxtype(
            const ConstAtomPtr& at1,
            const ConstAtomPtr& at2,
            const ConstAtomPtr& at3,
            const ConstAtomPtr& at4
        ) const;

        /**
        */
        TorsTors_EmpiricalHessian(
            const ConstSimpleInternalCoordinatePtr& coord1,
            const ConstSimpleInternalCoordinatePtr& coord2
        );

    public:

        /**
        */
        double getValue() const;
};

} //end namespace sc
