#ifndef gigide_coordinate_h_
#define gigide_coordinate_h_

#include <src/deftypes.h>
#include <math.h>
#include <vector>
#include <src/defines.h>
#include <src/gigmatrix.h>
#include <iostream>

#include <src/smartptr/src/set.h>
#include <src/smartptr/src/printstream.h>

namespace gigide {

class Molecule;
class Atom;
class SimpleInternalCoordinate;
class SymmetryInternalCoordinate;
class SymmetryOperation;
class PointGroup;

/** 
    @class InternalCoordinate

    Encapsulatesan internal coordinate.  This is capable of computing values and B std::vectors 
*/
class InternalCoordinate : public smartptr::Serializable {

    protected:
		/** 
            Each internal coordinate must be linked to a specific molecule 
        */
        ConstMoleculePtr mol_;

		/** 
            The B matrix associated with the internal coordinate.  This is NOT the complete B matrix
		    as given in Wilson, Decius, and Cross. Rows of bmatrix_ are the atoms and the columns are
			the x,y, and z components of the B vectors. In particular, the matrix is given by
            \f[
                B_{nj} = \frac{ds}{dx_{nj}}
            \f]
            where \f$s\f$ is the internal coordinate value and \f$x_{nj}\f$ is the the \f$j\f$th coordinate
            of the \f$n\f$th atom.
        */
		RectMatrixPtr bmatrix_;

		/** 
            A descriptive std::string for the type of coordinate. See @link<type> type() @endlink
        */
        std::string coordtype_;

		/** 
            A std::vector giving the characters of the coordinate under the various symmetry operations 
		    of the molecule.  The symmetry operations are ordered as in Cotton, 3rd ed for default Abelian
            symmetry.  For non-Abelian, the user inputs the order in the input file. */
        VectorPtr chars_;

        std::map<std::string, VectorPtr> subgroup_chars_;

        std::vector<ConstInternalCoordinatePtr> equiv_coords_;

    protected:
		/** 
            The initialize function computes the B matrix at the given molecular geometry and computes
		    the characters of the coordinate */
        void init();

        /** 
            Returns the derivative of the coordinate with respect to coordinate displacement.
            This is done by numerical differentiation.  This is typically use for validating
            complicated analytical formulas
            @param n The number of the atom to displace in the derivative
            @param coord The coordinate number (x=0,y=1,z=2) to displace
            @param disp The size of the coordinate displacement.  This should be chosen to maximize
                        accuracy by balancing machine precision and small displacement size.
            @return The numerical B vector derivative
        */
        double getUnitVectorDerivative(int n, int coord, double disp);

    public:

		/** 
            Parent constructor.  At this level, InternalCoordinate is still abstract and cannot
            be instantiated.
            @param mol The molecule the coordinate is based on
            @param name The name of the coordinate
        */
        InternalCoordinate(const ConstMoleculePtr& mol, gigstr name);

        InternalCoordinate(const ArchivePtr& arch);

        virtual ~InternalCoordinate();

        void serialize(const ArchivePtr& arch) const;

        /** 
            Static method for constructing delocalized internal coordinates from
            a set of simple internal coordinates.
            @param mol The molecule the coordinates are based on
            @param simples A vector of simple internal coordinates to use in constructing delocalized internals
            @param coords A vector of coordinates.  The delocalized internals will be appened to this vector upon return.
            @param tol The tolerance for throwing out small eigenvalues. Barring machine precision issues, the matrix
                       should be positive definite.
            @return The eigenvalues of the B*B^T matrix 
        */
        static VectorPtr
        addDelocalizedInternalCoordinates(
            const ConstMoleculePtr& mol,
            const smartptr::Set<ConstSimpleInternalCoordinatePtr>& simples,
            std::vector<InternalCoordinatePtr>& coords,
            double tol = 1e-8
        );

        static VectorPtr
        addDelocalizedInternalCoordinates(
            const ConstMoleculePtr& mol,
            const smartptr::Set<ConstSimpleInternalCoordinatePtr>& simples,
            std::vector<SymmetryInternalCoordinatePtr>& coords,
            double tol = 1e-8
        );

        /** 
            Static method for checking the numer of nonredundant internal coordinates, i.e.
            the completeness of the chosen set. This method will print a status summary
            to stdout.
            @param coords The set of internal coordinates to check completness of
            @param mol  The molecule the coordinates are based on
        */
        static VectorPtr
        checkCompleteness(
            const smartptr::Set<ConstInternalCoordinatePtr>& coords, 
            const ConstMoleculePtr& mol
        );


        /** 
            Compute the characters for the given internal coordinate 
            @throw ProgrammingError if called before molecular point group is computed and closed.
        */
        void computeCharacters();

        /**
            Test equivalence of internal coordinates.  This loops through all the symmetry elements
            of the molecule and does a numeric check to see if the one of the B vectors gets 
            mapped onto the other
            @param coord The internal coordinate to test equivalence to
            @param tol The tolerance for declaring two B vectors equivalent
            @return Whether the coordinates are equivalent
        */
        bool isEquivalent(const ConstInternalCoordinatePtr& coord) const;

        void addDegenerateCoordinate(const ConstInternalCoordinatePtr& coord);

        std::vector<ConstInternalCoordinatePtr> getDegenerateCoordinates() const {return equiv_coords_;}

        /**
            Accessor method for the molecule the coordinate is based on
            @return The molecule for this coordinate 
        */
        ConstMoleculePtr mol() const;

        /**
            Maps this coordinate's B matrix via the given symmetry operation
            @param op The symmetry operation
            @return The B matrix after the symmetry operation
        */
        RectMatrixPtr mapBMatrix(const ConstSymmetryOperationPtr& op) const;

        /**
            Print details about the internal coordinate (coordinate type, characters, value)
            @param os The output stream to print to (e.g. cout, stringstream, cerr)
        */
        virtual void printDetail(std::ostream& os = std::cout) const = 0;

        /**
        */
        static int
        getCoordinateNumber(
            const ConstInternalCoordinatePtr& coord,
            const smartptr::Set<ConstInternalCoordinatePtr>& coords
        );

        /**
        */
        static RectMatrixPtr formBMatrix(
            const smartptr::Set<ConstInternalCoordinatePtr>& coords
        );

		/** 
            Returns the B matrix associated with the molecule
			@return A matrix with rows corresponding to atoms and columns corresponding to the x,y,z disps of 
			        their B std::vectors.
		*/
        ConstRectMatrixPtr getBMatrix() const;

		/** 
            Returns the B std::vector associated with an individual atom 
			@param A RefPtr to an atom
			@return The B std::vector for the internal coordinate associated with that atom
		*/
        ConstVectorPtr getBVector(const ConstAtomPtr& atom) const;

		/** 
            Get the value of the coordinate at the current molecular geometry. In bohr/angstrom or radian depending on the coordinate.
            
            @return Internal coordinate value
		*/
        virtual double getValue() const = 0;

        /**
            Get a "human readable" value for the internal coordinate. This returns angstrom for lengths
            and degrees for angles.

            @return Value of internal coordinate in conventional units
        */
        virtual double getConvertedValue() const;

        /**
            Assuming molecules are equivalent except for geometries being displaced, computes a value
            for this internal coordinate based on a different molecule. The underlying #mol_ is not 
            changed upon exit.  This currently does not do any sanity checks to make the molecules
            are equivalent.

            @param mol The new molecule to get the value for
            @return The value of the internal coordinate in the other molecule
        */
        virtual double getValueForMolecule(const ConstMoleculePtr& mol) const = 0;

        /**
            Determines whether a value is valid.  Used for sanity checks on things like
            torsions and angles which have well-defined ranges.  
            By default, the value is assumed valid unless overwritten.

            @param val The value to check
            @return Whether the value is valid
        */
        virtual bool isValidValue(double val) const;

        /** 
            For torsion angles, etc, that are periodic, we need a way of canonicalizng the values
            to recognize that 2\f$\pi\f$ and 0 are the same, for example. By default, for things like regular angles, bonds, linx,
		    we just send back the original value.
		    @param  val The value to canonicalize
		    @return The canonicalized value, adjusted to meet a given choice of periodicity
		*/
        virtual double canonicalizeValue(double val) const;

        /** 
            Recomputes the B vector for the coordinate. This is typically called
            after changing the molecular coordinates. No "extended" recompute is done for symmetry internal
            coordinates. If the underlying simples are not recomputed, the symmetry internal
            coordinate will simply recompute to the exact same value.
        */
        virtual void recompute() = 0;

		/** 
            Returns the characters of the coordinates with symm ops ordered as in Cotton, 3rd ed 
            for default Abelian groups. For non-abelian, order is given by user input.

			@return The vector of characters for the internal coordinate
		*/
        ConstVectorPtr characters() const;

        double character(int cls) const;
        
        ConstVectorPtr subgroup_characters(gigstr subgroup) const;


        /**
            Returns a description of the connectivity for the coordinate, e.g. R(1-2) for bond length,
            A(1-2-3) for bond angle, etc.

            @return A string describing the coordinate connectivity
        */
        virtual std::string connectivityString() const = 0;

		/** Prints a description of the coordinate to standard out 
            @param os The ostream to print to, e.g. cout, cerr, stringstream
			@param includechar Whether to also print the symmetry operation characters of the coordinate
            @param includesubchars Wether to also print the characters for subgroups
		*/
        void print(std::ostream& os = std::cout, bool includechars = true, bool includesubchars = false) const;

		/** 
            Returns a description of the coordinate. <br>
           <b><i> STRE </i></b> bond length <br>
           <b><i> BEND </i></b> bond angle <br>
           <b><i> OUT  </i></b> out of plane bend <br>
           <b><i> TORS </i></b> torsion <br>
           <b><i> LINX </i></b> x linear bend <br>
           <b><i> LINY </i></b> y line bend <br>
           <b><i> LIN1 </i></b> fixed vector linear bend <br>
           <b><i> SYMM </i></b> symmetry internal coordinate <br>

			@return The coordinate type
		*/
        std::string type() const;

        /**
            Debug method for validating analytic formulas for new coordinates.  This computes
            B vector derivatives through numerical differentation and compares to the analytic
            result. Aborts on failed test.

            @param tol The tolerance for numeric and analytic results being equal
        */
        void testBVectors(double tol = 1e-10);

        virtual InternalCoordinatePtr copy(const ConstMoleculePtr& mol, const smartptr::Set<ConstSimpleInternalCoordinatePtr>& simples) const = 0;

        virtual InternalCoordinatePtr copy() const = 0;

        virtual void
        intder_normalize() = 0;


}; //end InternalCoordinate class declaration


/** 
    @class SimpleInternalCoordinate    

    Enapsulates a simple internal coordinate such as a bond length or angle 
*/
class SimpleInternalCoordinate : public InternalCoordinate
{
    protected:
		/** 
            Maps coordinate specific numbering to atom object. See #getAtom.
        */
        std::map<int, ConstAtomPtr> atom_map_;

        /**
            In methane, the bond angle H2-C1-H5 would have coordinate specific numbering
            1->H2, 2->C1, 3->H5 so connectivity vector would be <2,1,5>.
        */
        std::vector<int> connectivity_;

    protected:
		/** 
            Given the coordinate specific numbering, returns the atom number in the molecule.
            For example, in methane, the bond angle H2-C1-H5 would have coordinate specific numbering
            1->H2, 2->C1, 3->H5 so mol_number would map 1->2, 2-1, 3->5.

            @param num The coordinate specific number
            @return The atom number in the molecule
        */
        int mol_number(int num) const;

    private:
        /**
            If the molecule changes, reset atoms should be called to rebuild #atom_map_.
        */
        void resetAtoms();

    public:
		/** 
            Constructor
			@param connect	The connectivity of the internal coordinate. See #connectivity_
			@param mol		The molecule to build the internal coordinates on
            @param name     The name of the coordinate type
		*/
        SimpleInternalCoordinate(
            const std::vector<int>& connect, 
            const ConstMoleculePtr& mol,
            gigstr name
        );

        SimpleInternalCoordinate(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

        /**
            Returns a string representing the connectivity. See #connectivity_.
            @return Connectivity string
        */
        std::string connectivityString() const;

        /**
            See InternalCoordinate::printDetail
        */
        void printDetail(std::ostream& os = std::cout) const;

        void intder_normalize();

        /**
            Given a coordinate number relative to the coordinate, returns the corresponding atom.
            The numbering matches that of Wilson, Decius, and Cross. For example, in methane, the bond
            angle H2-C1-H5 would have coordinate specific numbering 1->H2, 2->C1, 3->H5.
            @param idx The atom number in the coordinate
            @return The atom corresponding to the number
        */
        ConstAtomPtr getAtom(int idx) const;

        double getValueForMolecule(const ConstMoleculePtr& mol) const;

        /**
            Fetches the dummy atom numbers for this molecule.  The dummy atom number is zero-based
            corresponding to the molecule numbering. By default, nothing is done.

            @param dummy_list The vector to append dummy atom numbers to.  Numbers are appended
                              and no content is cleared from the vector.
        */
        virtual void dummies(std::vector<int>& dummy_list) const;

        virtual bool matches(const ConstSimpleInternalCoordinatePtr& coord) const;

        /**
            Fetches the connectivity of the internal coordinate. See #connectivity_.

            @param connect The vector to put the connectivity into. Any content in the vector
                           is cleared upon entry to the function.
        */
        void connectivity(std::vector<int>& connect) const;

        virtual SimpleInternalCoordinatePtr simple_copy(const ConstMoleculePtr& mol) const = 0;

        InternalCoordinatePtr copy(const ConstMoleculePtr& mol, const smartptr::Set<ConstSimpleInternalCoordinatePtr>& simples) const;

        InternalCoordinatePtr copy() const; 

};

/** 
    @class BondLength

    Encapsulates a bond length 
*/
class BondLength : public SimpleInternalCoordinate {
    
    public:
        /**
            See SimpleInternalCoordinate::SimpleInternalCoordinate
        */
        BondLength(
            const std::vector<int>& connect, 
            const ConstMoleculePtr& mol
        );

        BondLength(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

        /**
        */
        double getValue() const;

        /**
        */
        void recompute();

        SimpleInternalCoordinatePtr simple_copy(const ConstMoleculePtr& mol) const;

        /**
            Static method for computing the distance between two atoms.
            @param at1 atom 1
            @param at2 atom 2
            @return The distance between atom 1 and atom 2
        */
        static double compute(const ConstAtomPtr& at1, const ConstAtomPtr& at2);

};

/** Encapsulates a bond angle */
class BondAngle: public SimpleInternalCoordinate {
    
    public:
        /**
            See SimpleInternalCoordinate::SimpleInternalCoordinate
        */
        BondAngle(
            const std::vector<int>& connect, 
            const ConstMoleculePtr& mol
        );

        BondAngle(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

        SimpleInternalCoordinatePtr simple_copy(const ConstMoleculePtr& mol) const;

        /**
        */
        double getValue() const;

        /**
        */
        void recompute();

        /**
            Override parent method. Checks if the value is between
            0 and \f$2\pi\f$.
            @param val The value to check
            @return If the value is valid
        */
        bool isValidValue(double val) const;

        /**
            Override parent method.
            @return The bond angle in degrees instead of radians.
        */
        double getConvertedValue() const;

        /**
            Static method for computing the angle between two bond vectors.
            @param at1 atom 1
            @param at2 atom 2
            @param at3 atom 3
            @return The bond angle between R(1-2) and R(1-3)
        */
        static double compute(const ConstAtomPtr& at1, const ConstAtomPtr& at2, const ConstAtomPtr& at3);
};



/** Encapsulates an out of plane bend.  Numbering follows Wilson, Decius, and Cross */
class OutOfPlaneBend: public SimpleInternalCoordinate {

    public:
		/** 
            See SimpleInternalCoordinate::SimpleInternalCoordinate
		*/
        OutOfPlaneBend(
            const std::vector<int>& connect, 
            const ConstMoleculePtr& mol
        );

        OutOfPlaneBend(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

        SimpleInternalCoordinatePtr simple_copy(const ConstMoleculePtr& mol) const;

        /**
        */
        double getValue() const;

        /**
            Override parent method. Checks if the value is between
            \f$\frac{-\pi}{2}\f$ and \f$\frac{\pi}{2}\f$
            @param val The value to check
            @return If the value is valid
        */
        bool isValidValue(double val) const;

        /**
        */
        void recompute();
};

/**
    @class PeriodicCoordinate

    Encapsulates a coordinate whose values are constrained to lie within a given periodic range.
*/
class PeriodicCoordinate : public SimpleInternalCoordinate {

        protected:
            /**
                The point at which the period of the torsion should reset. For example, \f$2\pi\f$ and 0 
                are the same value for torsions. The discontinuity might be set to \f$\frac{-\pi}{4}\f$
                so that all torsion values are computed as \f$\frac{-\pi}{4} < x < \frac{7\pi}{4}\f$. 
            */
            double discontinuity_;

            /**
                The period of the value
            */
            double period_;

        public:

            /**
                See SimpleInternalCoordinate::SimpleInternalCoordinate
                
                @param connect
                @param mol
                @param discontinuity See #discontinuity_
                @param period I hope you get this
                @param name
            */
            PeriodicCoordinate(
                const std::vector<int>& connect, 
                const ConstMoleculePtr& mol, 
                double discontinuity, 
                double period,
                gigstr name);

            PeriodicCoordinate(const ArchivePtr& arch);

            void serialize(const ArchivePtr& arch) const;

            /** 
                Take value and adapts it to the periodicity condition. See #discontinuity_.  The value is computed as
                \f[
                    y = (x - d) \bmod p + d
                \f]
                where x is the original value, d is the discontinuity, and p is the period.
                @param The non-periodic value
                @return The periodic value
            */
            double canonicalizeValue(double val) const;
};

/** 
    @class Torsion

    Encapsulates a torsion 
*/
class Torsion: public PeriodicCoordinate {

    public:
        /**
            See PeriodicCoordinate::PeriodicCoordinate
		*/
        Torsion(
            const std::vector<int>& connect, 
            const ConstMoleculePtr& mol, 
            double discontinuity = -PI/4
        );

        Torsion(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

        SimpleInternalCoordinatePtr simple_copy(const ConstMoleculePtr& mol) const;

        /**
        */
        double getValue() const;

        /**
        */
        void recompute();

        /**
            Override parent method.
            @return The bond angle in degrees instead of radians.
        */
        double getConvertedValue() const;

};

/** 
    @class LinX

    Encapsulates a LinX motion. Atoms 1,2,3 define an XZ plane, atoms 1,2 define the Z-axis and LinX is the X value
    on these axes.  It describes an in-plane motion when 3 atoms are nearly linear.
*/
class LinX: public SimpleInternalCoordinate {
    
    private:
        /**
            See LinX::LinX
        */
        BondAnglePtr angle_;

        /**
            See LinX::LinX
        */
        TorsionPtr torsion_;

    public:
		/** 
            Constructor. The constructor currently works by considering the value to be \f$ s = \cos \tau \sin \theta \f$
                         for the torsion \f$\tau(1,2,3,4)\f$ and the angle \f$\theta(1,2,3)\f$
                         
			@param connect The atom numbers defining the connectivity pattern.  Atoms 1,2,3 define an XZ plane,
                            atoms 1,2 define the Z-axis and LinX is the X value on these axes.
			@param mol		The molecule to build the coordinate off of
		*/
        LinX(
            const std::vector<int>& connect, 
            const ConstMoleculePtr& mol
        );

        LinX(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

        SimpleInternalCoordinatePtr simple_copy(const ConstMoleculePtr& mol) const;

        /**
        */
        double getValue() const;

        /**
        */
        void recompute();
};

/** 
    @class LinY

    Encapsulates a LinY motion. Atoms 1,2,3 define an XZ plane, atoms 1,2 define the Z-axis and LinY is the Y value
    on these axes describing an out-of-plane motion 
*/  
class LinY: public SimpleInternalCoordinate {
    
    private:
        /**
            See LinY::LinY
        */
        BondAnglePtr angle_;

        /**
            See LinY::LinY
        */
        TorsionPtr torsion_;

    public:
		/** 
            Constructor. The constructor currently works by considering the value to be \f$ s = \sin \tau \sin \theta \f$
                         for the torsion \f$\tau(1,2,3,4)\f$ and the angle \f$\theta(1,2,3)\f$
			@param connect The atom numbers defining the connectivity pattern.  Atoms 1,2,3 define the plane that Atom 4
							is shifted out of.
			@param mol		The molecule to base the coordinate on
		*/
        LinY(
            const std::vector<int>& connect, 
            const ConstMoleculePtr& mol
        );

        LinY(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

        SimpleInternalCoordinatePtr simple_copy(const ConstMoleculePtr& mol) const;

        /**
        */
        double getValue() const;

        /**
        */
        void recompute();
};


/**
    @class Lin1

    Encapsulates a Lin1 motion.  Atomss 1,2,3 define an XZ plane.  The dummy defines a direction perpendicular to the plane.
    The value of the coordinate is the displacement out of plane alone the perpendicular direction.
*/
class Lin1 : public SimpleInternalCoordinate {

    private:
        /**
            The unit vector defining the direction of the Lin1 displacement
        */
        VectorPtr ed_;

        /**
            Each molecule must keep track of a dummy atom representing the Lin1 vector.
            After the coordinate registers the dummy atom with the molecule, the molecule
            tells the coordinate what atom number has been assigned to the dummy atom.
        */
        int dummy_;

    public:
		/** 
            Constructor
			@param connect The atom numbers defining the connectivity pattern.  Atoms 1,2,3 define the plane that Atom 4
							is shifted out of.
			@param mol		The molecule to build the coordinate off of
		*/
        Lin1(
            const std::vector<int>& connect, 
            const ConstMoleculePtr& mol, 
            const std::vector<double>& vec
        );

        Lin1(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

        SimpleInternalCoordinatePtr simple_copy(const ConstMoleculePtr& mol) const;

        /**
        */
        double getValue() const;

        bool matches(const ConstSimpleInternalCoordinatePtr& coord) const;

        /**
        */
        void recompute();

        /**
        */
        ConstVectorPtr ed() const;

        /**
        */
        void dummies(std::vector<int>& dummy_list) const; 
};

/** 
    @class SymmetryInternalCoordinate

    Encapsulates a symmetry internal coordinate that is a linear combination of simple internals 
*/
class SymmetryInternalCoordinate : public InternalCoordinate
{
    private:
        /**
        */
        VectorPtr coeffs_;

        /**
        */
        smartptr::Set<ConstSimpleInternalCoordinatePtr> simples_;

    public:
        SymmetryInternalCoordinate(
            const ConstVectorPtr& coeffs, 
            const smartptr::Set<ConstSimpleInternalCoordinatePtr>& coords, 
            const ConstMoleculePtr& mol);

        SymmetryInternalCoordinate(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

        double getValueForMolecule(const ConstMoleculePtr& mol) const;

		/** Gets the coefficients definining the linear combination of simple internal coordinates. This method
			only returns the coefficients for the simple internals included in the symmetry internal 
			@return The coefficients of the symmetry internal
		*/
        ConstVectorPtr getCoefficients() const {return coeffs_;}

        /** Convert is ignored since computing things in terms of angstroms and degrees would give basically meaningless
			numbers since degrees would dominate the value */
        double getValue() const;

        /**
        */
        void recompute();

        void normalize();

        /**
        */
        void printDetail(std::ostream& os = std::cout) const;
	
		/** 
            Gets the coefficients defining the linear combination of simple internal coordinates.  This method
			returns the coefficients for all simple internals. The method basically just fills in zeros
			for all the simple internals that do not contribute.
			@param coords	The complete set of simple internal coordinates.  This must be given since
							the symmetry internal only keeps track of the contribution simples.
			@return The std::vector of coefficients
		*/
		VectorPtr getCoefficients(const smartptr::Set<ConstSimpleInternalCoordinatePtr>& coords) const;

		/** Translates the number of internal coordinates from the numbering for the given symmetry internal coordinate
			(which may only depend on 3 internals) to the numbering for the complete set of simple internals (which may be 6,9,12)
			@param simplenum The number (0,1,2...) of the simple coordinate in the linear combination
			@param coords	 The complete set of internal coordinates
			@return The number of the simple internal in the complete set based on zero counting
		*/
        int getSimpleCoordinateNumber(int simplenum, const smartptr::Set<ConstSimpleInternalCoordinatePtr>& coords) const;

        /**
            See @link<connectivityString> connectivityString() @endlink
        */
        std::string connectivityString() const;

        smartptr::Set<ConstSimpleInternalCoordinatePtr> getSimples() const;

        //InternalCoordinatePtr copy() const; 

        InternalCoordinatePtr copy() const;

        SymmetryInternalCoordinatePtr carbon_copy() const;

        InternalCoordinatePtr copy(const ConstMoleculePtr& mol, const smartptr::Set<ConstSimpleInternalCoordinatePtr>& simples) const;

        void intder_normalize();

};

class SymmetryOperation;
/**
    @class CoordinateSubspace

    Encapsulate a subspace containing a set of coordinates, usually used for computing properties of a set
    of degenerate coordinates, e.g. T1 irrep of tetrahedral.
*/
class CoordinateSubspace : public smartptr::Serializable {

    private:
        /**
            The coordinates in the subspace
        */
        std::vector<SymmetryInternalCoordinatePtr> coords_;

        /**
            The molecule the coordinates are based on
        */
        ConstMoleculePtr mol_;

        /**
            Flag indicating whether or not symmetry internals have been adapted
            to an abelian subgroup
        */
        bool madeabelian_;

    public:
        /**
            Initializes an empty coordinate subspace.

            @param mol The molecule the coordinate subspace is built on
        */
        CoordinateSubspace(const ConstMoleculePtr& mol);

        CoordinateSubspace(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

        /**
            Order = 1 for abelian irreps <br>
            Order = 2 for E irreps <br>
            Order = 3 for T irreps, etc... <br>
            @return The number of coordinates in the subspace
        */
        int order() const;

        /**
            Adds a coordinate to the subspace.  Currently no sanity checks are done to
            make sure the coordinate fits the subspace.

            @param coord The coordinate to add
        */
        void addCoordinate(
            const SymmetryInternalCoordinatePtr& coord
        );

        void print(std::ostream& os = std::cout) const;

        /**
            Fetches the ith coordinate in the subspace. This numbering is relative to
            the subspace, not the complete set of coordinates.

            @param i The coordinate number to fetch
            @return The ith coordinate
            @throw GigideException if coordinate number is greater the order() - 1
        */
        ConstSymmetryInternalCoordinatePtr getCoordinate(int i) const;

        /**
            Fetches the ith coordinate in the subspace. This numbering is relative to
            the subspace, not the complete set of coordinates.

            @param i The coordinate number to fetch
            @return The ith coordinate
            @throw GigideException if coordinate number is greater the order() - 1
        */
        SymmetryInternalCoordinatePtr getCoordinate(int i);

        /**
            Computes the character of a symmetry operation.  The characters are computed
            by choosing the coordinate number 0 as a representative for the subspace.

            @param oper The symmetry operation to compute the character of
            @return The character
            @throw GigideException Throw if subspace is empty

        */
        double
        character(
            const ConstSymmetryOperationPtr& oper
        ) const;

        /**
            Computes and prints the characters of a set of symmetry operations.
            @param opers The set of operations to compute characters for
            @param os The ostream to print to
            @throw GigideException Throws if subsapce is empty
        */
        void printCharacters(
            const smartptr::Set<ConstSymmetryOperationPtr>& opers,
            std::ostream& os = std::cout
        ) const;

        /**
            Computes the projection of a vector onto the coordinate subspace.
            This assume orthogonality of the vectors in the subspace.  This will give 
            an incorrect answer due to non-orthogonal overcounting otherwise. This may
            be fixed in future versions.

            @param bvec The vector to project onto
            @return The vector projection
        */
        VectorPtr
        projection(
            const ConstVectorPtr& bvec
        );

        void
        abelianify(
            const ConstPointGroupPtr& pg
        );

        
};


} //end namespace

#endif
