
#ifndef gigide_displacement_h
#define gigide_displacement_h

#include <map>
#include <vector>
#include <src/coordinates.h>
#include <src/deftypes.h>
#include <src/iter.h>
#include <src/serialize.h>

namespace gigide {

class TaylorTerm;
class Displacement;
class Derivative;

/** 
    @class DisplacementMapping

    Encapsulates a mapping between two equivalent displacements.  The contains information about
    the equivalent displacements and the symmetry operation relating them.
*/
class DisplacementMapping : public smartptr::Serializable
{
    private:
        /** 
            The displacement mapped onto by the symmetry operation 
        */
        DisplacementPtr disp_;
    
        /** 
            The 3x3 matrix describing the mapping of the X,Y,Z coordinates 
        */
        ConstRectMatrixPtr symmop_;

        /**
            The operator describing the permutation of atoms under the symmetry operation.
            See Molecule::getPermutationMatrix for the format. 
        */
        ConstRectMatrixPtr permop_;

    public:
        /**
            Constructor. 

            @param disp The displacement mapped onto by the symmetry operation
            @param symmop 3x3 matrix describing x,y,z transformation
            @param permop See Molecule::getPermutationMatrix
        */
        DisplacementMapping(
            const DisplacementPtr& disp, 
            ConstRectMatrixPtr symmop, 
            ConstRectMatrixPtr permop
        );

        DisplacementMapping(
            const ArchivePtr& parser 
        );

        void serialize(const ArchivePtr& writer) const;
        
        /**
            Get the displacement
            @return The displacement mapped onto
        */
        ConstDisplacementPtr getDisplacement() const;

        /**
            Assign an energy to the underlying displacement
        */
        void assignEnergy(double energy);

        /**
            Assign gradients to underlying displacement.  This requires a transformation
            of the gradients based on the symmetry operation.
            @param xyzgradients The xyz gradients to assign
        */
        void assignGradients(ConstRectMatrixPtr xyzgradients);

        /**
            Assign force constants to underlying displacement.  This requires a transformation
            of the force constants based on the symmetry operation.  The symmetry operation
            requires gradients and force constants to compute the correct answer.
            @param xyzfc
            @param xyzgradients
        */
        void assignForceConstants(ConstRectMatrixPtr xyzfc, ConstRectMatrixPtr xyzgrads);
};

/** 
    @class Displacement

    Encapsulates a displacement from the central molecule. 
*/
class Displacement : public smartptr::Serializable
{
    private:
        /** The list of integers defining the displacements. Suppose we have the coordinate sets
            R(O1-H2), R(O1-H3), A(H2-O1-H3). A displacement of 0.5 along R(O1-H2),
            0 along R(O1-H3), and 2.0 along A(H2-O1-H3) would be represented as
            would be <0.5, 0, 2.0> */
        std::vector<double> displacements_;

        /** The label for the displacement. The label is exactly the vector #displacements_ */
        std::string label_;

        /** The displacement matrix defining the displacement.  This is derived from the B-matrix */
        RectMatrixPtr displacement_matrix_;

        /** The std::vector of equivalent displacements.  This is included so that the energies of the equivalent
            displacements can be filled in later from this displacement's energy */
        std::vector<DisplacementMappingPtr> equiv_disps_;

        /** Whether or not an energy has been assigned yet to this displacement */
        bool energy_assigned_;
        /** Whether or not gradients have been assigned yet to this displacement */
        bool grad_assigned_;
        /** Whether or force constants have been assigned yet to this displacement */
        bool fc_assigned_;

        /** The energy at this displacement point.  Units here should match the units in the input file */
        double energy_;

        /** The gradients at this displacement point, if gradients are computed */
        ConstVectorPtr gradients_;

        /** The force constants at this displacement point, if force constants are computed */
        ConstRectMatrixPtr fc_;

        /** The center molecule of the force field */
        ConstMoleculePtr mol_;

        /** The molecule representing the displaced geometry */
        MoleculePtr dispmol_;

        /**
            The set of simple internal coordinates used to generate the displaced geometry
        */
        Set<ConstSimpleInternalCoordinatePtr> simples_;

        /** 
            The set of coordinates the displacement is built on 
        */
        Set<ConstInternalCoordinatePtr> coords_;

        /**
            In testing a given fitting function, the displacement will be a sum
            of different polynomial terms involving derivatives.  These polynomial
            terms are encapsulated by this vector.
        */
        std::vector<TaylorTermPtr> terms_;

        /**
            Whether or not the displacement is unique or if it is 
            equivalent to another displacement.
        */
        bool unique_;

        /**
            The "reference count" of the displacement.  This is the number of 
            derivative fitting functions which require this displacement.
        */
        int refcount_;

        /**
            The degree of the displacement.  This is just the sum of the vector #displacements_.
            This variable is necessary for quick validation of the equivalence of two displacements
        */
        double degree_;

        /**
            The x,y,z magnitude of a unit displacement along the coordinate.
        */
        double dispmag_;

        /**
            A string describing the displacement type. This is essentially label_ but with the zeroes
            removed and the values sorted.  Equivalent displacements will have equivalent displacement types.
        */
        std::string disptype_;

        void init();

    public:
        /** 
            Constructor.
            @param disps See #displacements_
            @param mol The center molecule in the displacement, not the displaced geometry!
            @param coords The set of internal coordinates
            @param simples The set of simple internal coordinates the internal coordinates are based on
        */
        Displacement(const std::vector<double>& disps,
                     const ConstMoleculePtr& mol,
                     const Set<ConstInternalCoordinatePtr>& coords,
                     const Set<ConstSimpleInternalCoordinatePtr>& simples);

        
        /** 
            Constructor for the zero displacement.
            @param mol The center molecule in the displacement, not the displaced geometry!
            @param coords The set of internal coordinates
            @param simples The set of simple internal coordinates the internal coordinates are based on
        */
        Displacement(const ConstMoleculePtr& mol,
                     const Set<ConstInternalCoordinatePtr>& coords,
                     const Set<ConstSimpleInternalCoordinatePtr>& simples);

        Displacement(const ArchivePtr& parser);

        /**
            Adds xml output describing the displacement
            @param writer The xml document to add the data to
        */
        void serialize(const ArchivePtr& writer) const;

        /**
            Validate that the generated geometry matches the displacement expected
        */
        void validate() const;


        /** 
            Return displacement numbers.  See #displacements_
            @return Displacement numbers
        */
        std::vector<double> getDispNumbers() const;

        /** 
            Returns the label for the displacement. See #label_
            @return The displacement label
        */
        std::string label() const;

        /**
            Clears the set of equivalent displacements.  See #equiv_disps_.
        */
        void clearMappings();

        /**
            Returns the displacement type. See #disptype_.
            @return Displacement type
        */
        std::string disptype() const;

        /**
            Send back the degree of the displacement. See #degree_.
            @return The displacement degree
        */
        double degree() const;

        /**
            Send back the x,y,z magnitude of the displacement. See #dispmag_.
            @return Displacement magnitude.
        */
        double dispmag() const;

        /**
            Send back whether values have been assigned to the displacement
            in the form of energies, gradients, or harmonic force constants.
            @return If values have been assigned to displacement
        */
        bool isAssigned() const;

        /** Returns whether an energy has been assigned to the displacement
            @return If an energy has been assigned
        */
        bool energyAssigned() const;

        /** 
            Links this displacement to another equivalent displacement.  When an energy gets assigned
            to this displacement, all linked displacements will be assigned the same value.
            See DisplacementMapping::permop_ and DisplacementMapping::symmop_.
            @param disp    The equivalent displacement
            @param symmop  The 3x3 matrix representing the x,y,z transformation
            @param permop  The permutation matrix for the atoms.
        */
        void addEquivalentDisplacement(const DisplacementPtr& disp,
                                       ConstRectMatrixPtr symmop,
                                       ConstRectMatrixPtr permop);

        void printEquivalentDisplacements(std::ostream& os = std::cout) const;

        /**
            Lets displacement know that derivative fitting function requires it.
        */
        void incrementRefcount();

        /**
            Lets displacement know that derivative fitting function no longer requires it.
        */
        void decrementRefcount();

        /**
            Gets the reference count.
            @return #refcount_
        */
        int getRefcount() const;

        /**
            Get whether this is a unique displacement or if it is equivalent to another displacement
            via a symmetry operation.
            @return Whether the displacement is unique
        */
        bool isUnique() const;

        /**
            Sets the uniqueness of the displacement.
            @param unique Whether the displacement is unique.
        */
        void setUnique(bool uniq);

        /**
            Returns the set of internal coordinates the displacement is based on.
            @return The internal coordinates
        */
        Set<ConstInternalCoordinatePtr> coords() const;

        /** Assigns an energy to this displacement.  This will also assign the same energy to all 
            equivalent displacements.
            @param energy The energy to assign.
        */
        void assignEnergy(double energy);

        /** Assigns the gradient. This will also assign gradients to all equivalent displacements. Though
            displacements may be equivalent, their orientation may be different.  In this case, the transformation
            between them must be determined by relating their displacement matrices.
            @param gradients The gradients to assign
            @param disp The displacement matrix for the assigning displacement.  This allows the displacement
                        to translate the gradients based on orientation differences.
        */
        void assignGradients(ConstRectMatrixPtr gradients);

        /**
            Given a set of force constants and gradients, transform the force constants
            to internal coordinates and store them.
            @param xyzfc The force constants
            @param xyzgrads The xyz gradients
        */
        void assignForceConstants(ConstRectMatrixPtr xyzfc, ConstRectMatrixPtr xyzgrads);

        /**
            Reset all the Taylor terms associated with this displacement.
            See TaylorTerm::reset.
        */
        void resetTaylorTerms();

        /**
            Have all Taylor term polynomials with order <= maxlevel
            accumulate the contributions from this displacement. See TaylorTerm.
        */
        void accumulateTerms(double coef, int maxlevel);

        /**
            Put all the nonzero Taylor term contributions to this displacement
            in the map based on a string label given by TaylorTerm::name()
            @param terms The map to put the terms in
            @param tol The tolerance for declaring a term to be zero
        */
        void addNonzeroTerms(std::map<std::string, TaylorTermPtr>& terms, double tol = 1e-6) const;

        /**
            Return a string describing the Taylor series representation of the energy
            at this displacement.
            @return The taylor series string
        */
        std::string taylorRepresentation() const;

        /** 
            Gets the energy associated with the displacement.
            @return The displacement's energy
            @throw GigideException If no energy has been assigned
        */
        double getEnergy() const;

        /**
            Figure out whether this is the central displacement
            @return If this is a zero displacement
        */
        bool isZeroDisplacement() const;

        /** 
            Gets the value for a given derivative. 
            @param deriv The derivative to get a value for
            @return The derivative value
            @throws GigideException Throws if no value is found
        */
        double getDerivativeValue(const ConstDerivativePtr& deriv) const;

        /**
            Assign the set of Taylor terms to this displacement. This will filter
            to only keep track of the nonzero terms.
            @param terms The terms to link to
        */
        void
        assignTerms(
            std::vector<TaylorTermPtr>& terms
        );

        /** 
            Prints a descriptive label 
            @param os The output stream to print to
        */
        void print(std::ostream& os = std::cout) const;

        /** 
            Sends back the displacement matrix associated with the displacment. This is just the
            the sum of the B vectors for the individual coordinates scaled by the displacement size.
            See InternalCoordinate::bmatrix_.
            @return The xyz B vector displacements for this displacements.
        */
        ConstRectMatrixPtr getDisplacementMatrix() const;
        
        /**
            Returns the molecule associated with the displaced geometry.
            @return The displaced molecule
        */
        ConstMoleculePtr getDisplacementMolecule() const;

        /** 
            Generate the displaced geometry for this displacement
            @param displacements The displacement sizes for the internal coordinates
        */
        void generateDisplacement(const std::vector<double>& displacements);

        /**
            Return the displacement extent of the ith coordinate in this displacement.
            This corresponds to multiples of the displacement size.  For example, suppose
            you are computing the first derivative of the second coordinate as
            \f[
                \frac{dE}{ds_2} = \frac{E(+\delta_2) - E(-\delta_2)}{2 delta_2}
            \f]
            If this is the \f$-\delta\f$ displacement, then the function would return <br>
            displacement(0) -> 0 <br>
            displacement(1) -> 0 <br>
            displacement(2) -> -1 <br>
            @return The displacement extent
        */
        double displacement(int i) const;

        /**
            Whether this displacement matches a given set of displacements. See #displacement.
            This is equivalent to loop all elements of the vector and checking if
            disps[i] == displacement(i).
            @param disps The set of displacements
        */
        bool matches(const std::vector<double>& disps) const;

        /**
            Generate a label for the set of displacements.  This is just a string representation of the vector.
            @param disps The set of displacements
            @return string label
        */
        static std::string label(const std::vector<double>& disps);

        /** 
            Generates a molecule with the followintg displacements in the internal coordinates.
            @param displacements The displacement sizes from the current molecular geometry
            @param mol A "guess" molecule for the final geometry
            @param coords The internal coordinates
            @param simples The simple internal coordinates the internals are based on
        */
        static MoleculePtr displaceGeometry(
            const std::vector<double>& displacements,
            const ConstMoleculePtr& mol,
            const Set<ConstInternalCoordinatePtr >& coords,
            const Set<ConstSimpleInternalCoordinatePtr >& simples
        );

        /** 
            Generates a molecule with the followintg displacements in the internal coordinates.
            @param displacements The displacement sizes from the current molecular geometry
            @param mol A "guess" molecule for the final geometry
            @param coords The internal coordinates
            @param simples The simple internal coordinates the internals are based on
        */
        static MoleculePtr displaceGeometry(
            const std::vector<double>& displacements,
            const ConstMoleculePtr& mol,
            const Set<ConstInternalCoordinatePtr>& coords,
            const Set<ConstSimpleInternalCoordinatePtr >& simples,
            std::vector<InternalCoordinatePtr>& newcoords,
            std::vector<SimpleInternalCoordinatePtr>& newsimples
        );

        /** 
            Generates a geometry with the given set of values  for the internal coordinates
            @param values The internal coordinate values of the generated geometry
            @param mol A "guess" molecule for the final geometry
            @param coords The internal coordinates
            @param simples The simple internal coordinates the internals are based on
        */
        static MoleculePtr generateGeometry(
            const std::vector<double>& values, 
            const ConstMoleculePtr& mol,
            const Set<ConstInternalCoordinatePtr>& coords,
            const Set<ConstSimpleInternalCoordinatePtr>& simples
        );

        /** 
            Generates a geometry with the given set of values  for the internal coordinates
            @param values The internal coordinate values of the generated geometry
            @param mol A "guess" molecule for the final geometry
            @param coords The internal coordinates
            @param simples The simple internal coordinates the internals are based on
        */
        static MoleculePtr generateGeometry(
            const std::vector<double>& values, 
            const ConstMoleculePtr& mol,
            const Set<ConstInternalCoordinatePtr>& coords,
            const Set<ConstSimpleInternalCoordinatePtr>& simples,
            std::vector<InternalCoordinatePtr>& newcoords,
            std::vector<SimpleInternalCoordinatePtr>& newsimples
        );

        /**
            Compute how many units along a given internal coordinate the molecule is displaced.
            @param mol The molecule to compute the displacement on
            @param coord The coordinate to compute the displacement along
            @param disp  The displacement size for the coordinate.  For example, if the internal coordinate value
                         is displaced 0.02 along the coordinate and disp is 0.1, then the method would return
                         2.0 as the displacement size.
            @return The displacement size as a multiple of disp
        */
        static double computeDisplacement(const ConstMoleculePtr& mol, const ConstInternalCoordinatePtr& coord, double disp);
};

/** 
    @class DisplacementIterator 

    Iterates through displacements in a force field
*/
class DisplacementIterator : public smartptr::Serializable
{
    public:
        typedef Iterator<DisplacementPtr, ConstDisplacementPtr> const_iterator; 
        typedef std::vector<DisplacementPtr>::iterator iterator;

    private:
        /** The highest level of derivative being done */
        int deriv_level_;

        /** The number of unique displacements */
        int nunique_;

        /** The complete set of TaylorTerms necessary for every displacement in the iterator */
        std::vector<TaylorTermPtr> terms_;

        /** The complete set of displacements */
        std::vector<DisplacementPtr> displacements_; 
        
        /** The set of internal coordinates the displacements correspond to */
        Set<ConstInternalCoordinatePtr> coords_;

        /** The set of simple internal coordinates the displacements correspond to */
        Set<ConstSimpleInternalCoordinatePtr> simples_;

        /** The molecule the displacements will be performed on */
        ConstMoleculePtr mol_;

        /**
            Read a set of displacements from a file and generate all possible displacement combinations
            of the the set of internal coordinates.  The file should be formatted as, for example <br><br>
            dispfile.txt <br>
            1 <br>
            -1 <br>
            1 1 <br>
            -1 -1 <br><br>

            For every possible combination of internal coordinates, this would generate all combinations with 
            a displacement of +1 in the coordinates, etc.  For water with 3 internal coordinates, the following
            displacements would be generated. <br><br>
            <1,0,0>; <0,1,0>; <0,0,1> <br>
            <-1,0,0>; <0,-1,0>; <0,0,-1> <br>
            <1,1,0>; <1,0,1>; <0,1,1> <br>
            <-1,-1,0>; <-1,0,-1>; <0,-1,-1> <br><br>

            @param filename The file to read the displacements form
            @param disp_vec The vector to the displacement in
        */
        void
        readDisplacements(
            const char* filename,
            std::vector<std::vector<double> >& disp_vec
        );

    public:
        /** Constructor
            @param coords The set of internal coordinates to build displacements for
            @param simples The set of simples internals the coords are based on
            @param mol    The molecule displacements will be built from
        */
        DisplacementIterator(
            const Set<ConstInternalCoordinatePtr>& coords, 
            const Set<ConstSimpleInternalCoordinatePtr>& simples, 
            const ConstMoleculePtr& mol
        );

        DisplacementIterator(const ArchivePtr& parser);

        void serialize(const ArchivePtr& writer) const;

        /** 
            Tries to find an equivalent displacement from the complete set of displacements
            @param disp The displacement to check
            @return If an equivalent displacement has been found
        */
        bool mapEquivalentDisplacements(const DisplacementPtr& disp);

        void mapEquivalentDisplacements();

        ~DisplacementIterator();

        /**
            Tests whether two displacements are equivalent by a set of symmetry operations
            @param refdisp The reference displacement
            @param disp The displacement to check equivalence to refdisp
            @param symmops The set of symmetry operations to use in checking equivalence
            @param permop If a symmetry operation is found, the reference is updated to the correct permutation
                           See DisplacementMapping::permop_;
            @param symmop If a symmetry operation is found, the reference is updated to the correct 3x3 matrix
                           See DisplacementMapping::symmop_;
        */
        bool testEquivalence(
            const ConstDisplacementPtr& refdisp,
            const ConstDisplacementPtr& disp,
            const Set<ConstSymmetryOperationPtr>& symmops,
            RectMatrixPtr& permop,
            RectMatrixPtr& symmop
        );


        /**
            Read and generate displacements from a dispfile. See readDisplacements.
            @param dispfile The file to read
        */
        void
        readDispFile(std::string dispfile);

        /**
            This has the side effect of assign Taylor terms to the displacement.

            @param disp The displacement to add to the iterator.
        */
        void addDisplacement(DisplacementPtr disp);

        /**
            Get a displacement matching the given displacement sizes. See #displacements_.
            If the displacement is not found, it is created and added to the iterator.
            Contrast with findDisplacement.

            @param disps The displacements sizes of the displacement to find
            @return The matching displacement
        */
        DisplacementPtr getDisplacement(std::vector<double>& disps);

        /**
            Find a displacement matching the given displacement sizes. See #displacements_.
            If the displacement is not found, it returns null. Contrast with getDisplacement.

            @param disps The displacements sizes of the displacement to find
            @return The matching displacement
        */
        DisplacementPtr findDisplacement(const std::vector<double>& disps);

        /**
            Compute the displacement of a given set of xyz coordinates from the central molecule
            as multiples of a given set of displacement sizes.  For example if the internal coordinate
            displacements are <0.02,0.0,0.02> and the displacement sizes are <0.01,0.01,0.02> then
            the computed displacements array would be <2.0, 0.0, 1.0>

            @param xyz 
        */
        std::vector<double> computeDisplacements(ConstRectMatrixPtr xyz, const std::vector<double>& dispsizes) const;

        /**
            Looks for the zero displacement.  Equivalent to calling findDisplacement with 
            an array of zeroes.
            @return The zero displacement if found.  Otherwise null.
        */
        DisplacementPtr findZeroDisplacement();

        /**
            @return If the iterator contains a zero displacement
        */
        bool hasZeroDisplacement();

        void addZeroDisplacement();

        /** 
            Starts the iteration through all the displacements 
        */
        const_iterator begin() const;

        iterator begin();

        /** 
            @return Whether the iteration is done
        */
        const_iterator end() const;

        iterator end();

        /** 
            @return The number of internal coordinates the displacement are based on
        */
        int ncoord() const;

        /** 
            @return The number of all displacements (not just unique ones) 
        */
        int ndisps() const;
    
        /**
            @return The number of unique displacements based on symmetry considerations
        */
        int nunique() const;

        /** 
            @return The level of derivative being computed
        */
        int level() const;

        void print(std::ostream& os = std::cout) const;

};


}

#endif
