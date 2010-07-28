#ifndef _intdif_molecule_h_
#define _intdif_molecule_h_

#include <vector>
#include <string>
#include <map>
#include <src/defines.h>
#include <src/keyword.h>
#include <src/coordinates.h>
#include <src/printstream.h>
#include <src/deftypes.h>
#include <src/serialize.h>
#include <iostream>

namespace gigide {

class SymmetryOperation;
class PointGroup;
class PointGroupClass;

/** 
    @class Atom
    Encapsulates an atom in a molecule 
*/
class Atom : public smartptr::Serializable {

    private:
		/** The mass of the atom */
        double  mass_;

		/** The charge of the atom */
        int charge_;

		/** The atomic symbol */
        std::string  symbol_;

		/** The vector defining the xyz coordinates */
        VectorPtr xyz_;

    public:
		/** 
            Constructor
			@param xyz A std::vector given the xyz coordinates
			@param atomic_symbol A symbol specifying the atom.  For now, only the standard H,He,C,etc. symbols are allowed
								 In future versions, non-standard D,C13,isotope symbols will be allowed.
		*/
        Atom(ConstVectorPtr xyz, gigstr atomic_symbol);

        Atom(const ArchivePtr& arch);

        void serialize(const ArchivePtr& arch) const;

		/** 
            Gets the xyz coordinates
			@return The std::vector of xyz coordinates
		*/
        ConstVectorPtr getXYZ() const;

		/** 
            Sets new xyz coordinates.  The assignment works by replacing the values in the
			original vector.  
			@param newxyz The coordinates to assign
		*/
        void setXYZ(ConstVectorPtr newxyz);

		/** 
            Scales the xyz coordinates by a certain value
			@param scfactor The scale factor for the scaling operation
		*/
        void scaleXYZ(double scfactor);

		/** 
			@return The atomic symbol as a string
		*/
        std::string symbol() const;

		/** 
			@return The mass of the atom
		*/
        double mass() const;

        void print(std::ostream& os = std::cout) const;

        /**
            @param disp The extent of the xyz displacement
        */
        void displaceXYZ(ConstVectorPtr disp);
};


/** 
    Enapsulates a Molecule 

    @class Molecule
*/
class Molecule : public smartptr::Serializable {

    private:
		/** 
            The std::vector of atoms in the molecule 
        */
        std::vector<AtomPtr> atoms_;

		/** 
            The xyz coordinates. This will be in in whatever units are specified by the
            bond units keyword. 
        */
        RectMatrixPtr xyz_;


        /**
            A mapping between atom pointer and the number of the atom.
            Zero based counting.
        */
        std::map<const Atom*, int> number_map_;

        /**
        */
        PointGroupPtr pg_;

        /**
            The set of simple internal coordinates for the molecule
        */
        std::vector<SimpleInternalCoordinatePtr> simples_;

        /* 
            Vector of dummy atoms 
        */
        std::vector<AtomPtr> dummies_;

        /**
            Whether or not the molecule is linear
        */
        bool islinear_;

	private:
		/** 
            Recenters the molecule to a center of mass origin 
        */
        void recenter();

		/** 
            Aligns the molecule on its principal axex 
        */
        void reorient();

		/** 
            Sets the xyz coordinates of the atoms after an update of the molecule xyz
			coordinates.  This must happen currently because the molecule and atom
			xyz coordinates are stored in two different locations.  This needs
			to change
		*/
        void resetAtomXYZ();

        /**
        */
        static bool statics_done;

        /**
            Whether or not the point group absolutely needs to be computed
        */
        static bool needpg_;

        /**
        */
        static void initialize_statics();

		/** 
            A map relating atomic symbols to atomic masses 
        */
        static std::map<std::string, double> MASSES;


        /**
        */
        void init();

    public:
		/** 
            Constructor
			@param xyz The matrix of xyz coordinates
			@param atomlist	A vector defining the atom types.  This must be in the same order as xyz
		*/
        Molecule(
            ConstRectMatrixPtr xyz, 
            const std::vector<std::string>& atomlist
        );

        /**
            smartptr::Serializable constructor
        */
        Molecule(
            const ArchivePtr& arch
        );

        Molecule(const ConstMoleculePtr& mol);

        void serialize(const ArchivePtr& arch) const;

        /**
            Compute the abelian point group corresponding to the molecule
        */
        void computeAbelianSymmetry();

        /**
            Compute the full point group of the molecule based on symmetry
            operations provided by the user in the input file.
        */
        void computePointGroup();

        void computePointGroup(const ArchivePtr& arch);

		/** 
            Figures out whether a molecule has the given symmetry operation
			@param op	A matrix specifying some symmetry operation
		*/
        bool hasSymmOp(ConstSymmetryOperationPtr op) const;

		/** 
            Given a RefPtr to an atom, figures out the number of the atom in the xyz coordinates
			based on zero counting
			@param atom A RefPtr to an atom
			@return	The atom number in the molecule
		*/
        int findAtomNumber(const ConstAtomPtr& atom) const;

        void print(std::ostream& os = std::cout) const;

		/** 
            Sets the xyz coordinates of the molecule. This assigns the values, not the matrix pointer.
			@param newxyz The xyz coordinates to assign
		*/
        void setXYZ(ConstRectMatrixPtr newxyz);

        /**
            @param op The symmetry operation to add
        */
        void addSymmetryOperation(const ConstSymmetryOperationPtr& op);

        /**
            Get the symmetry elements for the molecule

            @param symmops Vector to hold symmetry operations
        */
        void getSymmetryElements(std::vector<ConstSymmetryOperationPtr>& symmops) const;

		/** 
            Scales the xyz coordinates by a given value
			@param scfactor The scale factor
		*/
        void scaleXYZ(double scfactor);

        /**
            @return The point group
        */
        ConstPointGroupPtr getPointGroup() const;

        PointGroupPtr getPointGroup();

		/** 
            Gets the xyz coordinates
			@param unit The units desired, given as an enum value Bohr, Angstrom, or DefaultUnits
		*/
        ConstRectMatrixPtr getXYZ() const;

        /** 
            @param xyz The coordinate of the dummy atom
            @return The atom number based on zero counting as if it were a real atom.
                    This numbering therefore includes all atoms, dummies and real.
        */
        int register_dummy_atom(ConstVectorPtr xyz);

        /** 
            @return Number of dummy atoms
        */
        int ndummies() const;


		/** 
            Displaces the xyz coordinates by a certain amount. By
            necessity, this recomputes the point group.
			@param disp The matrix of xyz displacements
		*/
        void displace(ConstRectMatrixPtr disp);

        /**
            Displaces the xyz coordinates by a certain amount. By
            necessity, this recomputes the point group.
            @param n The atom number to displace
            @param disp The vector of xyz displacements
        */
        void displaceAtom(int n, ConstVectorPtr disp);

		/** 
            Get the atom object for a given atom number
			@param num The atom number, based on zero counting
			@return The atom object
		*/
        ConstAtomPtr getAtom(int num) const;

        ConstAtomPtr getDummy(int num) const;

        /**
            @param atom 
            @return The atom number based on zero counting
        */
        int getNumber(const ConstAtomPtr& atom) const;

		/** 
            @param units
            @param indent_level How many tabs to add before the geometry
            @return A neatly formatted string of the xyz coordinates. Atom labels not included.
		*/
        std::string getXYZString(gigstr units = angstrom, int indent_level = 0) const;

        /**
            @return Whether the molecule is linear
        */
        bool isLinear() const;

        /**
            @return The number of internal coordinates for the molecule. 3N-5 linear, 3N-6 nonlinear, 0 monoatomic.
        */
        int ninternals() const;

        /**
            @return The number of atoms
        */
        int natoms() const;

		/** 
            @return A copy of the molecule with no pointer links to the old one 
        */
        MoleculePtr copy() const;

        void
        recomputePointGroup(const ConstMoleculePtr& mol);

		/** 
            Gets the mass associated with a given atomic symbol 
			@param symbol The atomic symbol
			@return The corresponding mass
		*/
		static double getMass(gigstr symbol);

        
        /**
            Based on the user-specified units in the input file, returns a unit
            conversion for the specified units.
            @param units The units desired after multiplication by the conversion factor.
                         For example, if input units are "bohr" and "angstrom" is param,
                         function returns 0.529...
        */
        static double getConversion(gigstr units);

        virtual ~Molecule();
};

}

#endif
