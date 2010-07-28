
#ifndef gigide_qff_h
#define gigide_qff_h

#include <map>
#include <vector>
#include <string>
#include <src/molecule.h>
#include <src/coordinates.h>
#include <src/pyxml.h>
#include <src/keyword.h>
#include <src/taylor.h>
#include <src/fileio.h>
#include <src/fit.h>
#include <src/gigmatrix.h>
#include <src/deftypes.h>
#include <src/archive.h>

namespace gigide {

class Derivative;
class Displacement;
class DerivativeIterator;
class DisplacementIterator;

/** Encapsulates a force field computation */
class ForceField : public smartptr::Serializable
{
    private:
        /** The displacement sizes for each internal coordinate.  The units are set by the input file */
        std::vector<double> disp_sizes_;

        /** An iterator for generating and iterating derivatives of values */
        DerivativeIteratorPtr compute_iter_;

        /** A derivative that inclues the necessary derivaties for generating a certain
            level of robustness */
        DerivativeIteratorPtr fit_iter_;

        /** A derivative iterator for the provided values. This may do nothing
            if energies are provided, or it may do a lot if 2nd derivatives are provided. */
        DerivativeIteratorPtr value_iter_;

        /** The derivative iterator for all derivatives */
        DerivativeIteratorPtr deriv_iter_;

        /** An iterator for generating and iterating displacements */
        DisplacementIteratorPtr disp_iter_;

        /** The set of all internal coordinates */
        Set<ConstInternalCoordinatePtr> coords_;

        /** The set of all simple internals comprising the internal coordinates */
        Set<ConstSimpleInternalCoordinatePtr> simples_;

        /** The molecule for which a force field is being computed.  This is the undisplaced molecule */
        ConstMoleculePtr mol_;

        /** The force constant matrix... only applicable to when 2nd derivatives or higher are computed */
        SymmMatrixPtr fc_;

        /** The G matrix, describing the combination of two B matrix elements */
        SymmMatrixPtr G_;

        /** The gradients */
        VectorPtr grads_;

        /** The energy of the center molecule */
        double energy_;

        /** Depending on the type of derivatives being computed, we may or may not be able to throw out displacments */
        bool include_zeros_;

        int nderiv_;

        int nvalue_;

        bool getXMLGradients(const XMLParser& xml, RectMatrixPtr& gradients);
        bool getXMLEnergy(const XMLParser& xml, double& energy);
        bool getXMLForceConstants(const XMLParser& xml, RectMatrixPtr& fc);
        bool getXMLXYZ(const XMLParser& xml, RectMatrixPtr& xyz);

    private:
        void testDerivative(VectorPtr fit, DisplacementIteratorPtr dispit, int nderiv);
        std::vector<double> getDisps(std::vector<int> disp_numbers);

    public:
        /** Constructor.  The constructor builds a displacement and derivative iterator.
            @param center_mol The molecule at the center of the displacements
            @param coords     The set of internal coordinates to displace
            @param simples    The set of simple internals comprising the coordinates
        */
        ForceField(const ConstMoleculePtr& center_mol, 
                   const Set<ConstInternalCoordinatePtr>& coords,
                   const Set<ConstSimpleInternalCoordinatePtr>& simples);

        ForceField(const ArchivePtr& parser);

        /**
        */
        void
        serialize(const ArchivePtr& writer) const;

        /** Builds the B matrix from the internal coordinates */
        RectMatrixPtr buildB() const;

        static SymmMatrixPtr formGMatrix(
            const Set<ConstInternalCoordinatePtr>& coords
        );


        static void computeNormalModes(
            ConstSymmMatrixPtr F,
            ConstSymmMatrixPtr G,
            VectorPtr& freqs,
            RectMatrixPtr& evecs
        );

        /** Generates the displacements */
        void generateDisplacements();

        void writeDisplacementsToFile(std::string filename) const;

        void writeMoleculeToDispcart(
            std::ofstream& dispcart,
            const ConstMoleculePtr& mol,
            int dispnum
        ) const;

        void writeIntderFile(const ConstMoleculePtr& mol,
                             std::string filename,
                             IntderWriter::TransformationType type,
                             bool statpt) const;

        void assignEnergy(double energy);

        void assignGradients(ConstRectMatrixPtr xyzgrads);

        void assignForceConstants(ConstRectMatrixPtr xyzfc, ConstRectMatrixPtr xyzgrads);

        void writeAnharmFile() const;

        /** Builds the fitting matrix for computing derivatives from the finite differences */
        void buildFit();

        void buildNumericalFit();

        void buildFormulaFit();
    
        /** Builds the first and second derivative arrays */
        void buildArrays();

        /** Builds the derivative values by differentiating a given derivative value
            @param The specific derivative to undergo finite differences
        */
        void computeDerivatives(DerivativePtr value);

        void computeFrequencies();

        /** Returns the number of internal coordinates 
            @return The number of internal coordinates
        */
        int ncoords() const;

        int robustness(int derivlevel) const;
        
        int derivlevel_compute() const;

        /** Computes all the derivatives.  This has the precondition that all necessary values are in
            properly formatted output files
        */
        void compute();

        /** Read data for displacements from an xml file
            @param filename The name of the xml file to parse
        */
        void readXMLData(std::string filename = "data.xml");

        /** From an XML node representing a displacement, adds all the necessary data 
            @param node The XML node representing a displacement
            @param disp The displacement to store data in 
        */
        void readXMLDisplacement(const XMLParser& node, bool& extrazero);

        int nderiv() const;

        void validateFit();

        void print(std::ostream& os = std::cout) const;

        void validateRow(
            std::vector<TaylorTermPtr>& terms,
            ConstVectorPtr coefs
        );

        /** Gets all the values of a given derivative from all the computed displacements */
        ConstVectorPtr getDerivativeValues(const ConstDerivativePtr& value_deriv) const;

        ConstVectorPtr gradients() const;

        double energy() const;

        ConstSymmMatrixPtr fc() const;
};

} //end namespace sc

#endif
