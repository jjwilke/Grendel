
#ifndef intdif_fileio_h
#define intdif_fileio_h

#include <stdio.h>
#include <iostream.h>
#include <fstream>
#include <sstream>
#include <src/coordinates.h>
#include <src/molecule.h>

namespace gigide {

/** 
    Given a filename, it sends back a string with all of the file text
	@param filename The full path of the file
    @param comment The comment delimiter for ignoring lines
	@return The text of the file as a string minus comments
*/
std::string getFileText(gigstr filename, gigstr comment);

/**
    Get the the list of atoms and xyz coordinates. Converts coordinates
    based on input parameter.
    @param vector that will contain the atom list upon return
    @param filetext
    @param inputunits
    @param outputunits
    @return XYZ matrix
    @throw GigideException If not properly formatted matrix
    @throw
*/
RectMatrixPtr getXYZMatrix(std::vector<std::string>& atomlist,
                           gigstr filetext,
                           gigstr inputunits,
                           gigstr outputunits);

/**
    Given a block of text, read all numerical values and create a matrix
    @param filetext
    @param nrow
    @param ncol
    @return nrow x ncol matrix
    @throw GigideException If nrow * ncol does not match the total number of values
*/
RectMatrixPtr getMatrix(gigstr filetext, mdim_t nrow, mdim_t ncol);

/**
    Given a block of text, read all numerical values and create vector
    @param text
    @return Vector of values
*/
VectorPtr getVector(gigstr text);

class DerivativeIterator;
class ForceField;

/** Encapsulates an object for building output files */
class FileWriter
{
    protected:
        /** The contents of the file are stored in a std::string stream */
        std::stringstream contents_;

    public:
        /** Constructor */
        FileWriter();

        /** 
            Commits the contents accumulated so far to the given file. See #content_.
            @param filename The filename to print to
        */
        virtual void commit(gigstr filename);

        /** 
            Adds new contents to the file 
            @param newtext The text to add to the file
        */
        void add(gigstr newtext);

        /**
            Add a new line to the file.
        */
        void addEndl();

};

/**
    @class InputFileWriter
*/
class InputFileWriter : public FileWriter
{
    private:
        /**
            Create a force constant file matching ACES fcmfinal / PSI file15
            format.
            @param mol
            @param fc The set of force constants
            @param filename
        */
        static void makeFCFile(
            const ConstMoleculePtr& mol, 
            ConstRectMatrixPtr fc,
            gigstr filename
        );

    protected:
        ConstMoleculePtr mol_;

    public:
        /** 
            Constructor
            @param mol      The molecule for the input file
        */
        InputFileWriter(const ConstMoleculePtr& mol);

        /** 
            Adds a header for the file. By defalt, does nothing.
        */
        virtual void writeHeader();

        /** 
            Adds a footer for the file. By default, adds endl.
        */
        virtual void writeFooter();
        
        /**
            Finish off a section.  By default does nothing.
        */
        virtual void closeSection();

        /**
            Make a file containing the xyz coordinates and gradients.
            @param mol
            @param xyzgradients
        */
        static void makeFile11(const ConstMoleculePtr& mol, ConstRectMatrixPtr xyzgradients);

        /**
            Make a file containing the xyz force constants. This writes values
            consecutively with 3 values per row.
            @param mol
            @param fc
        */
        static void makeFile15(
            const ConstMoleculePtr& mol, 
            ConstRectMatrixPtr fc
        );

        /**
            Make a file containing the internal coordinate values
            and gradients
            @param mol
            @param coords
            @param gradients
        */
        static void makeFile12(
            const ConstMoleculePtr& mol, 
            const Set<ConstInternalCoordinatePtr>& coords,
            ConstVectorPtr gradients
        );

        /**
            Make a file contains the internal coordinate force constants. This writes values
            consecutively with 3 values per row.
            @param mol
            @param fc
        */
        static void makeFile16(
            const ConstMoleculePtr& mol, 
            ConstRectMatrixPtr fc
        );
};

/** 
    Encapsulates an object used for creating intdif input files 
*/
class IntderWriter : public InputFileWriter
{
        public:
            /**
                The type of tranformation for intder to perform 
            */
            typedef enum { TransformCartesians = 0, TransformFileInternals = 1, TransformInternals = 2 } TransformationType;

        private:
            /**
                Whether or not at a stationary point
            */
            bool statpt_;

            bool has_symms_;

            /**
            */
            Set<ConstSimpleInternalCoordinatePtr> simples_;

            /**
                This will be empty if internal coordinates are the same as simples
            */
            Set<ConstInternalCoordinatePtr> coords_;

            /**
                Iterator for derivatives to print
            */
            ConstDerivativeIteratorPtr iter_;

            /**
            */
            TransformationType type_;

            /**
                The level of derivative to transform
            */
            int nder_;

            /**
                Whether to include xyz coordinates in intder input file
            */
            bool includexyz_;
        
        public:
            /** 
                Constructor
                @param mol      The molecule upon which finite difference are being done
                @param statpt   Whether or not the object being treated is a stationary point
                @param type     The transformation type, whether cartesian to internals or internals to cartesians
                @param nder     The level of derivative being transform
                @param iter     The iterator for printing derivatives in the input file.  This may be null 
                                if no internal coordinate derivatives are to be printed
                @param simples
                @param coords  
                @param includexyz Whether to print xyz coordinates in input file.  If false, gets coordinates from file11.
            */
            IntderWriter(const ConstMoleculePtr& mol,
                         bool statpt,
                         TransformationType type,
                         int nder,
                         const ConstDerivativeIteratorPtr& iter,
                         const Set<ConstSimpleInternalCoordinatePtr>& simples,
                         const Set<ConstInternalCoordinatePtr>&  coords,
                         bool includexyz = true);

            /*
            */
            void commit(gigstr filename);

            /** Writes a derivative to the intder input file. */
            void writeDerivatives();

            static void clearFiles();

            /** Writes simple internal coordinates to the input file */
            void writeSimpleCoordinates();

            /** Writes a symmetry internal coordinate description to the input file */
            void writeSymmetryCoordinates();

            /** 
                Writes the xyz coordinates to the file.  No parameter is needed since this will come
                from the molecule passed to the constructor. 
            */
            void writeXYZ();

            /** 
                Writes the xyz coordinates to the file.  No parameter is needed since this will come
                from the molecule passed to the constructor. 
            */
            void writeMasses();

            /**
            */
            void writeHeader();

            /**
                Given a set of gradients and coordinates, perform all steps
                necessary to transform internal gradients to XYZ.
                @param intgrads The internal coordinate gradients
                @param mol
                @param simples
                @param coords
                @return The xyz gradients as an natom x 3 matrix
            */
            static RectMatrixPtr transformGradientsToCartesian(
                ConstVectorPtr intgrads,
                const ConstMoleculePtr& mol,
                const Set<ConstSimpleInternalCoordinatePtr>& simples,
                const Set<ConstInternalCoordinatePtr>& coords,
                bool validate = false
            );

            /**
                Given a set of gradients and coordinates, perform all steps
                necessary to transform xyz gradients to internals.
                @param xyzgrads natom x 3 matrix of gradients
                @param mol
                @param simples
                @param coords
                @return The internal gradients as a vector
            */
            static VectorPtr transformGradientsToInternals(
                ConstRectMatrixPtr xyzgrads,
                const ConstMoleculePtr& mol,
                const Set<ConstSimpleInternalCoordinatePtr>& simples,
                const Set<ConstInternalCoordinatePtr>& coords,
                bool validate = false
            );

            /**
                Given a set of force constants and coordinates, perform all steps
                necessary to transform internal fc to xyz fc.
                @param intfc ncoord x ncoord matrix of force constants
                @param intgrads Internal gradients. These are necessary if not at a stationary point.
                @param mol
                @param simples
                @param coords
                @return The xyz force constants as a matrix
            */
            static RectMatrixPtr transformForceConstantsToCartesian(
                ConstRectMatrixPtr intfc,
                ConstVectorPtr intgrads,
                const ConstMoleculePtr& mol,
                const Set<ConstSimpleInternalCoordinatePtr>& simples,
                const Set<ConstInternalCoordinatePtr>& coords,
                bool validate = false
            );

            /**
                Given a set of force constants and coordinates, perform all steps
                necessary to transform xyz fc to internal fc.
                @param xyzfc 
                @param xyzgrads Gradients are necessary if not at a stationary point.
                @param mol
                @param simples
                @param coords
                @return The internal force constants as a matrix
            */
            static RectMatrixPtr transformForceConstantsToInternals(
                ConstRectMatrixPtr xyzfc,
                ConstRectMatrixPtr xyzgrads,
                const ConstMoleculePtr& mol,
                const Set<ConstSimpleInternalCoordinatePtr>& simples,
                const Set<ConstInternalCoordinatePtr>& coords,
                bool validate = false
            );

            void writeFooter();

        private:
            void closeSection();

};

/**
    @clas AnharmWriter
*/
class AnharmWriter : public InputFileWriter
{
    private:
        /** The force field object that anharm is to use */
        ConstForceFieldPtr qff_;

    public:
        /**
            @param mol
            @param qff
        */
        AnharmWriter(const ConstMoleculePtr& mol,
                     const ConstForceFieldPtr& qff);

        /**
            @param filename
        */
        void commit(gigstr filename);

        /**
        */
        void writeHeader();

        /**
        */
        void writeResonanceCutoffs();

        /**
        */
        void writeXYZ();
};

}

#endif
