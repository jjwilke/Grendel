
#ifndef gigide_input_h
#define gigide_input_h
#include <src/qff.h>
#include <vector>
#include <src/coordinates.h>
#include <src/keyword.h>
#include <src/defines.h>
#include <src/geometry.h>
#include <src/gigmatrix.h>

namespace gigide {


/**
    Encapsulates an input file for a particular program.

    @class InputFile 
*/
class InputFile : public smartptr::Countable {
    
    private:
        /**
        */
        std::string filetext_;

    public:
        /**
            @param text The input file text
        */
        InputFile(gigstr text);

        /**
            Gets a keyword iterator for a particular section.  This current assumes that sections appear as <br>
            $section_name <br>
                keyname = keyval <br>
                ... <br>
            $end

            @param section_name Currently case-sensitive
            @return A keyword iterator for all of the input. Iterator is empty if
                    section does not exist.
        */
        KeywordIteratorPtr getKeywords(gigstr section_name) const;

        /**
            Tests whether a section exists.
            @param section_name Current case-sensitive
            @return Whether the section exists
        */
        bool hasSection(gigstr section_name) const;

        /**
            Returns the text of a section. Does not include 
            opening and closing $tags. See getKeywords.
            @param section_name Case-sensitive
            @return The text of the section.  Empty string if section does not exist.
        */
        std::string getSection(gigstr section_name) const;

        /**
            Retrieve the list of all included sections.
            @param names The vector that will hold the section names
        */
        void getSectionNames(std::vector<std::string>& names) const;

        /**
            Get the geometry from a particular section. If the keyword "angstrom"
            or "bohr" is found in the section, the geometry is assumed to have those units.
            The geometry is return with the units specified in the "bond units" keyword.

            @param atomlist Will hold the list of atoms found
            @param section_name
            @return The XYZ coordinates
        */
        RectMatrixPtr getGeometry(std::vector<std::string>& atomlist, gigstr section_name = "geometry") const;

};

/**
    Encapsulates input file operations specific to Gigide

    @class GigideInputFile
*/
class GigideInputFile : public InputFile {

    private:
        std::map<std::string, SimpleInternalCoordinatePtr> label_map_;

        /**
            Loop through the definition of symmetry internal coordinates
            and add them to the set.

            @param name The name of the coordinate to add
            @param keyval The keyword value associated with the coordinate name
            @param mol
            @param coords The coordinate set to append the new coordinate to
            @param simples The set of simple internals to base the symm coords on
            @param label_map The mapping between simple coord names and simple coord objects
        */
        void
        appendCoordinates(
            gigstr name,
            const KeywordValuePtr& keyval,
            const ConstMoleculePtr& mol,
            std::vector<InternalCoordinatePtr>& coords,
            std::map<std::string, SimpleInternalCoordinatePtr>& label_map
        );

        /**
            Read a list of symmetry internal coordinates

            @param keyset The keyword iterator.  Each symm coord is defined by a single keyword value.
            @param label_map The mapping between simple coord names and simple coord objects
            @param simples The set of simple internals to base the symm coords on
            @param symm_coords The coordinate set to append the new coordinate to
            @param mol
        */
        void
        readSymmetryCoordinates(
            KeywordIteratorPtr keyset,
            std::map<std::string, SimpleInternalCoordinatePtr>& label_map,
            std::vector<InternalCoordinatePtr>& symm_coords,
            const ConstMoleculePtr& mol
        );

        /**
            Read a list of symmetry internal coordinates

            @param mol
            @param simples The set of simple internals to base the symm coords on
            @param label_map The mapping between simple coord names and simple coord objects
            @param coords The coordinate set to append the new coordinates to
        */
        void
        addSymmetryInternalCoordinates(
            const ConstMoleculePtr& mol,
            std::vector<std::string>& simples,
            std::map<std::string, SimpleInternalCoordinatePtr>& label_map,
            std::vector<InternalCoordinatePtr>& coords
        );

        /**
            In defining geometry and symmetry features, xyz points in space can
            be defined.  This reads any xyz points defined in $points.

            @param mol
            @param points The map that will hold point names and point values
        */
        void
        readXYZPoints(
            const ConstMoleculePtr& mol,
            std::map<std::string, XYZPointPtr>& points
        );

        /**
            In defining geometry and symmetry features, axes in space can
            be defined.  This reads any axes points defined in $axes.

            @param mol
            @param points The map of previously defined xyz points
            @param axes The map that will hold the defined axes
        */
        void
        readAxes(
            const ConstMoleculePtr& mol,
            std::map<std::string, XYZPointPtr>& points,
            std::map<std::string, AxisPtr>& axes
        );

        /**
            Read the symmetry operations defined in a particular section.

            @param mol
            @param points Previously defined xyz pionts
            @param axes Previoulsy defined axes
            @param symmops The vector holding new symmetry operations
            @param section The text of the section defining the symm op
        */
        void
        readSymmetryOperations(
            const ConstMoleculePtr& mol,
            std::map<std::string, XYZPointPtr>& points,
            std::map<std::string, AxisPtr>& axes,
            std::map<std::string, SymmetryOperationPtr>& symmops
        );

        void
        readSubgroups(
            const ConstMoleculePtr& mol,
            std::map<std::string, SymmetryOperationPtr>& symmops
        );

        void
        geometryError(
            gigstr value,
            gigstr spec
        );
    
    public:
        /**
            @param text
        */
        GigideInputFile(gigstr text);

        /**
            Get the keyword set defined by the $options section.

            @return keyword set for $options ... $end section
        */
        KeywordSetPtr getOptions() const;

        /**
            Validate the input file by verifying syntax and making sure
            only allowed section names are given.

            @throw If input file error is found
        */
        void validateInputFile();

        /**
            Read symmetry operations and notify molecule object of symmetry.

            @param mol
        */
        void
        readSymmetryOperations(
            const ConstMoleculePtr& mol
        );

        /**
            @param mol
            @param simples The vector that will hold simple internal coordinates found
        */
        void
        getSimpleCoordinates(
            MoleculePtr mol,
            std::vector< SimpleInternalCoordinatePtr >& simples
        );

        /**
            @param mol
            @param simples 
            @param final_coords The vector that will hold internal coordinates.  This may end
                                up being the same as simples.
        */
        void
        getSymmetryCoordinates(
            const ConstMoleculePtr& mol,
            std::vector<SimpleInternalCoordinatePtr>& simples,
            std::vector<InternalCoordinatePtr>& final_coords
        );

        /**
            @param atomlist
            @return geometry
        */
        RectMatrixPtr getGeometry(std::vector<std::string>& atomlist);
};

} //end namespace
#endif
