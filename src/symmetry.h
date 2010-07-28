
#ifndef gigide_symmetry_h
#define gigide_symmetry_h

#include <vector>
#include <string>
#include <map>
#include <src/defines.h>
#include <src/keyword.h>
#include <src/coordinates.h>
#include <src/gigmatrix.h>
#include <src/archive.h>

namespace gigide {

class Molecule;
class XYZPoint;
/**
    @class SymmetryOperation

    Class encapsulating symmetry operation
*/

class Rotation;
class ImproperRotation;
class Reflection;
class Inversion;
class IdentityElement;

class SymmetryOperation : public smartptr::Serializable {
    
    public:
        /** Enum giving the type of symmetry operation */
        typedef enum {rotation, improper_rotation, reflection, inversion, identity} OperationType;

        /**
            @param mol
            @param type
        */
        SymmetryOperation(
            OperationType type
        );

        SymmetryOperation(
            const ArchivePtr& parser
        );

        void
        serialize(
            const ArchivePtr& ptr
        ) const;

        /**
            @return The operation type
        */
        OperationType type() const;

        /**
            @return The 3x3 orthogonal matrix
        */
        ConstRectMatrixPtr matrix() const;

		/** 
            Gets the permutation matrix associated with a given symmetry operation.
            Suppose we have methane with the numbering C1,H2,H3,H4,H5. If we have
            a C3 operation that maps 2->3, 3->4, 4->2.  If the operation permutes
            (i->j) then \f$P_{ij} = 1\f$ in the permutation matrix.  Otherwise,
            the matrix element is zero.  For the symmetry operation given above,
            the matrix is
            \f{eqnarray*}{
                1 & 1\\
                1 & 0
            \f}
			@param oper The symmetry operation
            @throw If the symmetry operation is not valid for the molecule
		*/
        RectMatrixPtr getPermutationMatrix(const ConstMoleculePtr& mol) const;

        /**
            @return The name of the symmetry operation
        */
        std::string name() const;

        /**
            @return A description of the symmetry operation
        */
        virtual std::string description() const;

        /**
            Perform numerical check on equivalence of two operations.

            @param oper
            @return Whether two operations are equivalent.
        */
        bool isEquivalent(const ConstSymmetryOperationPtr& oper) const;

        /**
            Create a new symmetry operation by conjugation.  In particular, if this operation is \f$M\f$
            and the parameter is \f$T\f$, this computes \f$ T M T^{-1}\f$.
            @param oper The transformation T
            @param oper The conjugated matrix \f$ T M T^{-1}\f$
        */
        SymmetryOperationPtr conjugate(const ConstSymmetryOperationPtr& oper) const;

        /**
            @param os The stream to print to
        */
        void print(std::ostream& os = std::cout) const;

        SymmetryOperationPtr inverse() const;

        SymmetryOperationPtr multiply(const SymmetryOperationPtr& r) const;

        /**
            @param type
            @return A descriptive string for the operation type
        */
        static std::string opname(OperationType type);
        
        /**
            Return 
        */
        static SymmetryOperationPtr identity_op();
        static SymmetryOperationPtr c2x_op();
        static SymmetryOperationPtr c2y_op();
        static SymmetryOperationPtr c2z_op();
        static SymmetryOperationPtr inversion_op();
        static SymmetryOperationPtr sigmaxy_op();
        static SymmetryOperationPtr sigmaxz_op();
        static SymmetryOperationPtr sigmayz_op();

        static void init_statics();

        static SymmetryOperationPtr buildOperation(
            ConstRectMatrixPtr m
        );

    protected:
        /**
            @param vec
            @return A nice-looking string of the axis
        */
        static std::string axis_str(const VectorPtr& vec);

        /**
            see type();
        */
        OperationType type_;

        /**
            see name();
        */
        std::string name_;

        /**
            The 3x3 orthogonal matrix defining the transformation of axes
        */
        RectMatrixPtr matrix_;

        static SymmetryOperationPtr identity_op_;
        static SymmetryOperationPtr c2x_op_;
        static SymmetryOperationPtr c2y_op_;
        static SymmetryOperationPtr c2z_op_;
        static SymmetryOperationPtr inversion_op_;
        static SymmetryOperationPtr sigmaxy_op_;
        static SymmetryOperationPtr sigmaxz_op_;
        static SymmetryOperationPtr sigmayz_op_;
        static bool initdone_;

        static RotationPtr
        buildRotation(
            ConstVectorPtr revals,
            ConstRectMatrixPtr evecs,
            ConstRectMatrixPtr m
        );

        static ReflectionPtr
        buildReflection(
            ConstVectorPtr revals,
            ConstRectMatrixPtr evecs,
            ConstRectMatrixPtr m
        );

        static ImproperRotationPtr
        buildImproperRotation(
            ConstVectorPtr revals,
            ConstRectMatrixPtr evecs,
            ConstRectMatrixPtr m
        );


};

/**
    @class ImproperRotation
*/
class ImproperRotation : public SymmetryOperation {
    
    private:
        /**
            The axis of rotation
        */
        VectorPtr axis_;

        /**
            Defines the type of rotation, e.g if order is 4, this
            represents an \f$S_4\f$ operation.
        */
        int order_;
        
        int exponent_;

        void init();

    public:
        ImproperRotation(
            ConstVectorPtr axis,
            int order,
            int exponent = 1
        );

        ImproperRotation(
            const ArchivePtr& parser
        );

        void serialize(const ArchivePtr& writer) const;

        std::string description() const;

        int order() const;

        int exponent() const;

        ConstVectorPtr axis() const;

};

/**
    @class Rotation
*/
class Rotation : public SymmetryOperation {

    private:
        /**
            The axis of rotation
        */
        VectorPtr axis_;

        /**
            Defines the type of rotation, e.g if order is 4, this
            represents an \f$S_4\f$ operation.
        */
        int order_;

        int exponent_;

        void init();

    public:
        /**
            @param mol
            @param axis The axis of rotation
            @param order See #order_
            @param exponent Defines the extent of rotation. If exponent is 3 and order is 4, e.g.
                            then this defines \f$S_4^3\f$.
        */
        Rotation(
            ConstVectorPtr axis,
            int order,
            int exponent = 1
        );

        Rotation(
            const ArchivePtr& parser
        );

        void serialize(const ArchivePtr& writer) const;

        /**
            @return See #order_;
        */
        int order() const;

        int exponent() const;

        const VectorPtr& axis() const;

        std::string description() const;


};

/**
    @class Reflection
*/
class Reflection : public SymmetryOperation {

    private:
        /** The axis normal to the plane */
        VectorPtr axis_;

        void init(ConstVectorPtr axis);

    public:
        /**
            @param mol
            @param axis The normal vector to the plane
        */
        Reflection(
            ConstVectorPtr axis
        );

        Reflection(
            const ArchivePtr& parser
        );

        void serialize(const ArchivePtr& writer) const;

        /**
            Constructs an axis based on 3 points in a plane.

            @param mol
            @param pt1 The first point in the plane
            @param pt2 The first point in the plane
            @param pt3 The first point in the plane
        */
        Reflection(
            const XYZPointPtr& pt1,
            const XYZPointPtr& pt2,
            const XYZPointPtr& pt3
        );

        std::string description() const;

        const VectorPtr& axis() const;

};

/**
    @class Inversion
*/
class Inversion : public SymmetryOperation {
        
    private:
        void init();

    public:
        /**
            @param mol
        */
        Inversion();

        Inversion(
            const ArchivePtr& parser
        );

        void serialize(const ArchivePtr& writer) const;

};


/**
    @class IdentityElement
*/
class IdentityElement : public SymmetryOperation {
    
    private:
        void init();

    public:
        /**
            @param mol
        */
        IdentityElement();

        IdentityElement(
            const ArchivePtr& parser
        );

        void serialize(const ArchivePtr& writer) const;
};

/**
    @class PointGroupClass
*/
class PointGroupClass : public smartptr::Serializable {

    private:
        /**
            The type of symmetry operation contained in the class, e.g. rotation
        */
        SymmetryOperation::OperationType type_;

        /**
            The set of symmetry operations in the class
        */
        Set<ConstSymmetryOperationPtr> members_;

        /**
        */
        ConstMoleculePtr mol_;

    public:
        /**
            @param mol
            @param type
        */
        PointGroupClass(
            const ConstMoleculePtr& mol,
            SymmetryOperation::OperationType type
        );

        
        PointGroupClass(
            const ArchivePtr& parser
        );

        void serialize(const ArchivePtr& writer) const;

        /**
            @return The number of symmetry operations in the class
        */
        int order() const;

        /**
            @return The operation type
        */
        SymmetryOperation::OperationType type() const;

        /**
            No sanity checks are currently done

            @param oper
        */
        void addOperation(const ConstSymmetryOperationPtr& oper);

        /**
            @param oper
            @return Whether the operation already exists
        */
        bool hasOperation(const ConstSymmetryOperationPtr& oper) const;


        /**
            @return A representation symmery operation for the class. Null ptr if empty.
        */
        ConstSymmetryOperationPtr getRepresentative() const;

        /**
            @param symmops The vector to put the current classes symmetry operations into
        */
        void getSymmetryElements(
            std::vector<ConstSymmetryOperationPtr>& symmops
        ) const;

#if 0
        /**
            Conjuagte all element of the current class by the given operator
            and add them to the class if new.

            @param oper
        */
        void addConjugate(SymmetryOperationPtr oper);
#endif

        /**
            @param os
        */
        void print(std::ostream& os = std::cout) const;

        bool
        testMembership(
            const ConstSymmetryOperationPtr& op,
            const Set<ConstSymmetryOperationPtr>& symmops
        ) const;
};

/**
    @class PointGroup
*/
class PointGroup : public smartptr::Serializable {

    private:
        /**
            The set of all classes in the point groups
        */
        std::vector<PointGroupClassPtr> classes_;

        /**
        */
        std::vector<ConstSymmetryOperationPtr> symmops_;

        ConstMoleculePtr mol_;

        /**
            A string identifying the point group
        */
        std::string pgstr_;

        /**
            Whether the point group has been "closed", i.e. now new element
            can be created as products of current elements.
        */
        bool valid_;

        /**
            Whether or not classes have been computed
        */
        bool classescomputed_;

        bool isabelian_;

        /**
        */
        static bool initdone_;
        /**
        */
        static std::map<std::string, int> symmaxes_;
        /**
        */
        static std::map<std::string, int> symmplanes_;

        /**
            The set of coordinate subspaces grouped by degeneracy
        */
        std::vector<CoordinateSubspacePtr> subspaces_;

        /**
        */
        std::vector<PointGroupPtr> subgroups_;

    private:
        static void
        init_statics();

        void
        formClasses();

        /**
            Add a new class containing the symmetry operation
        */
        void addClass(const ConstSymmetryOperationPtr& op);

        bool
        findClass(const ConstSymmetryOperationPtr& op);

        bool computeAbelianness();

    public:
        /**
            @param mol
            @param title Descriptive name (e.g. c2v) for point group
        */
        PointGroup(
            const ConstMoleculePtr& mol, 
            std::string title = ""
        );

        PointGroup(
            const ArchivePtr& parser
        );

        virtual ~PointGroup();

        void serialize(const ArchivePtr& writer) const;

        /**
            @return The name of the point group
        */
        std::string name() const {return pgstr_;}

        /**
            @return The number of elements in the group
        */
        int order() const;

        /**
            @return The number of classes
        */
        int nclasses() const;

        /**
            Adds a new operation to the group.  This creates a new point group class.
            Future versions should check if the class already exists.

            @param op
        */
        void addOperation(const ConstSymmetryOperationPtr& op);

        /**
            @param os
        */
        void printClasses(std::ostream& os = cout) const;

        /**
            @param os
        */
        void print(std::ostream& os = cout) const;

        /**
            @return Whether this corresponds to an abelian point group
        */
        bool isAbelian() const;

        /**
            Create all possible combinations of symmetry elements
            and verify that the groups is closed under multiplication.
        */
        void formClosedGroup();

        void
        formAbelianClasses();

        void close();

        /** 
            Check to see whether the operation exists in the point group.  The check
            is done by numerically checking the 3x3 transformation matrices to see if they
            are equivalent within a numerical cutoff.
            @param op The symmetry operation to check for
            @return If operation is part of point group
            @throw Programming error if called before formClosedGroup is called to finish
                   initialization.
        */
        bool hasOperation(const ConstSymmetryOperationPtr& op) const;

        bool hasOperation(ConstRectMatrixPtr op) const;

        /**
            @param classvec The vector to hold classes
            @throw Programming error if called before formClosedGroup is called to finish
                   initialization.
        */
        void getClasses(std::vector<PointGroupClassPtr>& classvec) const;

        /**
            @param symmops The vector to hold the symm ops
            @throw Programming error if called before formClosedGroup is called to finish
                   initialization.
        */
        void getSymmetryElements(
            std::vector<ConstSymmetryOperationPtr >& symmops
        ) const;

        /**
            Computes whether the set of characters contains a contribution
            from the totally symmetric irrep.  This assumes the vector and the list
            of point group classes are aligned.
            @param chars
            @return If vec has totally symmetry irrep contribution
        */
        bool isTotallySymmetric(const VectorPtr& chars) const;

        void recompute();

        void
        getSubspaces(
            std::vector<CoordinateSubspacePtr>& subspaces
        ) const;

        void
        formBasis(
            const Set<ConstSymmetryInternalCoordinatePtr>& internals,
            std::vector<CoordinateSubspacePtr>& subspaces
        ) const;

        void
        formBasis(
            const Set<ConstSimpleInternalCoordinatePtr>& simples
        );

        Set<PointGroupPtr>
        getSubgroups() const;

        void
        addSubgroup(
            PointGroupPtr pg
        );


};

class Irrep {

    private:
        PointGroupPtr pg_;

        VectorPtr chars_;
    
    public:
        Irrep();
};

class CharacterTable {
};

} //end namespace

#endif
