
#ifndef gigide_derivative_h
#define gigide_derivative_h

#include <map>
#include <vector>
#include <src/coordinates.h>
#include <src/deftypes.h>
#include <src/iter.h>
#include <src/fit.h>

namespace gigide {

/** 
    @class Derivative

    Encapsulates a derivative 
*/
class TaylorTerm;
class DerivativeIterator;

class Derivative : public smartptr::Serializable
{
    public:
        /**
            enum describing unique derivative types.  See #typelabel_.
        */
        typedef enum {Dertype_i, Dertype_ii, Dertype_iii, Dertype_iiii,
                      Dertype_ij, Dertype_iij, Dertype_iiij,Dertype_iijj,
                      Dertype_ijk, Dertype_iijk,
                      Dertype_ijkl} DerivativeType;

    private:
        /** 
            The individual derivative levels for the coordinates in this derivative.
			For example, if you have the derivative 
            \f[
                \frac{d^4 E}{dx_1 dx_2 dx_4^2}
            \f]
			<1 1 0 2> to indicate the levels
		*/
        std::vector<int> levels_;

        /**
            See #levels_. Same vector with zeroes removed.
        */
        std::map<int,int> nonzero_levels_;

        /**
            Stores a string map to describe the derivative.  This is useful in recognizing
            derivatives as being equivalent. As an example,
            the derivative
            \f[
                \frac{d^4 E}{dx_1 dx_2 dx_4^2}
            \f]
            would be the map <br>
            map["i"] = 1 <br>
            map["j"] = 1 <br>
            map["k"] = 2 <br>
            to identify the derivative as having the structure
            \f[
                \frac{d^4 E}{dx_i dx_j dx_k^2}
            \f]
            The following derivatve
            \f[
                \frac{d^4 E}{dx_1 dx_2 dx_3 dx_4}
            \f]
            would have the map <br>
            map["i"] = 1 <br>
            map["j"] = 1 <br>
            map["k"] = 1 <br>
            map["l"] = 1 <br>
        */
        std::map<std::string, int> indices_;

        /**
            A string label representing the derivative. See #indices_. For example, the derivative
            \f[
                \frac{d^4 E}{dx_1 dx_2 dx_4^2}
            \f]
            would have typelabel "ijkk"
        */
        std::string typelabel_;

        std::vector<DerivativePtr> equiv_derivs_;

        /** 
            The taylor coefficient that would be attached to this derivative. For example,
            the derivative
            \f[
                \frac{d^4 E}{dx_1^3 dx_2}
            \f]
            would have the coefficient
            \f[
                \frac{1}{3! 1!}
            \f]
		*/
        double coeff_;

        /** 
            The level of the derivative, i.e. 1st,2nd,3rd,etc... 
        */
        int level_;

        /**
            The molecule the derivative is based on.
        */
        ConstMoleculePtr mol_;

		/** 
            The set of coordinates the derivatives are based on 
        */
        smartptr::Set<ConstInternalCoordinatePtr> coords_;

		/**
            Depending on symmetry to be determined later, the derivative might be rigorously zero 
        */
        bool nonzero_;

        /**
            The fitting function for computing the derivative
        */
        FitPtr fit_;

        /**
            The value of the derivative.  This is assumed zero until assigned later by #fit_.
        */
        double deriv_value_;

        /**
            Whether or not the value has been assigned
        */
        bool assigned_;

        /**
            List of letters defining abstract derivative types. The vector is therefore
            <"i","j","k",l",...>.  See #typelabel_.
        */
        static std::vector<std::string> letters_;

        /**
            Standard init check
        */
        static bool initdone_;

        /**
            The derivative type.  See #DerivativeType and #typelabel_;
        */
        DerivativeType dertype_;

        /**
            Whether the derivative is unique, or if it is equivalent to another derivative
            and therefore does not need to be computed.
        */
        bool unique_;

        /**
            Maps #typelabel_ to the appropriate enum #DerivativeType
        */
        static std::map<std::string, DerivativeType> typemap_;

    private:
        /**
            Init static variables
        */
        static void init_statics();

    public:
		/** Constructor
			@param levels The individual derivative levels of a set of coordinates. See #levels_.
            @param coords   The set of internal coordinates the derivative is built on
            @mol The molecule
		*/
        Derivative(
            const std::vector<int>& levels,
            const smartptr::Set<ConstInternalCoordinatePtr>& coords,
            const ConstMoleculePtr& mol
        );

        Derivative(const XMLArchivePtr& parser);

        void serialize(const XMLArchivePtr& writer) const;

        /** 
            Computes the coefficient of the derivative in the Taylor expansion for a given set of displacements.
            For example, for the set of displacements <0.75,1.0,0,0.5> for the derivative
            \f[
                \frac{d^4 E}{dx_1 dx_2 dx_4^2}
            \f]
            we would compute
            \f[
                C = \frac{1}{1! 1! 0! 3!} \times 0.75 \times 1.0 \times (0.5)^2
            \f]
            @param disps The vector of displacements
            @return The Taylor coefficient
        */
        double taylorCoeff(std::vector<double>& disps) const;

        /** 
            Computes the distplacement coefficient of the derivative in the Taylor expansion for a given set of displacements.
            For example, for the set of displacements <0.75,1.0,0,0.5> for the derivative
            \f[
                \frac{d^4 E}{dx_1 dx_2 dx_4^2}
            \f]
            we would compute
            \f[
                C = 0.75 \times 1.0 \times (0.5)^2
            \f]
            @param disps The vector of displacements
            @return The Taylor coefficient
        */
        double displacementCoeff(std::vector<double>& disps) const;

        /**
            Builds a string representing the polynomial term associated
            with the derivative. For the example derivative
            \f[
                \frac{d^4 E}{dx_1 dx_2 dx_4^2}
            \f]
            the string would be "x1x2(x4)^2"
        */
        std::string polynomialString() const;

        bool matches(const std::vector<int>& levels) const;

        /** 
            Returns a descriptive label for the derivatiave. For the example
            \f[
                \frac{d^4 E}{dx_1 dx_2 dx_4^2}
            \f]
            The label would be "1 1 0 2".  See #levels_.
            @return The label
        */
        std::string label() const;

        /**
            Prints details about the levels and value of the derivative. 
            See #levels_ and #value_.
            @param os The ouput stream to print to
        */
        void printDetail(std::ostream& os = std::cout) const;


        void mapEquivalentDerivatives(ConstDerivativeIteratorPtr iter);

        /**
            Sends back the coordinate number of a given letter index.  See #typelabel_.
            For example, the derivative 
            \f[
                \frac{d^4 E}{dx_1 dx_2 dx_4^2}
            \f]
            would have type label "ijkk" and index("k") would return 4.
            @param The string in #typelabel_
            @return The coordinate number associated with the letter
        */
        int index(const std::string& str) const;

        void
        addEquivalentDerivative(DerivativePtr equiv);

        /** 
            Prints a description of the derivative 
            @param os The output stream to print to
        */
        void print(std::ostream& os = std::cout) const;

        /**     
            Returns whether the derivative is rigorously nonzero by symmetry
            @return If the derivative is nonzero by symmetry
        */
        bool nonzero() const;

        /**     
            Returns whether the derivative is rigorously zero by symmetry
            @return If the derivative is zero by symmetry
        */
        bool zero() const;

        /**
            Assign a fitting function to the derivative.
            @param fit The fitting function
        */
        void assignFit(FitPtr fit);

        /**
            Given a level of robustness, validate the fitting function
            to make sure the formula is correct. See Fit::assignFit.
        */
        void validateFit(int robustness);

        /** 
            Returns the levels of the derivative.  See #levels_.
            @return The derivative levels
        */
        std::vector<int> levels() const;

        /** 
            Returns the level of the derivative, i.e. 1st, 2nd, 3rd
            @return Degree of derivative
        */
        int level() const;

        /** 
            Returns a different representation of the derivative. For the derivative
            \f[
                \frac{d^4 E}{dx_1 dx_2 dx_4^2}
            \f]
            The levels representation is <1 1 0 2>. The indices representation would be <1 2 4 4>
            @return The index representation of the derivative
        */
        std::vector<int> indices() const;

        /**
            See #typelabel_
            @return The typelabel of the derivative
        */
        std::string typelabel() const;

        /**
        */
        double compute(ConstDerivativePtr value) const;

        /** 
            Gets the value of the derivative.
            @return Derivative value
            @throw GigideException If the value has not been assigned
        */
        double value() const;

        /** 
            smartptr::Sets the value of the deriavtive
            @param val The value of the derivative
        */
        void setValue(double val);

        /**
            Get the the dertype. See #DerivativeType
            @return The DerivativeType num 
        */
        DerivativeType dertype() const;

        /**
            Print the set of equivalent derivatives
            @param os The output stream to print to
        */
        void printEquivalentDerivatives(std::ostream& os = std::cout) const;

        /**
            smartptr::Set the uniqueness of the derivative. Derivative is
            assumed unique unless otherwise set.
            @param uniq
        */
        void setUnique(bool uniq);

        /**
            Returns whether or not the derivative is unique
            @return The uniqueness of the derivative
        */
        bool isUnique() const;

        /**
            Build a label for the derivative. See label()
            @param levels See #levels_
        */
        static std::string label(const std::vector<int>& levels);

};


/** 
    @class DerivativeIterator

    Iteration 
*/
class DerivativeIterator : public smartptr::Serializable
{
    public:
        typedef Iterator<DerivativePtr, ConstDerivativePtr> const_iterator;

        typedef std::vector<DerivativePtr>::iterator iterator;


    private:
        /** The highest level of derivative */
        int deriv_level_;

        /** When iterating, zero derivatives may or may not be included and this may vary depending on the level
            For example, you may wish to include zero 2nd derivatives, but only non-zero 3rd derivatives. In the iteration,
            all derivatives that are zero and of degree higher than nonzero level will not be included */
        int nonzero_level_;

        /** The complete set of derivatives */
        std::vector<DerivativePtr> derivs_;

        /** The set of internal coordinates */
        smartptr::Set<ConstInternalCoordinatePtr> coords_;

        /** The molecule the derivatives are based on */
        ConstMoleculePtr mol_;

    private:
        /** Generates all possible unique derivative combinations. These are generated in the level representation.
            See Derivative::levels_ for the vector format.
            @param level The level of derivative to generate (i.e. 1st, 2nd, etc...)
            @param combos The vector to place the final results in
        */
        void generateDerivativeCombos(int level, std::vector<std::vector<int> >& combos) const;

        /** Derivatives are generated by recursion.
            @param level The level of derivative combinations being generated, i.e. (1,2,3...)
            @param quanta The number of derivative quanta that have so far been used
                          in the recursion. For example, the vector <1,1,0,0> has quanta 2,
                          <1,0,2,0> has quanta 3 and so forth.
            @param start_index The coordinate number to start distributing quanta on
            @param all_combos The set of all derivative combinations
            @param current_combo The current derivative being built in the recursion. This should copy
                                 so that the reference is reset.  This vector needs to be independent.
        */
        void recurseDerivativeCombos(int level, int quanta, int start_index,
                                     std::vector< std::vector<int> >& all_combos,
                                     std::vector<int> current_combo) const;

    public:
        /** 
            Constructor
            @param level The level of derivative to be computed (1st, 2nd, etc...)
            @param coords The set of coordinates to generate the derivatives for
            @param mol
        */
        DerivativeIterator(
            int level,
            const smartptr::Set<ConstInternalCoordinatePtr>& coords,
            const ConstMoleculePtr& mol
        );

        DerivativeIterator(const XMLArchivePtr& parser);

        void serialize(const XMLArchivePtr& writer) const;

        ~DerivativeIterator();

        /** 
        */
        int nderivs() const;

        /** 
        */
        const_iterator begin() const;

        iterator begin();

        const_iterator end() const;

        iterator end();

        /** 
            Gets the highest level of derivative included in the iteration
            @return The highest level of derivative
        */
        int level() const;

        /**
            Print the details of all derivatives in the iteration. See #nderivs.
            @param os The output stream to print to
        */
        void printDetail(std::ostream& os = std::cout) const;

        /**
            Fetch a derivative matching the levels.  See Derivative::level_ for format.
            @levels The set of levels defining the derivative
            @return The derivative matching the levels
        */
        DerivativePtr getDerivative(const std::vector<int>& levels) const;

        /**
            Fetch a derivative matching the indices.  See Derivative::indices_ for format.
            @levels The set of indices defining the derivative
            @return The derivative matching the indices
        */
        DerivativePtr getDerivativeFromIndices(const std::vector<int>& indices) const;

        /** 
            Find a combination derivative.  We have three different derivative iterators.  The main iterator
            iterates through all the final values that will be computed.  The compute iterator is the iterator
            for computed derivatives.  The value iterator iterates through the provided values.  In other words, suppose
            you wish to compute 3rd derivatives by displaced gradients.  The main iterator iterates all 3rd derivatives,
            the compute iterator iterates all 2nd derivatives, and the value iterator iterates all first derivatives. Given
            a compute derivative and a value derivative, this method returns the combined derivative in the main iterator.
            For example, the computed derivative might be <1,0,1> and the value derivative might be <1,0,0> so the
            combination derivative would be <2,0,1>
            @param comp_deriv The computed derivative
            @param value_deriv The value that was differentiated
            @return The mixed derivative
        */
        DerivativePtr findCombination(const ConstDerivativePtr& comp_deriv, const ConstDerivativePtr& value_deriv) const;

        int ncoords() const;


};

}

#endif
