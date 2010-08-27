#ifndef gigide_taylor_h
#define gigide_taylor_h

#include <vector>
#include <iostream>
#include <src/deftypes.h>

namespace gigide {

class Displacement;
class DisplacementIterator;
class Derivative;
class DerivativeIterator;

/**
    @class TaylorTerm

    Each TaylorTerm encapsulates a polynomial in the Taylor expansion of the energy of a given displacement.
    Each TaylorTerm contains a polynomial, prefactor, and derivative component. See #prefactor_ and #indices_.
*/
class TaylorTerm : public smartptr::Countable {

    private:
        /**
            The coefficient of the given polynomial in an expansion.
        */
        double coef_;

        /**
            The indices defining this term in the Taylor series.  For example, the term
            \f[
                \frac{1}{1! 1! 0! 2!} \frac{d^4E}{dx_1 dx_2 dx_4^2} (x_1 - x_1^0) (x_2 - x_2^0) (x_4 - x_4^0)^2
            \f]
            would be represented by <1,2,4,4>. This is the index representation. See Derivative::indices().
        */
        std::vector<int> indices_;

        /**
            The Taylor prefactor associated with the taylor Term. See #indices_.
        */
        double prefactor_;

        DerivativePtr deriv_;

    
    public:
        /**
            Constructor.
            Each displacement energy can be expanded as a Taylor series in 
            a set of Taylor series terms \f$\{T_i\}\f$
            \f[
                E = \sum_i T_i = \sum_i c_i D_i P_i
            \f]
            where \f$c_i\f$ is the coffiecient of the ith term, \f$D_i\f$ is the ith
            derivative, and \f$P_i\f$ is ith polynomial.  See #indices_ for an
            example.  This class represent a single term \f$T_i\f$ in the expansion.
            @param indices The indices defining the term. See Derivative::indices()
        */
        TaylorTerm(
            const std::vector<int>& indices
        );

        /**
        */
        void
        assignDerivative(
            DerivativeIteratorPtr deriter
        );

        /**
            Accumulate the contributions to this Taylor term that would come from
            a displacement with a given coefficient in a fitting function.
            @param disp The displacement contributing
            @param coef The coefficient of the displacement in the fitting function
        */
        void
        accumulate(
            ConstDisplacementPtr disp,
            double coef
        );

        /**
            Return the level of the underlying derivative in the term. See #indices_.
            @return The derivative level
        */
        int level() const;

        /**
            Compute the coefficient of this Taylor term in the energy series
            expansion of a given displacement.  The coefficient includes
            the contribution from the standard #prefactor_.
            @param disp The displacement to compute the coef of
            @return The coefficient in the expansion
        */
        double
        coef(ConstDisplacementPtr disp) const;

        /**
            Compute the value of the taylor series term.  This is coef()
            times the value of #deriv_
        */
        double
        value(ConstDisplacementPtr disp) const;

        /**
            Returns a string representation of the polynomial term associated with a displacement.
            @param disp The displacement to get the polynomial for
            @param string The string representation of the polynomail
        */
        std::string
        polynomialString(ConstDisplacementPtr disp) const;

        /**
            Return the coefficient of the Taylor term
            after accumulating contributions from several displacements.
            @return The current coefficient of the Taylor term
        */
        double coef() const;

        /**
            Reset the coefficient of the current Taylor term.
        */
        void
        reset();

        /**
            Return a string representation of the term in the index
            representation.  See #indices_.
            @return A descriptive string for the Taylor term
        */
        std::string
        name() const;

        /**
            Print a description of the taylor term
            @param os The output tream to print to
        */
        void
        print(std::ostream& os = std::cout) const;

        DerivativePtr deriv() const;

        /**
            Filter out all the nonzero contributions to a given displacement.

            @param terms The set of all terms
            @param nozero_terms The vector in which all nonzero_terms will be placed
            @param disp The displacement to get the nonzero_terms frmo
            @param tol The tolerance for declaring a displacement to be zero
        */
        static void
        getNonzeroTerms(
            const std::vector<TaylorTermPtr>& terms,
            std::vector<TaylorTermPtr>& nonzero_terms,
            ConstDisplacementPtr disp,
            double tol = 1e-6
        );

        /**
            Generate the set of all unique Taylor terms for a range of derivatives.

            @param terms The vector to append the terms to
            @param start_depth The lowest level of derivative
            @param stop_depth The highest level of derivative
            @param ncoords The number of coordinates
        */
        static void
        generateTerms(
            std::vector<TaylorTermPtr>& terms,
            int start_depth,
            int stop_depth,
            int ncoords
        );

        /**
            Recursion method for generating terms.  See generateTerms. 

            @param terms The vector to contain the created terms
            @param current_depth The number of spots that have been filled.  For example,
                                 if doing 3rd derivatives, you might have already placed
                                 all index possibilities in the first spot, and will have
                                 two index spots left to place.
            @param maxdepth The total number of index spots to fill.
            @param ncoords The number of coordinates
        */
        static void
        appendTerms(
            std::vector<std::vector<int> >& terms,
            int current_depth,
            int maxdepth,
            int ncoords
        );

};

class TaylorSeriesEnergy : public smartptr::Countable {

    private:
        DisplacementPtr disp_;

        std::vector<TaylorTermPtr> terms_;

        double energy_;

        TaylorSeriesEnergy(
            const std::vector<TaylorTermPtr>& terms,
            DisplacementPtr disp,
            double E_0
        );

    public:
        static void
        buildEnergyApproximations(
            std::vector<TaylorSeriesEnergyPtr>& energies,
            DisplacementIteratorPtr dispiter,
            DerivativeIteratorPtr deriter
        );

        void print(std::ostream& os = std::cout) const;

        double get_energy() const;

        double error() const;

};

} //end namespace

#endif
