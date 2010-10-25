
#ifndef gigide_fit_h
#define gigide_fit_h

#include <vector>
#include <src/deftypes.h>

#include <src/smartptr/src/printstream.h>

namespace gigide {

class TaylorTerm;
class DerivativeIterator;
class Derivative;
class Displacement;
class DisplacementIterator;

class FitPoint : public smartptr::Serializable {

    private:
        double coef_;
        DisplacementPtr disp_;

    public:
        FitPoint(
            const DisplacementPtr& disp,
            double coef
        );

        FitPoint(
            const XMLArchivePtr& parser
        );

        void serialize(const XMLArchivePtr& writer) const;

        double
        contribution(
            const ConstDerivativePtr& value
        ) const;

        void
        addNonzeroTerms(std::map<std::string, TaylorTermPtr>& terms) const;

        void
        accumulateTerms(int maxlevel);

        std::string name() const;

        void
        reset();

        ~FitPoint();
        
};

class Fit;
class FitFactory : public smartptr::Countable {

    public:
        virtual FitPtr getFit(const ConstDerivativePtr& deriv, const DisplacementIteratorPtr& iter) = 0;
};

/**
    @class Fit

    A class representing a derivative fit to a set of finite displacements
*/
class Fit : public smartptr::Serializable {
    
    private:
        std::vector<FitPointPtr> points_;
        static std::map<int, std::map<std::string, FitFactoryPtr> > factories_;
        static bool initdone_;
        static void init_statics();

    protected:

        void addPoint(
            DisplacementPtr disp,
            double coef
        );

    public:
        Fit(
            DisplacementIteratorPtr iter,
            VectorPtr coefs
        );

        Fit();

        Fit(
            const XMLArchivePtr& parser
        );

        void serialize(const XMLArchivePtr& writer) const;

        void init(
            DisplacementIteratorPtr iter,
            VectorPtr coefs
        );

        void
        reset();

        void
        print(std::ostream& os = std::cout) const;

        double
        computeDerivative(
            ConstDerivativePtr value
        ) const;

        /**
            Using the default, internally defined fitting functions, assign a fitting
            function to the derivative.
            @param deriv The derivative to assign the fit to
            @param iter The displacement iterator that will provide the finite displacement points
            @param robustness The level of contaminating derivative contributions that will be zero. For example,
                              if fitting 2nd derivatives and robustness is 2, any contaminating contribution
                              from 3rd and 4th derivatives will be zero.  However, 5th derivatives will still
                              contribute.
            @throw GigideException If no default fit has been programmed for requested derivative and robustness
        */
        static void
        assignFit(
            DerivativePtr deriv,
            DisplacementIteratorPtr iter,
            int robustness
        );

        void
        printNonzeroTerms(std::ostream& os = std::cout) const;

        void
        computeTerms(int maxlevel);
};

} //end namespace

#endif
