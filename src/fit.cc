
#include <src/qff.h>
#include <src/fit.h>
#include <src/taylor.h>
#include <src/utilities.h>
#include <src/derivative.h>
#include <src/displacement.h>
#include <src/exception.h>

using namespace gigide;
using namespace std;
using namespace smartptr;

SerialDeclare(FitPoint);
SerialDeclare(Fit);

class IndexI {
    public:
        static string symbol(){return "i";}
};

class IndexJ {
    public:
        static string symbol(){return "j";}
};

class IndexK {
    public:
        static string symbol(){return "k";}
};

class IndexL {
    public:
        static string symbol(){return "l";}
};

template<
    class FitType,
    class I
>
class CreateFit_i : public FitFactory {

    public:
        FitPtr getFit(const ConstDerivativePtr& deriv, const DisplacementIteratorPtr& iter);

};

template<
    class FitType,
    class I,
    class J
>
class CreateFit_ij : public FitFactory {

    public:
        FitPtr getFit(const ConstDerivativePtr& deriv, const DisplacementIteratorPtr& iter);

};

template<
    class FitType,
    class I,
    class J,
    class K
>
class CreateFit_ijk : public FitFactory {

    public:
        FitPtr getFit(const ConstDerivativePtr& deriv, const DisplacementIteratorPtr& iter);

};

template<
    class FitType,
    class I,
    class J,
    class K,
    class L
>
class CreateFit_ijkl : public FitFactory {

    public:
        FitPtr getFit(const ConstDerivativePtr& deriv, const DisplacementIteratorPtr& iter);

};

#define CONHEAD(clsname) public Fit { public: clsname(
#define CONTAIL , const DisplacementIteratorPtr& iter); };
#define FIT1(clsname) class clsname : CONHEAD(clsname) int i CONTAIL
#define FIT2(clsname) class clsname : CONHEAD(clsname) int i, int j CONTAIL
#define FIT3(clsname) class clsname : CONHEAD(clsname) int i, int j, int k CONTAIL
#define FIT4(clsname) class clsname : CONHEAD(clsname) int i, int j, int k, int l CONTAIL

FIT1(Fit_i);
FIT1(Fit_i_d3); //3rd order robustness
FIT1(Fit_ii);
FIT1(Fit_ii_d3);
FIT1(Fit_iii);
FIT1(Fit_iiii);
FIT2(Fit_ij);
FIT2(Fit_ij_d3);
FIT2(Fit_iij);
FIT2(Fit_iiij);
FIT2(Fit_iijj);
FIT3(Fit_ijk);
FIT3(Fit_iijk);
FIT4(Fit_ijkl);

FitPoint::FitPoint(
    const DisplacementPtr& disp,
    double coef
) : disp_(disp), coef_(coef)
{
    SetRuntime(FitPoint);
    disp_->incrementRefcount();
}

FitPoint::FitPoint(const XMLArchivePtr& arch)
{
    SetRuntime(FitPoint);
    serial_load(coef);
    serial_load(disp);
}

void
FitPoint::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(coef);
    serial_save(disp);
}

double
FitPoint::contribution(
    const ConstDerivativePtr& value
) const
{
    double derval = disp_->getDerivativeValue(value) * coef_;

    return derval;
}

void
FitPoint::addNonzeroTerms(
    map<string, TaylorTermPtr >& terms
) const
{
    disp_->addNonzeroTerms(terms);
}

void
FitPoint::reset()
{
    disp_->resetTaylorTerms();
}

void
FitPoint::accumulateTerms(int maxlevel)
{
    disp_->accumulateTerms(coef_, maxlevel);
}

string
FitPoint::name() const
{
    stringstream sstr;
    int num,denom;
    sstr << stream_printf("%18.14f", coef_);
    sstr << disp_->label();
    return sstr.str();
}

FitPoint::~FitPoint()
{
    disp_->decrementRefcount();
}

Fit::Fit() //blank for use with subclasses
{
    SetRuntime(Fit);
    init_statics();
}

Fit::Fit(
    DisplacementIteratorPtr iter,
    VectorPtr coefs
)
{
    SetRuntime(Fit);
    init_statics();
    init(iter, coefs);
}

Fit::Fit(const XMLArchivePtr& arch)
{
    SetRuntime(Fit);
    serial_load(points);
}

void
Fit::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(points);
}

void
Fit::init(
    DisplacementIteratorPtr iter,
    VectorPtr coefs
)
{
    int n = 0;
    for (DisplacementIterator::iterator it(iter->begin()); it != iter->begin(); ++it, ++n)
    {
        double coef = coefs.get_element(n);
        if ( fabs(coef) < 1e-6 ) //do not include
            continue;

        addPoint(*it, coef);
    }
}

void
Fit::addPoint(
    DisplacementPtr disp,
    double coef
)
{
    FitPointPtr fitpoint = new FitPoint(disp, coef);
    points_.push_back(fitpoint);
}

double
Fit::computeDerivative(
    ConstDerivativePtr value
) const
{
    double val = 0;
    vector<FitPointPtr>::const_iterator it;
    for (it = points_.begin(); it != points_.end(); ++it)
    {
        val += (*it)->contribution(value);
    }
    return val;
}

void
Fit::print(ostream& os) const
{
    if (points_.size() == 0) //no points
        return;

    stringstream sstr;
    vector<FitPointPtr>::const_iterator it;
    sstr << points_[0]->name();
    for (it = points_.begin() + 1; it != points_.end(); ++it)
        sstr << " + " << (*it)->name();

    os << sstr.str() << endl;
}

void
Fit::init_statics()
{
    if (initdone_)
        return;

    //second order robustness formulae
    factories_[2]["i"] = new CreateFit_i<Fit_i,IndexI>;
    factories_[2]["ii"] = new CreateFit_i<Fit_ii,IndexI>;
    factories_[2]["iii"] = new CreateFit_i<Fit_iii,IndexI>;
    factories_[2]["iiii"] = new CreateFit_i<Fit_iiii,IndexI>;
    factories_[2]["ij"] = new CreateFit_ij<Fit_ij,IndexI,IndexJ>;
    factories_[2]["iij"] = new CreateFit_ij<Fit_iij,IndexI,IndexJ>;
    factories_[2]["ijj"] = new CreateFit_ij<Fit_iij,IndexJ,IndexI>;
    factories_[2]["iiij"] = new CreateFit_ij<Fit_iiij,IndexI,IndexJ>;
    factories_[2]["ijjj"] = new CreateFit_ij<Fit_iiij,IndexJ,IndexI>;
    factories_[2]["iijj"] = new CreateFit_ij<Fit_iijj,IndexI,IndexJ>;
    factories_[2]["ijk"] = new CreateFit_ijk<Fit_ijk,IndexI,IndexJ,IndexK>;
    factories_[2]["iijk"] = new CreateFit_ijk<Fit_iijk,IndexI,IndexJ,IndexK>;
    factories_[2]["ijjk"] = new CreateFit_ijk<Fit_iijk,IndexJ,IndexI,IndexK>;
    factories_[2]["ijkk"] = new CreateFit_ijk<Fit_iijk,IndexK,IndexI,IndexJ>;
    factories_[2]["ijkl"] = new CreateFit_ijkl<Fit_ijkl,IndexI,IndexJ,IndexK,IndexL>;

    //third order robustness formulae
    factories_[3]["i"] = new CreateFit_i<Fit_i_d3,IndexI>;
    factories_[3]["ii"] = new CreateFit_i<Fit_ii_d3,IndexI>;
    factories_[3]["ij"] = new CreateFit_ij<Fit_ij_d3,IndexI,IndexJ>;

    initdone_ = true;
}

void
Fit::reset()
{
    vector<FitPointPtr>::iterator it;
    for (it = points_.begin(); it != points_.end(); ++it)
        (*it)->reset();
}

void
Fit::computeTerms(int maxlevel)
{
    vector<FitPointPtr>::iterator it;
    for (it = points_.begin(); it != points_.end(); ++it)
        (*it)->accumulateTerms(maxlevel);
}

void
Fit::printNonzeroTerms(ostream& os) const
{
    vector<FitPointPtr>::const_iterator it;
    map<string, TaylorTermPtr > terms;
    for (it = points_.begin(); it != points_.end(); ++it)
        (*it)->addNonzeroTerms(terms);

    map<string, TaylorTermPtr >::iterator itt;
    for (itt = terms.begin(); itt != terms.end(); ++itt)
        itt->second->print(os);
}

void
Fit::assignFit(
    DerivativePtr deriv,
    DisplacementIteratorPtr iter,
    int robustness
)
{
    string typelabel = deriv->typelabel();
    if (typelabel == "") //just assign a null fit, this is the center displacement
    {
        deriv->assignFit(new Fit);
        return;
    }

    FitPtr fit =  factories_[robustness][typelabel]->getFit(deriv, iter);
    if (fit.get() == NULL)
    {
        stringstream sstr;
        sstr << "No fitting function for derivative for robustness " << robustness << endl;
        deriv->printDetail(sstr);
        except(sstr.str().c_str());
    }
    deriv->assignFit(fit);
}

#include <src/i.fit>
#include <src/i_d3.fit>
#include <src/ii.fit>
#include <src/ii_d3.fit>
#include <src/iii.fit>
#include <src/iiii.fit>
#include <src/ij.fit>
#include <src/ij_d3.fit>
#include <src/iij.fit>
#include <src/iiij.fit>
#include <src/iijj.fit>
#include <src/ijk.fit>
#include <src/iijk.fit>
#include <src/ijkl.fit>

template <
    class FitType,
    class I
>
FitPtr 
CreateFit_i<FitType,I>::getFit(const ConstDerivativePtr& deriv, const DisplacementIteratorPtr& iter)
{
    return new FitType(
                deriv->index(I::symbol()), 
                iter);
}

template <
    class FitType,
    class I,
    class J
>
FitPtr 
CreateFit_ij<FitType,I,J>::getFit(const ConstDerivativePtr& deriv, const DisplacementIteratorPtr& iter)
{
    return new FitType(
        deriv->index(I::symbol()), 
        deriv->index(J::symbol()),
        iter);
}

template <
    class FitType,
    class I,
    class J,
    class K
>
FitPtr 
CreateFit_ijk<FitType,I,J,K>::getFit(const ConstDerivativePtr& deriv, const DisplacementIteratorPtr& iter)
{
    return new FitType(
                deriv->index(I::symbol()), 
                deriv->index(J::symbol()), 
                deriv->index(K::symbol()), 
                iter);
}

template <
    class FitType,
    class I,
    class J,
    class K,
    class L
>
FitPtr 
CreateFit_ijkl<FitType,I,J,K,L>::getFit(const ConstDerivativePtr& deriv, const DisplacementIteratorPtr& iter)
{
    return new FitType(
                deriv->index(I::symbol()), 
                deriv->index(J::symbol()), 
                deriv->index(K::symbol()), 
                deriv->index(L::symbol()), 
                iter);
}

bool Fit::initdone_ = false;
map<int, map<string, FitFactoryPtr> > Fit::factories_;
