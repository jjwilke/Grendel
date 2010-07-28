
#include <src/taylor.h>
#include <src/utilities.h>
#include <src/qff.h>
#include <src/displacement.h>
#include <src/derivative.h>

using namespace gigide;
using namespace std;

TaylorTerm::TaylorTerm(
    const vector<int>& indices
)
    : indices_(indices), coef_(0)
{
    vector<int> pvec;
    int lastindex = -1;
    int count = 1;
    for (int i=0; i < indices_.size(); ++i)
    {
        int nextindex = indices_[i];
        if (nextindex == lastindex)
        {
            count++;
        }
        else
        {
            pvec.push_back(count);
            count = 1;
        }

        lastindex = nextindex;
    }
    //add the last one
    pvec.push_back(count);

    prefactor_ = 1.0;
    for (int i=0; i < pvec.size(); ++i)
    {
        prefactor_ /= factorial(pvec[i]);
    }
}

string
TaylorTerm::polynomialString(
    ConstDisplacementPtr disp
) const
{
    stringstream sstr;
    int num,denom;
    double c = coef(disp);
    bool found = getFraction(c,num,denom);
    if (found)
        sstr << stream_printf("%d/%d",num,denom);
    else
        sstr << stream_printf("%8.4f",c);

    for (int i=0; i < indices_.size(); ++i)
    {
        int index = indices_[i];
        sstr << stream_printf("x%d", index);
    }
    return sstr.str();
}

double
TaylorTerm::coef(
    ConstDisplacementPtr disp
) const
{
    double coef = prefactor_;
    vector<int>::const_iterator it;
    for (it = indices_.begin(); it != indices_.end(); ++it)
        coef *= disp->displacement(*it);
    return coef;
}

void
TaylorTerm::getNonzeroTerms(
    vector<TaylorTermPtr>& terms,
    vector<TaylorTermPtr>& nonzero_terms,
    ConstDisplacementPtr disp,
    double tol
)
{
    vector<TaylorTermPtr >::iterator it;
    for (it = terms.begin(); it != terms.end(); ++it)
    {
        TaylorTermPtr term = *it;
        double coef = term->coef(disp);
        if ( fabs(coef) > tol )
            nonzero_terms.push_back(term);
    }
}

void
TaylorTerm::generateTerms(
    std::vector< TaylorTermPtr >& terms,
    int start_depth,
    int stop_depth,
    int ncoords
)
{
    vector< vector<int> >::iterator it;

    //build the taylor terms for all derivatives > 0
    for (int n=start_depth; n <= stop_depth; ++n)
    {
        vector< vector<int> > vecs;
        appendTerms(vecs, 0, n, ncoords);
        
        for (it = vecs.begin(); it != vecs.end(); ++it)
            terms.push_back(new TaylorTerm(*it));
    }
}

double
TaylorTerm::coef() const
{
    return coef_ * prefactor_;
}

void
TaylorTerm::appendTerms(
    vector<vector<int> >& terms,
    int current_depth,
    int maxdepth,
    int ncoords
)
{
    if (current_depth == maxdepth) return;

    if (current_depth == 0) //no terms yet
    {
        for (int n=0; n < ncoords; ++n)
        {
            vector<int> vec; vec.push_back(n);
            terms.push_back(vec);
        }

    }
    else
    {
        vector< vector<int> >::iterator it;

        vector< vector<int> > newvec;
        for (int n=0; n < ncoords; ++n)
        {
            vector< vector<int> > copy = terms;
            for (it = copy.begin(); it != copy.end(); ++it)
            {
                vector<int> next = *it;
                bool noncanonical = false;
                for (int i=0; i < next.size(); ++i)
                {
                    if (next[i] > n) //non-canonical
                        noncanonical = true;
                }

                if (noncanonical)
                {
                    continue;
                }
                else
                {
                    next.push_back(n);
                    newvec.push_back(next);
                }
            }
        }

        terms = newvec; //transfer terms
    }
    //now descend the ranks
    appendTerms(terms, current_depth + 1, maxdepth, ncoords);
}

void
TaylorTerm::accumulate(
    ConstDisplacementPtr disp,
    double coef
)
{
    vector<int>::iterator it;
    int polynomial = 1;
    for (it = indices_.begin(); it != indices_.end(); ++it)
    {
        polynomial *= disp->displacement(*it);
    }

    coef_ += coef * polynomial;
}

void
TaylorTerm::print(ostream& os) const
{
    stringstream sstr;
    sstr << stream_printf("%4.2f T(", coef());
    if (indices_.size() > 0)
    {
        sstr << indices_[0];
        for (int i=1; i < indices_.size(); ++i)
            sstr << stream_printf(",%d", indices_[i]);
    }
    sstr << ")";

    os << sstr.str() << endl;
}

string
TaylorTerm::name() const
{
    stringstream sstr;
    sstr << "T(";
    if (indices_.size() > 0)
    {
        sstr << indices_[0];
        for (int i=1; i < indices_.size(); ++i)
            sstr << stream_printf(",%d", indices_[i]);
    }
    sstr << ")";

    return sstr.str();
}

int
TaylorTerm::level() const
{
    return indices_.size();
}

void
TaylorTerm::reset()
{
    coef_ = 0;
}
