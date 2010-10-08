#include <src/derivative.h>
#include <src/symmetry.h>
#include <src/utilities.h>
#include <src/permutation.h>
#include <src/exception.h>
#include <src/molecule.h>
#include <src/timer.h>

using namespace gigide;
using namespace std;
using namespace smartptr;

SerialDeclare(Derivative);
SerialDeclare(DerivativeIterator);

DerivativeIterator::DerivativeIterator(
    int level,
    const Set<ConstInternalCoordinatePtr>& coords,
    const ConstMoleculePtr& mol
) : coords_(coords), mol_(mol)
{
    SetRuntime(DerivativeIterator);

    this->deriv_level_ = level;
    this->coords_ = coords;
    mol_ = mol;
    
    vector< vector<int> > deriv_combinations;
    for (int i=0; i <= deriv_level_; ++i)
    {
        generateDerivativeCombos(i, deriv_combinations);
    }

    //we have now generated the "quanta" of the derivatives, i.e.
    //the vector <1,0,1,0> says we have a 1 quanta of 1, 1 quanta of 3
    //which would be the mixed derivative d4E/d1x d3y, single derivative
    //which respect to x and a third derivative with respect to y
    //now we should go through and translate these "quanta" into actual derivatives
    for (int num=0; num < deriv_combinations.size(); num++)
    {
        //this only holds non-zero quanta
        vector<int> quanta_set = deriv_combinations[num];
        //this holds all quanta, including zero for the uninvolved
        vector<int> all_quanta;
        //add the non-zero quanta
        int numcoords_done=0;
        int total_quanta=0;
        for (int quanta_level=0; quanta_level < quanta_set.size(); ++quanta_level)
        {
            int num_diff_quanta = quanta_set[quanta_level];
            for (int i=0; i < num_diff_quanta; ++i)
            {
                all_quanta.push_back(quanta_level + 1);
                numcoords_done++;
            }
            total_quanta += num_diff_quanta;
        }
        //and finally fill out with zeros;
        for (int quanta=numcoords_done; quanta < coords_.size(); ++quanta)
        {
            all_quanta.push_back(0);
        }
        
        //we only want to include this derivative if the number of different indices
        //does not exceed the number of coordinates... i.e. we can't compute a mixed fourth derivative
        //of three internal coordinates
        if (total_quanta <= coords_.size())
        {
            vector< vector<int> > final_combos;
            PermutationGeneratorIntPtr permgen(new PermutationGenerator<int>(all_quanta, coords_.size()));
            permgen->generatePermutations(final_combos);

            for (int comb_num=0; comb_num < final_combos.size(); comb_num++)
            {
                vector<int> combo = final_combos[comb_num];
                DerivativePtr new_deriv = new Derivative(combo, coords_, mol_);
                derivs_.push_back(new_deriv);
                //we want a way of accessing derivatives by their label
                string label = Derivative::label(combo);
            }
        }
    }
}

DerivativePtr
DerivativeIterator::getDerivativeFromIndices(const vector<int>& indices) const
{
    vector<int> levels(coords_.size(), 0);
    vector<int>::const_iterator it(indices.begin());
    for ( ; it != indices.end(); ++it)
        ++levels[*it];

    return getDerivative(levels);
}

DerivativePtr
DerivativeIterator::getDerivative(const vector<int>& levels) const
{
    //iterate the map to find the displacement
    for (vector<DerivativePtr>::const_iterator it = derivs_.begin(); it != derivs_.end(); ++it)
    {
        DerivativePtr deriv = *it;
        if (deriv->matches(levels))
        {
            return deriv;
        }
    }
    //nothing found
    return NULL;
}

DerivativeIterator::DerivativeIterator(const XMLArchivePtr& arch)
    : Serializable(arch)
{
    SetRuntime(DerivativeIterator);

    serial_load(deriv_level);

    serial_load(nonzero_level);

    serial_load(coords);

    serial_load(derivs);

    serial_load(mol);
}

int
DerivativeIterator::ncoords() const
{
    return coords_.size();
}

void
DerivativeIterator::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(deriv_level);
    serial_save(nonzero_level);
    serial_save(coords);
    serial_save(derivs);
    serial_save(mol);
}

DerivativeIterator::~DerivativeIterator()
{
}

DerivativeIterator::const_iterator
DerivativeIterator::begin() const
{
    return derivs_.begin();
}

DerivativeIterator::iterator
DerivativeIterator::begin()
{
    return derivs_.begin();
}

DerivativeIterator::const_iterator
DerivativeIterator::end() const
{
    return derivs_.end();
}

DerivativeIterator::iterator
DerivativeIterator::end()
{
    return derivs_.end();
}

int
DerivativeIterator::level() const
{
    return deriv_level_;
}

DerivativePtr
DerivativeIterator::findCombination(const ConstDerivativePtr& deriv, const ConstDerivativePtr& value) const
{
    vector<int> deriv_levels = deriv->levels();
    vector<int> value_levels = value->levels();
    
    if (deriv_levels.size() != value_levels.size())
    {
        except("Invalid derivative combination.  Different number of internal coordinates for each.");
    }

    vector<int> combo;
    for (int i=0; i < deriv_levels.size(); i++)
        combo.push_back(deriv_levels[i] + value_levels[i]);

    DerivativePtr combo_deriv = getDerivative(combo);
    return combo_deriv;
}

void
DerivativeIterator::printDetail(ostream& os) const
{
    for (const_iterator it(begin()); it != end(); ++it)
        (*it)->printDetail(os);
}

void
DerivativeIterator::generateDerivativeCombos(int level, vector< vector<int> >& all_combos) const
{
    //begin with a "blank" combination
    vector<int> current_combo(level, 0);
    recurseDerivativeCombos(level, 0, 0, all_combos, current_combo);
}

void
DerivativeIterator::recurseDerivativeCombos(int level, int quanta, int index_level,
                                            vector< vector<int> >& combos, vector<int> current_combo) const
{
    if      (quanta == level) //we have filled out the current level
    {
        //push back apparently pushes back by reference, so we must copy first
        combos.push_back(current_combo);
    }
    else if (index_level >= level) //we have run over the number of allowed indices
        return;
    else if (quanta > level) //we have overshot
        return;
    else if (quanta < level) //keep going on the recursion
    {
        int index = 0;
        int current_quanta = quanta;
        while (current_quanta <= level)
        {
            current_quanta = quanta + index*(index_level+1);
            current_combo[index_level] = index; 
            recurseDerivativeCombos(level, current_quanta, index_level+1, combos, current_combo);
            index++;
        }
    }
}

int 
DerivativeIterator::nderivs() const
{
    return derivs_.size();
}

string
Derivative::label() const
{
    return label(levels_);
}

string
Derivative::label(const vector<int>& levels)
{
    ostringstream label;
    //add the first index
    if (levels.size())
        label << stream_printf("%2d", levels[0]);
    //add the rest
    for (int i=1; i < levels.size(); i++)
        label << stream_printf(" %2d", levels[i]);
    return label.str();
}

string
Derivative::polynomialString() const
{
    stringstream sstr;
    for (int i=0; i < levels_.size(); ++i)
    {
        int count = levels_[i];
        if (count == 1)
            sstr << stream_printf("x%d", i);
        else if (count > 1)
            sstr << stream_printf("(x%d)^%d", i, count);
    }
    return sstr.str();
}

Derivative::Derivative(
    const vector<int>& levels, 
    const Set<ConstInternalCoordinatePtr>& coords,
    const ConstMoleculePtr& mol
) :
    levels_(levels), coords_(coords), mol_(mol)
{
    SetRuntime(Derivative);
    level_ = 0;
    //the number of indices involved, i.e. you can have a 4th derivative with 3 indices, d4E/dxi2 dxj dxk
    int num_indices = 0;

    //in computing the permutational symmetry of the derivative, we need to keep track of double counting
    int prefactor = 1;
    for (int coord=0; coord < levels.size(); coord++)
    {
        int coord_level = levels[coord];
        level_ += coord_level;
        if (coord_level) 
        {
            num_indices++;
            prefactor *= factorial(coord_level);
        }
    }

    //the taylor coefficient is 1/n! * the permutational symmetry number of the derivative
    coeff_ = 1.0 / prefactor;

    //for how many quanta of this coordinate we have, multiply character array
    bool nonzero_abelian = true;
    Set<PointGroupPtr> subgroups = mol_->getPointGroup()->getSubgroups();
    for (int s=0; s < subgroups.size(); ++s)
    {
        PointGroupPtr subpg = subgroups[s];
        VectorPtr subchars(subpg->nclasses());
        for (int i=0; i < subchars.n(); ++i)
        {
            double elem = 1.0;
            for (int num=0; num < coords.size(); ++num)
            {
                ConstInternalCoordinatePtr coord = coords[num];
                double nextchar = coord->subgroup_characters(subpg->name())[i];
                int level = levels[num];
                for (int l=0; l < level; ++l)
                     elem *= nextchar;
            }
            subchars.set_element(i, elem);
        }
        bool nonzero = subpg->isTotallySymmetric(subchars);
        if (!nonzero) //subgroup symmety makes this zero
        {
            nonzero_abelian = false;
            break;
        }
    }

    VectorPtr fullchars(mol_->getPointGroup()->nclasses());
    for (int i=0; i < fullchars.n(); ++i)
    {
        double elem = 1.0;
        for (int num=0; num < coords.size(); ++num)
        {
            ConstInternalCoordinatePtr coord = coords[num];
            double nextchar = coord->character(i);
            int level = levels[num];
            for (int l=0; l < level; ++l)
                 elem *= nextchar;
        }
        fullchars.set_element(i, elem);
    }
    bool nonzero_full = mol_->getPointGroup()->isTotallySymmetric(fullchars);

    //at this point, we have completed the character multiplications for each coordinate
    //and we have them stored in char_array. See if this char_array belong to the totally symmetric irrep
    //only ignore if we explicitly say to ignore zero derivs
    bool keepzeros = !(KeywordSet::getKeyword("usezero")->getValueBoolean());
    nonzero_ = (nonzero_abelian && nonzero_full) || keepzeros;


    init_statics();

    //compute the type indices
    vector<int>::iterator it;
    vector<string>::iterator itstr = letters_.begin();
    stringstream sstr;
    string letter = *itstr;
    int n = 0;
    for (it = levels_.begin(); it != levels_.end(); ++it, ++n)
    {
        int count = *it;
        for (int c=0; c < count; ++c)
            sstr << letter;
        
        if (count > 0)
        {
            indices_[letter] = n;
            ++itstr;
            letter = *itstr;

            nonzero_levels_[n] = count;
        }
    }
    typelabel_ = sstr.str();

    dertype_ = typemap_[typelabel_];

    unique_ = true;
    assigned_ = false;
}

Derivative::Derivative(const XMLArchivePtr& arch)
    : Serializable(arch)
{
    SetRuntime(Derivative);

    serial_load(coeff);
    serial_load(level);
    serial_load(nonzero);
    serial_load(assigned);
    serial_load(unique);
    serial_load(typelabel);
    serial_load_enum(dertype);
    serial_load(levels);
    serial_load(nonzero_levels);
    serial_load(indices);
    serial_load(mol);
    serial_load(fit);
    serial_load(coords);
    serial_load(equiv_derivs);
    serial_load(deriv_value);

    init_statics();
}

void
Derivative::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(coeff);
    serial_save(level);
    serial_save(nonzero);
    serial_save(assigned);
    serial_save(unique);
    serial_save(typelabel);
    serial_save_enum(dertype);
    serial_save(levels);
    serial_save(nonzero_levels);
    serial_save(indices);
    serial_save(mol);
    serial_save(fit);
    serial_save(coords);
    serial_save(equiv_derivs);
    serial_save(deriv_value);
}

void
Derivative::addEquivalentDerivative(
    DerivativePtr equiv
)
{
    equiv->setUnique(false);

    equiv_derivs_.push_back(equiv);
}

void
Derivative::mapEquivalentDerivatives(
    ConstDerivativeIteratorPtr iter
)
{
    //only diagonal force constants are supported right now
    if (nonzero_levels_.size() != 1)
        return;

    int num = nonzero_levels_.begin()->first;
    int level = nonzero_levels_.begin()->second;

    vector<ConstInternalCoordinatePtr> equiv_coords = coords_[num]->getDegenerateCoordinates();

    for (int n=0; n < equiv_coords.size(); ++n)
    {
        int coordnum = InternalCoordinate::getCoordinateNumber(equiv_coords[n], coords_);
        vector<int> levels(coords_.size(), 0);
        levels[coordnum] = level;
        DerivativePtr deriv = iter->getDerivative(levels);
        addEquivalentDerivative(deriv);
    }

}

void
Derivative::printDetail(ostream& os) const
{
    for (int i=0; i < levels_.size(); ++i)
    {
        int quanta = levels_[i];
        if (quanta == 0)
            continue;
        os << "Quanta " << quanta << ": "; coords_[i]->print(os);
    }
    cout << stream_printf("Nonzero? %d", nonzero_) << endl;
    os << stream_printf("Value = %14.8f", deriv_value_) << endl;
}

void
Derivative::setUnique(bool uniq)
{
    //I am a parent and cannot be overridden
    if (equiv_derivs_.size() != 0)
        return;

    unique_ = uniq;
}

bool
Derivative::isUnique() const
{
    return unique_;
}

Derivative::DerivativeType
Derivative::dertype() const
{
    return dertype_;
}

void
Derivative::setValue(double val)
{
    if (assigned_) return;

    deriv_value_ = val;
    assigned_ = true;

    for (int i=0; i < equiv_derivs_.size(); ++i)
        equiv_derivs_[i]->setValue(val);

}

int
Derivative::level() const
{
    return level_;
}

double
Derivative::value() const
{
    if (!assigned_)
        except("Derivative accessed before value assigned");

    return deriv_value_;
}

int
Derivative::index(const std::string& str) const
{
    map<string, int>::const_iterator it = indices_.find(str);
    if (it == indices_.end())
        return 0;

    return it->second;
}

void
Derivative::init_statics()
{
    if (initdone_)
        return;

    letters_.push_back("i");
    letters_.push_back("j");
    letters_.push_back("k");
    letters_.push_back("l");
    letters_.push_back("m");
    letters_.push_back("n");
    letters_.push_back("o");
    letters_.push_back("p");

    typemap_["i"] = Dertype_ii;
    typemap_["ii"] = Dertype_ii;
    typemap_["iii"] = Dertype_iii;
    typemap_["iiii"] = Dertype_iiii;
    typemap_["ij"] = Dertype_ij;
    typemap_["iij"] = Dertype_iij;
    typemap_["ijj"] = Dertype_iij;
    typemap_["iijj"] = Dertype_iijj;
    typemap_["ijjj"] = Dertype_iiij;
    typemap_["iiij"] = Dertype_iiij;
    typemap_["ijk"] = Dertype_ijk;
    typemap_["iijk"] = Dertype_iijk;
    typemap_["ijjk"] = Dertype_iijk;
    typemap_["ijkk"] = Dertype_iijk;
    typemap_["ijkl"] = Dertype_ijkl;

    initdone_ = true;
}

string
Derivative::typelabel() const
{
    return typelabel_;
}

void
Derivative::validateFit(int robustlevel)
{
    //if robust = 2, we want to check +1 terms
    //if robust = 3, we want to check +1,+2 terms
    //robust = 2 means first polluting terms come two orders higher
    //which means we expect +1 order derivatives to be zero
    //robust = 3 means first polluting terms come three orders higher
    //which means we expect +1,+2 order derivatives to be zero
    int extraterms = robustlevel - 1;
    cout << "Validating fit for "; print();
    fit_->reset();
    fit_->computeTerms(level() + extraterms);
    fit_->printNonzeroTerms();

    bool printfit = KeywordSet::getKeyword("printfit")->getValueBoolean();
    if (printfit)
        fit_->print();
}

double
Derivative::compute(ConstDerivativePtr value) const
{
    return fit_->computeDerivative(value);
}

void
Derivative::assignFit(
    FitPtr fit
)
{
    fit_ = fit;
}

void
Derivative::print(ostream& os) const
{
    os << label();
    if (!nonzero_)
        os << " = 0";
    os << endl;
}

vector<int>
Derivative::indices() const
{
    vector<int> deriv_indices;
    for (int coord=0; coord < levels_.size(); coord++)
    {
        for (int deg=0; deg < levels_[coord]; deg++)
        {
            deriv_indices.push_back(coord); 
        }
    }
    return deriv_indices;
}

double
Derivative::taylorCoeff(vector<double>& disps) const
{
    double newcoeff = coeff_;
    for (int i=0; i < levels_.size(); ++i)
    {
        int deriv_level = levels_[i];
        double dispsize = disps[i];
        for (int j=0; j < deriv_level; ++j)
            newcoeff *= dispsize;
    }

    return newcoeff;
}

double
Derivative::displacementCoeff(
    vector<double>& disps
) const
{
    double newcoeff = 1.0;
    for (int i=0; i < levels_.size(); ++i)
    {
        int deriv_level = levels_[i];
        double dispsize = disps[i];
        for (int j=0; j < deriv_level; ++j)
            newcoeff *= dispsize;
    }

    return newcoeff;
}

bool 
Derivative::nonzero() const
{
    return nonzero_;
}

bool
Derivative::zero() const
{
    return !nonzero_;
}

vector<int>
Derivative::levels() const
{
    return levels_;
}

bool
Derivative::matches(const vector<int>& levels) const
{
    if (levels.size() != levels_.size()) return false;

    for (int i=0;  i < levels.size(); ++i)
    {
        if (levels[i] != levels_[i])
            return false;
    }

    return true;
}

