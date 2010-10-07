#include <config.h>
#include <src/displacement.h>
#include <src/keyword.h>
#include <src/exception.h>
#include <src/fileio.h>
#include <src/defines.h>
#include <src/taylor.h>
#include <src/permutation.h>
#include <src/symmetry.h>
#include <src/derivative.h>
#include <src/utilities.h>

#define ZERODISP 1e-3


using namespace gigide;
using namespace smartptr;
using namespace std;


SerialDeclare(Displacement);
SerialDeclare(DisplacementIterator);
SerialDeclare(DisplacementMapping);

DisplacementIterator::DisplacementIterator(
    const Set<ConstInternalCoordinatePtr>& coords, 
    const Set<ConstSimpleInternalCoordinatePtr>& simples,
    const ConstMoleculePtr& mol,
    const std::vector<double> dispsizes
) :
    mol_(mol), simples_(simples), coords_(coords), dispsizes_(dispsizes)
{
    SetRuntime(DisplacementIterator);

    int nderiv = KeywordSet::getKeyword("nderiv")->getValueInteger();
    int nvalue = KeywordSet::getKeyword("nvalue")->getValueInteger();
    this->deriv_level_ = nderiv - nvalue;

    int ncomp = nderiv - nvalue;
    int nrobust = KeywordSet::getKeyword("robustness")->getValueInteger() - 1;
    TaylorTerm::generateTerms(
        terms_,
        0,
        ncomp + nrobust, //validate robustness
        coords_.size()
    );

    nunique_ = 0;

    bool usedispfile = KeywordSet::getKeyword("usedispfile")->getValueBoolean();
    if (!usedispfile)
        return;

    except("Dispfile not yet implemented.  Please turn off usedispfile");

    string basedir = GIGDATADIR;
    switch(deriv_level_)
    {
        case 0:
            basedir += "/derivs_0th";
            break;
        case 1:
            basedir += "/derivs_1st";
            break;
        case 2:
            basedir += "/derivs_2nd";
            break;
        case 3:
            basedir += "/derivs_3rd";
            break;
        case 4:
            basedir += "/derivs_4th";
            break;
        default:
            except(stream_printf("Invalid derivative level %d. Please try again.\n", deriv_level_));
    }
    readDispFile(basedir);
}

DisplacementIterator::DisplacementIterator(const XMLArchivePtr& arch)
    : Serializable(arch)
{
    SetRuntime(DisplacementIterator);
    serial_load(deriv_level);
    serial_load(nunique);
    serial_load(mol);
    serial_load(displacements);
    serial_load(coords);
    serial_load(simples);
}

void
DisplacementIterator::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(deriv_level);
    serial_save(nunique);
    serial_save(mol);
    serial_save(displacements);
    serial_save(coords);
    serial_save(simples);
}

DisplacementIterator::const_iterator
DisplacementIterator::begin() const
{
    return displacements_.begin();
}

DisplacementIterator::iterator
DisplacementIterator::begin()
{
    return displacements_.begin();
}

DisplacementIterator::const_iterator
DisplacementIterator::end() const
{
    return displacements_.end();
}

DisplacementIterator::iterator
DisplacementIterator::end()
{
    return displacements_.end();
}

void
DisplacementIterator::print(ostream& os) const
{
    os << stream_printf("Displacement iterator with %d displacements", displacements_.size());
    for (const_iterator it(begin()); it != end(); ++it)
    {
        os << endl;
        (*it)->print(os);
    }
}

void
DisplacementIterator::readDispFile(
    string dispfile
)
{
    vector< vector<double> > short_disp_list;
    readDisplacements(dispfile.c_str(), short_disp_list);

    //from the short list, generate the full displacmements by filling in zeros
    for (int dispnum=0; dispnum < short_disp_list.size(); dispnum++)
    {
        vector<double> disp = short_disp_list[dispnum];
        //only construct a set of displacements from this if the lenght does
        //not exceed the number of coordinates
        if ( disp.size() <= coords_.size() )
        {
            vector<double> full_disp = disp; //this copies
            //put in the zeros of the displacement
            for (int i=disp.size(); i < coords_.size(); i++)
                full_disp.push_back(0);
			
            vector< vector<double> > disp_combos;
            PermutationGeneratorDoublePtr permgen( new PermutationGenerator<double>(full_disp, coords_.size()));
            permgen->generatePermutations(disp_combos);

            for (int dispnum=0; dispnum < disp_combos.size(); dispnum++)
            {
                vector<double> disp_sizes = disp_combos[dispnum];
                DisplacementPtr new_disp = new Displacement(disp_sizes, mol_, coords_, simples_, dispsizes_);
                displacements_.push_back(new_disp);
            }
        }
    }
}

void
DisplacementIterator::addDisplacement(DisplacementPtr disp)
{
    disp->assignTerms(terms_);
    displacements_.push_back(disp);
}

void
DisplacementIterator::readDisplacements(
    const char* filename,
    vector< vector<double> >& disp_vec
)
{
    ifstream myfile(filename);
    string line;
    if (myfile.is_open())
    {
        while (!myfile.eof())
        {
            getline(myfile, line);
            stringstream sstr(line);
            double i;
            vector<double> disp;
            while (sstr >> i)
            {
                disp.push_back(i);
            }
            disp_vec.push_back(disp);
        }
    }
    //delete the last... not sure why it adds an extra one
    disp_vec.pop_back();
}

DisplacementPtr
DisplacementIterator::getDisplacement(vector<double>& disps)
{
    DisplacementPtr disp = findDisplacement(disps);
    if (disp.get() != NULL)
        return disp;

    disp = new Displacement(disps, mol_, coords_, simples_, dispsizes_);

    addDisplacement(disp);
    return disp;
}

DisplacementPtr
DisplacementIterator::findDisplacement(const vector<double>& disps)
{
    //iterate the map to find the displacement
    for (iterator it(begin()); it != end(); ++it)
    {
        DisplacementPtr disp = *it;
        if (disp->matches(disps))
        {
            return disp;
        }
    }
    //nothing found
    return NULL;
}

int 
DisplacementIterator::ncoord() const
{
    return coords_.size();
}

int
DisplacementIterator::nunique() const
{

    return nunique_;
}

int
DisplacementIterator::ndisps() const
{
    return displacements_.size();
}

int
DisplacementIterator::level() const
{
    return deriv_level_;
}

DisplacementPtr
DisplacementIterator::findZeroDisplacement()
{
    vector<double> zeros(coords_.size(), 0);
    DisplacementPtr check = findDisplacement(zeros);
    return check;
}

bool
DisplacementIterator::hasZeroDisplacement()
{
    return findZeroDisplacement().get() != NULL;
}

void
DisplacementIterator::addZeroDisplacement()
{
    vector<double> zeros(coords_.size(), 0);
    DisplacementPtr new_disp = new Displacement(zeros, mol_, coords_, simples_, dispsizes_);
    displacements_.push_back(new_disp);
}

void
DisplacementIterator::mapEquivalentDisplacements()
{
    int ntot = displacements_.size();
    int n = 1;
    for (iterator it(begin()); it != end(); ++it, ++n)
    {
        cout << stream_printf("Finding equivalent displacements for %d out of %d", n, ntot) << endl;
        DisplacementPtr disp = *it;
        mapEquivalentDisplacements(disp);
    }

    nunique_ = 0;
    vector<DisplacementPtr>::const_iterator it;
    for (it = displacements_.begin(); it != displacements_.end(); ++it)
    {
        if ( (*it)->isUnique() )
            ++nunique_;
    }
}

bool
DisplacementIterator::mapEquivalentDisplacements(
    const DisplacementPtr& refdisp
)
{
	int debug = KeywordSet::getKeyword("symmetry debug")->getValueInteger();

    if (debug >= 1)
    {
        cout << "Mapping equivalent displacements for "; refdisp->print();
    }

    vector<ConstSymmetryOperationPtr> symm_elems;
    vector<DisplacementPtr>::iterator it;
    mol_->getSymmetryElements(symm_elems);
    RectMatrixPtr symm_op, perm_op;

    bool foundmatch = false;

    //now map all equivalent displacements
    for (it = displacements_.begin(); it != displacements_.end(); ++it)
    {
        DisplacementPtr disp = *it;
        if (debug >= 2)
        {
            cout << "\tTesting equivalence to "; disp->print();
        }
        
        bool equiv = testEquivalence(refdisp, disp, symm_elems, perm_op, symm_op);
        if (equiv)
        {
            refdisp->addEquivalentDisplacement(disp, symm_op, perm_op);
            foundmatch = true;
        }
    }

    //if we have gone through all the displacements and found no match, send back null
    return foundmatch;
}

bool
DisplacementIterator::testEquivalence(
    const ConstDisplacementPtr& refdisp,
    const ConstDisplacementPtr& disp,
    const Set<ConstSymmetryOperationPtr>& symm_elems,
    RectMatrixPtr& perm_op,
    RectMatrixPtr& symm_op
)
{
    bool degree_mismatch = fabs(disp->degree() - refdisp->degree()) > 1e-8; 
    bool type_mismatch = disp->disptype() != refdisp->disptype();
    bool mag_mismatch = fabs(disp->dispmag() - refdisp->dispmag()) > 1e-8;

	int debug = KeywordSet::getKeyword("symmetry debug")->getValueInteger();

    int exp = KeywordSet::getKeyword("symmetry tolerance")->getValueInteger();
    double tol = pow(10, -exp);

#if 0
    vector<double> incs(refdisp->getIncrements()); 
    bool printstuff = (fabs(incs[7] + 2.0) < 1e-8) && (fabs(incs[8] + 1.0) < 1e-8);
    if (printstuff)
    {
        refdisp->print();
        debug = 2;
    }
#endif

    if (debug >= 2)
    {
        disp->print();
        cout << stream_printf("\tDegree mismatch=%d Type mismatch=%d Magnitude mismatch=%d",
                              degree_mismatch,
                              type_mismatch,
                              mag_mismatch) << endl;
    }

    if (degree_mismatch || type_mismatch || mag_mismatch)
        return false;

    if (refdisp == disp) //don't map same displacements
        return false;



    ConstRectMatrixPtr refBmat = refdisp->getDisplacementMatrix(); 
    ConstRectMatrixPtr dispBmat = disp->getDisplacementMatrix(); 
    //loop through the symmetry operations to see if any of them work
    for (int opnum=0; opnum < symm_elems.size(); opnum++)
    {
        ConstSymmetryOperationPtr symm_elem = symm_elems[opnum];
        if (symm_elem->type() == SymmetryOperation::identity)
            continue;  //don't use identity

        if (debug >= 2)
            cout << "\t\tTesting equivalence through " << symm_elem->name() << endl;

        RectMatrixPtr permutation = symm_elem->getPermutationMatrix(mol_);
        RectMatrixPtr temp = refBmat * symm_elem->matrix();
        RectMatrixPtr transBmat = permutation * temp;

        //if the displacement matrices for the two are the same, return the displacement
        if (equals(dispBmat, transBmat, tol))
        {
            if (debug)
                cout << "\t\tFound equivalent displacement!" << endl;

            perm_op = permutation.copy();
            symm_op = symm_elem->matrix().copy();

            return true;
        }
    }

    return false;
}


vector<double>
DisplacementIterator::computeIncrements(
    ConstRectMatrixPtr xyz
)
{
    MoleculePtr dispmol = mol_->copy();
    dispmol->setXYZ(xyz);
    vector<double> increments;
    for (int i=0; i < coords_.size(); i++)
    {
        ConstInternalCoordinatePtr coord = coords_[i];
        double disp = dispsizes_[i];
        double increment = computeIncrement(dispmol, coord, disp);
        increments.push_back(increment);
    }
    
    int debug = KeywordSet::getKeyword("geometry debug")->getValueInteger();
    if (debug >= 2)
    {
        cout << "Center molecule" << endl;
        cout << mol_->getXYZString() << endl;

        cout << "Displacement molecule" << endl;
        cout << dispmol->getXYZString() << endl << endl << endl;
    }
    

    return increments;
}

double
DisplacementIterator::computeIncrement(const ConstMoleculePtr& mol, const ConstInternalCoordinatePtr& refcoord, double disp)
{
    double refval = refcoord->getValue();
    double newval = refcoord->getValueForMolecule(mol);
    double dispnum = (newval - refval) / disp;

    int debug = KeywordSet::getKeyword("geometry debug")->getValueInteger();
    if (debug >= 2)
    {
        cout << stream_printf("Reference value: %18.14f Displacement value: %18.14f Displacement: %18.14e",
                              refval,
                              newval,
                              dispnum) << endl;
    }

    return dispnum;

}

DisplacementIterator::~DisplacementIterator()
{
    //have the displacement mappings clear themselves
    DisplacementIterator::iterator it;
    for (it = begin(); it != end(); ++it)
        (*it)->clearMappings();
}

Displacement::Displacement(const vector<double>& increments, 
                           const ConstMoleculePtr& mol,
                           const Set<ConstInternalCoordinatePtr>& coords,
                           const Set<ConstSimpleInternalCoordinatePtr>& simples,
                           const std::vector<double> dispsizes
)
  : mol_(mol), coords_(coords), simples_(simples), increments_(increments), dispsizes_(dispsizes)
{
    SetRuntime(Displacement);

    energy_assigned_ = false;
    grad_assigned_ = false;
    fc_assigned_ = false;
    energy_ = 0.0;
    refcount_ = 0;
    unique_ = true;

    init();
}

Displacement::Displacement(const ConstMoleculePtr& mol,
                           const Set<ConstInternalCoordinatePtr>& coords,
                           const Set<ConstSimpleInternalCoordinatePtr>& simples,
                           const std::vector<double> dispsizes
)
    : mol_(mol), coords_(coords), simples_(simples), increments_(coords.size(), 0), dispsizes_(dispsizes)
{
    SetRuntime(Displacement);

    energy_assigned_ = false;
    grad_assigned_ = false;
    fc_assigned_ = false;
    energy_ = 0.0;
    refcount_ = 0;
    unique_ = true;

    init();
}

Displacement::Displacement(const XMLArchivePtr& arch)
    : Serializable(arch)
{
    SetRuntime(Displacement);
    serial_load(energy);
    serial_load(energy_assigned);
    serial_load(grad_assigned);
    serial_load(fc_assigned);
    serial_load(unique);
    serial_load(refcount);
    serial_load(degree);
    serial_load(dispmag);
    serial_load(disptype);
    serial_load(dispsizes);
    serial_load(gradients);
    serial_load(fc);
    serial_load(increments);
    serial_load(dispmol);
    serial_load(displacement_matrix);
    serial_load(mol);
    serial_load(equiv_disps);
    serial_load(coords);
    serial_load(simples);
    serial_load(label);
}

void
Displacement::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);

    serial_save(energy);
    serial_save(energy_assigned);
    serial_save(grad_assigned);
    serial_save(fc_assigned);
    serial_save(unique);
    serial_save(refcount);
    serial_save(degree);
    serial_save(dispmag);
    serial_save(disptype);
    serial_save(dispsizes);
    serial_save(gradients);
    serial_save(fc);
    serial_save(increments);
    serial_save(dispmol);
    serial_save(displacement_matrix);
    serial_save(mol);
    serial_save(equiv_disps);
    serial_save(coords);
    serial_save(simples);
    serial_save(label);
}
        

void
Displacement::init()
{
    label_ = label(increments_);

    //build the "displacement" matrix for this displacments
    //this just corresponds to the B vectors scaled by the displacement extent
    //start from the first coordinates b matrix
    displacement_matrix_ = coords_[0]->getBMatrix() * increments_[0];
    for (int i=1; i < coords_.size(); i++)
    {
        double increment = increments_[i];
        if ( fabs(increment) > 1e-6 ){
            ConstRectMatrixPtr disp = coords_[i]->getBMatrix();
            displacement_matrix_.accumulate(disp * increment);
        }
    }

    dispmag_ = 0;
    mdim_t size = displacement_matrix_.nrow() * displacement_matrix_.ncol();
    const double* block = displacement_matrix_.data();
    const double* blockptr = block;
    for (mdim_t i=0; i < size; ++i, ++blockptr)
    {
        double val = *blockptr;
        dispmag_ += val * val;
    }

    degree_ = 0;
    vector<double>::const_iterator it;
    vector<double> nonzero_disps;
    for (it = increments_.begin(); it != increments_.end(); ++it)
    {
        double val = *it;
        degree_ += fabs(val);

        if (val != 0)
            nonzero_disps.push_back(fabs(val));
    }

    sort(nonzero_disps.begin(), nonzero_disps.end());
    
    stringstream sstr;
    for (it = nonzero_disps.begin(); it != nonzero_disps.end(); ++it)
        sstr << stream_printf(" %4.2f", *it);
    disptype_ = sstr.str();

    if (!dispsizes_.size())
        dispsizes_ = GigideKeyword::getDisplacementSizes(coords_.size());
}

void
Displacement::validate() const
{
    int debug = KeywordSet::getKeyword("geometry debug")->getValueInteger();
    double tol = pow(10, -1 * KeywordSet::getKeyword("geometry tolerance")->getValueInteger());

    if (debug)
        cout << stream_printf("\t        %20s   %20s   %20s", "Reference", "Displacement","Step") << endl;  

    for (int i=0; i < coords_.size(); ++i)
    {
        ConstInternalCoordinatePtr coord = coords_[i];
        double val = coord->getValueForMolecule(dispmol_);
        double refval = coord->getValue();
        double absdisp = val - refval;
        double disp = absdisp / dispsizes_[i];

        if (debug)
            cout << stream_printf("\td[%d] -> %20.14f   %20.14f   %20.14f", i, refval, val, disp) << endl;

        //this should match displacements
        double error = fabs(disp - increments_[i]) * dispsizes_[i]; //extra * is for renormalization
        if (error > tol)
        {
            cout << stream_printf("Dispsize = %8.4e Displacement = %18.12f", dispsizes_[i], increments_[i]) << endl;
            cerr << stream_printf("\t        %20s   %20s   %20s", "Reference", "Displacement","Step") << endl;  
            cerr << stream_printf("\td[%d] -> %20.14f   %20.14f   %20.14f", i, refval, val, absdisp) << endl;
            cerr << stream_printf("Tolerance: %8.4e  Error: %8.4e", tol, error) << endl;
            stringstream sstr;
            sstr << "Generated displacement is not valid" << endl;
            sstr << label() << endl;
            except(sstr.str());
        }
    }
}


void
Displacement::clearMappings()
{
    equiv_disps_.clear();
}

bool
Displacement::matches(
    const vector<double>& increments
) const
{
    if (increments.size() != increments_.size())
        return false;

    vector<double>::const_iterator itme, itother;
    for (itme = increments_.begin(), itother = increments.begin();
         itme != increments_.end();
         ++itme, ++itother)
    {
        if ( fabs(*itme - *itother) > ZERODISP ) //diff too big
        {
            return false;
        }
    }

    //all match
    return true;
}

string
Displacement::label() const
{
    return label_;
}

bool
Displacement::energyAssigned() const
{
    return energy_assigned_;
}

Set<ConstInternalCoordinatePtr> 
Displacement::coords() const
{
    return coords_;
}

vector<double>
Displacement::getIncrements() const
{
    return increments_;
}

double
Displacement::getEnergy() const
{
    if (!energy_assigned_)
        except("No energy assigned yet to displacement. Cannot getEnergy.");

    return energy_;
}

ConstRectMatrixPtr
Displacement::getDisplacementMatrix() const
{
    return displacement_matrix_;
}

ConstMoleculePtr 
Displacement::getDisplacementMolecule() const
{
    return dispmol_;
}

void
Displacement::resetTaylorTerms()
{
    vector<TaylorTermPtr >::iterator it;
    for (it = terms_.begin(); it != terms_.end(); ++it)
        (*it)->reset();
}

void
Displacement::addNonzeroTerms(
    map<string, TaylorTermPtr>& terms,
    double tol
) const
{
    vector<TaylorTermPtr>::const_iterator it;
    for (it = terms_.begin(); it != terms_.end(); ++it)
    {
        TaylorTermPtr term = *it;
        double coef = term->coef();
        if ( fabs(coef) > tol )
        {
            terms[term->name()] = term;
        }
    }
}

void
Displacement::accumulateTerms(double coef, int maxlevel)
{
    vector<TaylorTermPtr >::const_iterator it;
    for (it = terms_.begin(); it != terms_.end(); ++it)
    {
        if ( (*it)->level() <= maxlevel )
            (*it)->accumulate(this, coef);
    }
}

void
Displacement::assignTerms(
    vector<TaylorTermPtr >& terms
)
{
    TaylorTerm::getNonzeroTerms(
        terms,
        terms_,
        this
    );
}

string
Displacement::taylorRepresentation() const
{
    if (terms_.size() == 0)
        return "";

    vector<TaylorTermPtr>::const_iterator it;
    stringstream sstr;
    sstr << terms_[0]->polynomialString(this);
    for (it = terms_.begin() + 1; it != terms_.end(); ++it)
    {
        sstr << " + " << (*it)->polynomialString(this);
    }
    return sstr.str();
}

double
Displacement::increment(int i) const
{
    return increments_[i];
}

double
Displacement::displacement(int i) const
{
    return dispsizes_[i] * increments_[i];
}

void
Displacement::print(ostream& os, bool include_e) const
{
    os << "[";
    int i=0;
    for ( ; i < increments_.size() - 1; ++i)
    {
         os << stream_printf("%5.2f, ", increments_[i]);
    }
    os << stream_printf("%5.2f]", increments_[i]);
    if (energy_assigned_ && include_e)
        os << stream_printf(" E = %18.12f", energy_);
}

bool
Displacement::isAssigned() const
{
    return energy_assigned_ || fc_assigned_ || grad_assigned_; 
}

bool
Displacement::isUnique() const
{
    return unique_;
}

void
Displacement::setUnique(bool uniq)
{
    //my original status cannot be overwritten
    //I have subtasks
    if (equiv_disps_.size() != 0)
        return; 

    unique_ = uniq;
}


double
Displacement::dispmag() const
{
    return dispmag_;
}

string
Displacement::disptype() const
{
    return disptype_;
}

double
Displacement::degree() const
{
    return degree_;
}

string
Displacement::label(const vector<double>& disps)
{
    stringstream sstr;
    sstr << "<" << stream_printf("%5.3f", disps[0]);
    for (int i=1; i < disps.size(); i++)
        sstr << stream_printf(",%5.3f", disps[i]);
    sstr << ">";
    return sstr.str();
}

int
Displacement::getRefcount() const
{
    return refcount_;
}

void
Displacement::incrementRefcount()
{
    refcount_++;
}

void
Displacement::decrementRefcount()
{
    refcount_--;
}


void
Displacement::assignEnergy(double energy)
{
    if (energy_assigned_)
        return; //done


    //set the energy for this displacement
    energy_ = energy;
    energy_assigned_ = true;

    cout << stream_printf("\tAssigned energy to %s", label().c_str()) << endl;

    //pass the information along to all equivalent displacements
    for (int i=0; i < equiv_disps_.size(); i++)
        equiv_disps_[i]->assignEnergy(energy);
}

bool
Displacement::isZeroDisplacement() const
{
    for (int i=0; i < increments_.size(); i++)
    {
        if (increments_[i] != 0)
            return false;
    }
    return true;
}

void
Displacement::assignGradients(ConstRectMatrixPtr xyzgradients)
{
    if (grad_assigned_)
        return;
    
    gradients_ = IntderWriter::transformGradientsToInternals(xyzgradients, dispmol_, simples_, coords_);
    grad_assigned_ = true;

    for (int i=0; i < equiv_disps_.size(); i++)
    {
        equiv_disps_[i]->assignGradients(xyzgradients);
    }

}

void
Displacement::assignForceConstants(ConstRectMatrixPtr xyzfc, ConstRectMatrixPtr xyzgrads)
{
    if (fc_assigned_)
        return;

    fc_ = IntderWriter::transformForceConstantsToInternals(xyzfc, xyzgrads, dispmol_, simples_, coords_);
    fc_assigned_ = true;

    for (int i=0; i < equiv_disps_.size(); i++)
    {
        equiv_disps_[i]->assignForceConstants(xyzfc, xyzgrads);
    }

}

double
Displacement::getDerivativeValue(const ConstDerivativePtr& deriv) const
{
	int debug = KeywordSet::getKeyword("fit debug")->getValueInteger();
    if (debug >= 2)
    {
        cout << "Reading values from displacement ";
        print();
        if (gradients_.nonnull()) gradients_.print("Gradients");
        if (fc_.nonnull()) fc_.print("Force constants");
    }

    vector<int> indices = deriv->indices();
    if (indices.size() == 0) //this just needs the energy
    {
        return energy_;
    }
    else if (indices.size() == 1) //needs a gradient
    {
        int index = indices[0];
        if (gradients_.null())
        {
            except("Gradients not assigned to necessary point. Cannot continue.");
        }
        return gradients_.get_element(index);
    }
    else if (indices.size() == 2) //needs a force constant
    {
        int row = indices[0];
        int col = indices[1];

        if (fc_.null())
        {
            except("Force constants not assigned to necessary point. Cannot continue.");
        }

        return fc_.get_element(row, col);
    }

    return 0; //get rid of annoything compiler flags
}

void
Displacement::generateDisplacement()
{
    vector<double> disp = getIncrements();
    vector<double> displacements;
    for (int i=0; i < disp.size(); i++)
    {
        double increment = disp[i];
        double dispsize = dispsizes_[i];
        double totaldisp = increment * dispsize;
        displacements.push_back(totaldisp);
    }
    dispmol_ = displaceGeometry(displacements, mol_, coords_, simples_);
}

MoleculePtr
Displacement::displaceGeometry(
    const vector<double>& displacements,
    const ConstMoleculePtr& mol,
    const Set<ConstInternalCoordinatePtr>& coords,
    const Set<ConstSimpleInternalCoordinatePtr>& simples
)
{
    vector<InternalCoordinatePtr> blank1;
    vector<SimpleInternalCoordinatePtr> blank2;
    return displaceGeometry(displacements, mol, coords, simples, blank1, blank2);
}

MoleculePtr
Displacement::displaceGeometry(
    const vector<double>& displacements,
    const ConstMoleculePtr& mol,
    const Set<ConstInternalCoordinatePtr>& coords,
    const Set<ConstSimpleInternalCoordinatePtr>& simples,
    std::vector<InternalCoordinatePtr>& newcoords,
    std::vector<SimpleInternalCoordinatePtr>& newsimples
)
{
    vector<double> desired_values;
    int num_internals = coords.size();
    //figure out the desired values of the internal coordinates that we want
    for (int coordnum=0; coordnum < displacements.size(); coordnum++)
    {
        double disp = displacements[coordnum];
        ConstInternalCoordinatePtr coord = coords[coordnum];
        double current_value = coord->getValue();
        double desired_value = current_value + disp;
        desired_values.push_back(desired_value);
    }
    return generateGeometry(desired_values, mol, coords, simples, newcoords, newsimples);
}
MoleculePtr
Displacement::generateGeometry(
    const vector<double>& desired_values, 
    const ConstMoleculePtr& mol,
    const Set<ConstInternalCoordinatePtr>& main_coords,
    const Set<ConstSimpleInternalCoordinatePtr>& main_simples
)
{
    vector<InternalCoordinatePtr> blank1;
    vector<SimpleInternalCoordinatePtr> blank2;
    return generateGeometry(desired_values, mol, main_coords, main_simples, blank1, blank2);
}

MoleculePtr
Displacement::generateGeometry(
    const vector<double>& desired_values, 
    const ConstMoleculePtr& mol,
    const Set<ConstInternalCoordinatePtr>& main_coords,
    const Set<ConstSimpleInternalCoordinatePtr>& main_simples,
    std::vector<InternalCoordinatePtr>& coords,
    std::vector<SimpleInternalCoordinatePtr>& simples
)
{
    //create the new molecule
    MoleculePtr dispmol = mol->copy();

    //create the copies to be used in generating the geometry
    for (int i=0; i < main_coords.size(); ++i)
        simples.push_back(main_simples[i]->simple_copy(dispmol));

    for (int i=0; i < main_coords.size(); ++i)
        coords.push_back(main_coords[i]->copy(dispmol, simples));


    int num_internals = coords.size();

    int tolerance_exp = KeywordSet::getKeyword("geometry tolerance")->getValueInteger();
    double tolerance = pow(10, -tolerance_exp);
    int maxiter = KeywordSet::getKeyword("geometry maxit")->getValueInteger();
    int print_frequency = KeywordSet::getKeyword("geometry print frequency")->getValueInteger();
    int debug = KeywordSet::getKeyword("geometry debug")->getValueInteger();

    //this is the total displacement summed over all coordinates
    RectMatrixPtr total_disp_matrix(dispmol->natoms(), 3);

    //check the coordinates to make sure the values are valid
    for (int coordnum=0; coordnum < num_internals; coordnum++)
    {
        ConstInternalCoordinatePtr coord = coords[coordnum];
        double newval = desired_values[coordnum];
        if (!coord->isValidValue(newval))
        {
            stringstream sstr;
            sstr << stream_printf("%12.8f is not a valid value for the given coordinate", newval) << endl;
            coord->print(sstr);
            except(sstr.str());
        }
    }

    //now start the iterations until the molecule has converged to what we want
    double error = tolerance + 1;
    int iteration = 0;
    while (error > tolerance)
    {
    if (debug >= 2 && (iteration % print_frequency == 0) )
        cout << dispmol->getXYZString() << endl;
        //zero the disp array for the current step
        total_disp_matrix.zero();
        error = 0.0;
        for (int coordnum=0; coordnum < num_internals; coordnum++)
        {
            ConstInternalCoordinatePtr coord = coords[coordnum];
            //this is the bvector displacement for just the given coordinate
            RectMatrixPtr coord_disp_matrix = total_disp_matrix.clone(); 
            //figure out how far away the current value
            //is from the desired value
            double desired_value = coord->canonicalizeValue(desired_values[coordnum]);
            double current_value = coord->getValue();
            double diff = desired_value - current_value;
            //add that difference to the error
            error += fabs(diff);
            //this will be the shift induced by the bvectors;
            double shift = 0.0;
            if (debug >= 2 && (iteration % print_frequency == 0) )
                cout << stream_printf("%20s", coord->type().c_str())
                     << " Current Value: " << stream_printf("%18.14f", current_value)
                     << " Desired Value: " << stream_printf("%18.14f", desired_value)
                     << " Error:         " << stream_printf("%8.4e", diff) << endl;
            for (int natom=0; natom < dispmol->natoms(); natom++)
            {
                ConstAtomPtr next_atom = dispmol->getAtom(natom);
                //get the bvector for this coordinate for this atom
                ConstVectorPtr bvec = coord->getBVector(next_atom);
                //how far we displace along the bvector depends on how far away
                //we are from the desired value
                double bvec_mag = bvec.norm();
                //this is how much the current bvector will change the internal coordinate
                double bvec_shift = bvec_mag * bvec_mag;
                //add that contribution to the total shift
                shift += bvec_shift;
                //assign the bvector to the complete displacement matrix
                coord_disp_matrix.assign_row(bvec, natom);
            }
            if (debug >= 3 && (iteration % print_frequency == 0) )
                coord_disp_matrix.print("B Matrix");


            if ( fabs(shift) < 1e-8 )
            {
                stringstream err;
                err << stream_printf("Coordinate has B matrix magnitude %12.6e.  Cannot normalize displacement.", shift) << endl;
                coord->printDetail(err);
                except(err.str());
            }

            //scale the coordinate displacement by the appropriate amount
            double scale_factor = diff / shift * SHIFT_DAMP; //test what happens with a "damping" vector
            if (debug >= 5 && (iteration % print_frequency == 0) && fabs(diff) > 1e-2)
            {
                cout << stream_printf("Difference: %12.6e", diff) << endl;
                cout << stream_printf("B-Vector Magnitude: %12.6e", shift) << endl;
                cout << stream_printf("Scale factor: %12.6e", scale_factor) << endl;
                coord->printDetail();
            }

            coord_disp_matrix.scale(scale_factor);
            //and accumulate the displacement for this coordinate into the total displacement
            total_disp_matrix.accumulate(coord_disp_matrix);
        }
        if (debug >= 3 && (iteration % print_frequency == 0) )
            total_disp_matrix.print("Displacement");
        dispmol->displace(total_disp_matrix);
        //now, we should go back and have the coordinates recompute themselves with the updated xyz coordinates
        for (int i=0; i < simples.size(); i++)
            simples[i]->recompute();
        for (int i=0; i < coords.size(); i++)
            coords[i]->recompute();
        //and increment the iteration number
        iteration++;
        if (debug && (iteration % print_frequency == 0) )
            cout << "Iteration:  " << stream_printf("%4d", iteration) << " Error: " << stream_printf("%8.4e", error) << " Tolerance: " << stream_printf("%8.4e", tolerance)  << endl;
        
        if (iteration >= maxiter)
        {
            stringstream sstr;
            sstr << stream_printf("Maximum number of iterations exceeded. Final error %8.4e", error) << endl;
            sstr << "Final geometry at end of iterations:" << endl;
            sstr << dispmol->getXYZString() << endl;
            except(sstr.str());
        }
    }
    if (debug)
        cout << endl << endl;

    //recompute point group
    dispmol->recomputePointGroup(mol);

    return dispmol;
}

void
Displacement::addEquivalentDisplacement(
    const DisplacementPtr& equiv, 
    ConstRectMatrixPtr symm_op, 
    ConstRectMatrixPtr perm_op
)
{
    equiv->setUnique(false);

    vector<DisplacementMappingPtr>::iterator it;
    for (it = equiv_disps_.begin(); it != equiv_disps_.end(); ++it)
    {
        if ( (*it)->getDisplacement() == equiv )
            return; //we already have this 
    }
    
    //we don't have this
    DisplacementMappingPtr map = new DisplacementMapping(equiv, symm_op, perm_op);  
    equiv_disps_.push_back(map);
}

void
Displacement::printEquivalentDisplacements(ostream& os) const
{
    os << "Equivalent displacements for "; print(os);
    vector<DisplacementMappingPtr>::const_iterator it;
    for (it = equiv_disps_.begin(); it != equiv_disps_.end(); ++it)
    {
        os << "\t"; (*it)->getDisplacement()->print(os);
    }
    cout << endl;
}

DisplacementMapping::DisplacementMapping(
    const DisplacementPtr& disp, 
    ConstRectMatrixPtr symmop, 
    ConstRectMatrixPtr permop
) : disp_(disp), symmop_(symmop), permop_(permop)
{
    SetRuntime(DisplacementMapping);
}

DisplacementMapping::DisplacementMapping(const XMLArchivePtr& arch)
{
    SetRuntime(DisplacementMapping);
    serial_load(disp);
    serial_load(symmop);
    serial_load(permop);
}

void
DisplacementMapping::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(disp);
    serial_save(symmop);
    serial_save(permop);
}

ConstDisplacementPtr
DisplacementMapping::getDisplacement() const
{
    return disp_;
}

void 
DisplacementMapping::assignEnergy(double energy)
{
    disp_->assignEnergy(energy);
}

void
DisplacementMapping::assignGradients(ConstRectMatrixPtr xyzgradients)
{
    RectMatrixPtr newmat = xyzgradients * symmop_;
    RectMatrixPtr newgrads = permop_ * newmat;
    disp_->assignGradients(newgrads);
}

void
DisplacementMapping::assignForceConstants(ConstRectMatrixPtr xyzfc, ConstRectMatrixPtr xyzgrads)
{
    //construct the transformation matrix based on permop and symmop
    int natoms = xyzfc.nrow() / 3;
    RectMatrixPtr tform(xyzfc.nrow(), xyzfc.nrow());

    int rowoffset = 0;
    for (int i=0; i < natoms; ++i, rowoffset += 3)
    {
        int coloffset = 0;
        for (int j=0; j < natoms; ++j, coloffset += 3)
        {
            double permval = permop_.get_element(j,i); //need transpose

            if ( fabs(permval) < 1e-8 ) //nothing to include
                continue;

            //assign the symmetry operation as the subblock
            tform.assign_subblock(symmop_, rowoffset, coloffset);
        }
    }

    RectMatrixPtr transformed_fc = tform.t() * xyzfc * tform;
    RectMatrixPtr transformed_grads = permop_ * xyzgrads * symmop_;
    disp_->assignForceConstants(transformed_fc, transformed_grads);
}

