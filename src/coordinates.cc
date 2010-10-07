#include <src/coordinates.h>
#include <src/coordinates.h>
#include <src/molecule.h>
#include <src/symmetry.h>
#include <src/utilities.h>
#include <src/exception.h>
#include <math.h>

#include <src/smartptr/src/ref.h>

#define X 0
#define Y 1
#define Z 2

using namespace std;
using namespace gigide;
using namespace smartptr;

SerialDeclare(BondLength);
SerialDeclare(BondAngle);
SerialDeclare(OutOfPlaneBend);
SerialDeclare(Torsion);
SerialDeclare(LinX);
SerialDeclare(LinY);
SerialDeclare(Lin1);
SerialDeclare(SymmetryInternalCoordinate);
SerialDeclare(CoordinateSubspace);

InternalCoordinate::InternalCoordinate(const ConstMoleculePtr& mol, gigstr name)
    : mol_(mol), coordtype_(name)
{
}

InternalCoordinate::InternalCoordinate(const XMLArchivePtr& arch)
    : Serializable(arch)
{
    serial_load(mol);
    serial_load(bmatrix);
    serial_load(chars);
    serial_load(coordtype);
}

InternalCoordinate::~InternalCoordinate()
{
}

void
InternalCoordinate::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(mol);
    serial_save(bmatrix);
    serial_save(chars);
    serial_save(coordtype);
    serial_save(subgroup_chars);
}

bool
SimpleInternalCoordinate::matches(const ConstSimpleInternalCoordinatePtr& coord) const
{
    return coordtype_ == coord->type() && connectivityString() == coord->connectivityString();
}

string
InternalCoordinate::connectivityString() const
{
    return "";
}

string
InternalCoordinate::type() const
{
    return coordtype_;
}

bool
InternalCoordinate::isValidValue(double val) const
{
    return true;
}

double
InternalCoordinate::canonicalizeValue(double val) const
{
    return val;
}

ConstRectMatrixPtr
InternalCoordinate::getBMatrix() const
{
    return bmatrix_;
}

ConstVectorPtr
InternalCoordinate::characters() const
{
    return chars_;
}

double
InternalCoordinate::character(int cls) const
{
    return chars_[cls];
}

ConstVectorPtr
InternalCoordinate::subgroup_characters(gigstr subgroup) const
{
    map<string, VectorPtr>::const_iterator it = subgroup_chars_.find(subgroup);

    if (it == subgroup_chars_.end())
        except(stream_printf("%s is not a valid subgroup", subgroup.c_str()));

    return it->second;
}

void
InternalCoordinate::addDegenerateCoordinate(const ConstInternalCoordinatePtr& coord)
{
    equiv_coords_.push_back(coord);
}

void
InternalCoordinate::print(ostream& os, bool includechar, bool includesubchars) const
{
    string connect = connectivityString();
    os << stream_printf("%8s %10s %18.12f", coordtype_.c_str(), connect.c_str(), getConvertedValue());
    if (includechar)
    {
        os << " Characters = [";
        for (int i=0; i < chars_.n(); ++i)
            os << stream_printf("%8.2f", chars_[i]);
        os << " ]";
    }

    if (includesubchars)
    {
        map<string, VectorPtr>::const_iterator it;
        for (it = subgroup_chars_.begin(); it != subgroup_chars_.end(); ++it)
        {
            os << endl << stream_printf("\t%5s: [", (it->first).c_str());
            VectorPtr v = it->second;
            for (int i=0; i < v.n(); ++i)
                os << stream_printf("%8.2f", v[i]);
            os << " ]";
        }
    }
}

void
InternalCoordinate::init()
{
    int natoms = mol_->natoms();
    bmatrix_ = RectMatrixPtr(natoms, XYZ_DIM); 

    recompute();

    if (chars_.null()) //might be read in by XML constructor
        chars_ = new Vector(mol_->getPointGroup()->nclasses());
}

bool
InternalCoordinate::isEquivalent(const ConstInternalCoordinatePtr& coord) const
{
    double exp = KeywordSet::getKeyword("symmetry tolerance")->getValueDouble();
    double tol = pow(10, -exp);

    ConstRectMatrixPtr myBmat = getBMatrix();
    ConstRectMatrixPtr otherBmat = coord->getBMatrix();
    vector<ConstSymmetryOperationPtr> symmops;
    vector<ConstSymmetryOperationPtr>::const_iterator it;
    mol_->getSymmetryElements(symmops);
    for (it = symmops.begin(); it != symmops.end(); ++it)
    {
        ConstSymmetryOperationPtr symmop = *it;
        RectMatrixPtr permop = symmop->getPermutationMatrix(mol_); 
        RectMatrixPtr transBmat = permop * otherBmat * symmop->matrix();
        if (equals(transBmat, myBmat, tol))
            return true;
    }

    return false;
}

VectorPtr
InternalCoordinate::checkCompleteness(
    const Set<ConstInternalCoordinatePtr>& coords, 
    const ConstMoleculePtr& mol
)
{
    if (coords.size() == 0)
    {
        except("There are no coordinates! Surely this must be wrong.");
    }

    RectMatrixPtr Bmat = formBMatrix(coords);
    SymmMatrixPtr BB_T(Bmat.nrow());
    BB_T.accumulate_symmetric_product(Bmat);
    RectMatrixPtr evecs;
    VectorPtr evals;
    BB_T.eigen(evals, evecs);

    int n_nonzero = 0;
    //make this not depend on the order of eigenvalues
    for (int i=0; i < coords.size(); ++i)
    {
        if ( fabs(evals[i]) > 1e-8 ) ++n_nonzero;
    }

    if (n_nonzero != mol->ninternals())
    {
        stringstream sstr;
        sstr << stream_printf("Have %d nonredundant internal coordinates, but need %d!",
                                   n_nonzero, mol->ninternals()) << endl;
        evals.print("BB^T evals", sstr);
        for (int i=0; i < coords.size(); ++i)
        {
            coords[i]->printDetail(sstr);
        }

        string msg = sstr.str();
        except(msg);
    }

    cout << "There are a total of " << n_nonzero << " nonredundant internal coordinates in your set,"
        << " which is correct. Good job!" << endl;

    return evals;
}

VectorPtr
computeChars(
    const ConstPointGroupPtr& pg,
    const ConstVectorPtr& bvec
)
{
    VectorPtr chars = new Vector(pg->nclasses());

    vector<PointGroupClassPtr> classes;
    vector<CoordinateSubspacePtr> subspaces; 

    vector<PointGroupClassPtr>::iterator itclass;
    vector<CoordinateSubspacePtr>::iterator itspace;

    pg->getClasses(classes);
    pg->getSubspaces(subspaces);

    int num = 0;
    for (itclass = classes.begin(); itclass != classes.end(); ++itclass, ++num)
    {
        PointGroupClassPtr cls = *itclass;
        ConstSymmetryOperationPtr symm_op = cls->getRepresentative();
        if (symm_op.get() == NULL)
        {
            except("Null symmetry operation.  This should not be.");
        }

        double character = 0;

        double onecheck = 0;
        for (itspace = subspaces.begin(); itspace != subspaces.end(); ++itspace)
        {
            CoordinateSubspacePtr subspace = *itspace;
            double subspace_char = subspace->character(symm_op);
            VectorPtr proj = subspace->projection(bvec);
            double magnitude = proj.dot(bvec);
            character += subspace_char * magnitude;

            onecheck += magnitude;
        }
        chars.set_element(num, character);
    }

    return chars;
}

void
InternalCoordinate::computeCharacters()
{
    VectorPtr bvec = getBMatrix()->toVector();
    bvec.normalize(); //normalize for use in projections

    chars_ = computeChars(mol_->getPointGroup(), bvec);

    Set<PointGroupPtr> subgroups = mol_->getPointGroup()->getSubgroups();
    for (int i=0; i < subgroups.size(); ++i)
    {
        PointGroupPtr subgroup = subgroups[i];
        VectorPtr subchars = computeChars(subgroup, bvec);
        subgroup_chars_[subgroup->name()] = subchars;
    }
}

double
InternalCoordinate::getConvertedValue() const
{
    return getValue();
}

ConstVectorPtr 
InternalCoordinate::getBVector(const ConstAtomPtr& atom) const
{
    int number = mol_->getNumber(atom);
    return bmatrix_.get_row(number);
}

SymmetryInternalCoordinate::SymmetryInternalCoordinate(
    const ConstVectorPtr& coeffs, 
    const Set<ConstSimpleInternalCoordinatePtr>& coords,
    const ConstMoleculePtr& mol) 
    : InternalCoordinate(mol, "SYMM"), 
    simples_(coords) 
{
    SetRuntime(SymmetryInternalCoordinate);
    coeffs_ = coeffs.copy();
    init();
}

SymmetryInternalCoordinate::SymmetryInternalCoordinate(const XMLArchivePtr& arch)
    : InternalCoordinate(arch)
{
    SetRuntime(SymmetryInternalCoordinate);
    serial_load(simples);
    serial_load(coeffs);
}

void
SymmetryInternalCoordinate::serialize(const XMLArchivePtr& arch) const
{
    InternalCoordinate::serialize(arch);
    serial_save(simples);
    serial_save(coeffs);
}

void
SymmetryInternalCoordinate::intder_normalize()
{
    normalize();
}

void
SymmetryInternalCoordinate::normalize()
{
    coeffs_.normalize();
    recompute();
}

Set<ConstSimpleInternalCoordinatePtr>
SymmetryInternalCoordinate::getSimples() const
{
    return simples_;
}

string
SymmetryInternalCoordinate::connectivityString() const
{
    return ""; //no string
}

double
SymmetryInternalCoordinate::getValueForMolecule(const ConstMoleculePtr& mol) const
{
    double value = 0.0;
    for (int i=0; i < simples_.size(); ++i)
    {
        double coeff = coeffs_[i];
        double simple_value = simples_[i]->getValueForMolecule(mol);
        value += coeff * simple_value;
    }
    return value;
}

int
SymmetryInternalCoordinate::getSimpleCoordinateNumber(
    int simplenum, 
    const Set<ConstSimpleInternalCoordinatePtr>& coords
) const
{
    for (int i=0; i < coords.size(); ++i)
    {
        if (coords[i].get() == simples_[simplenum].get())
            return i;
    }
    return -1;
}

VectorPtr
SymmetryInternalCoordinate::getCoefficients(
    const Set<ConstSimpleInternalCoordinatePtr>& coords
) const
{
	VectorPtr coefs(coords.size());
	for (int i=0; i < coeffs_.n(); ++i)
	{
		int j = getSimpleCoordinateNumber(i, coords);
		coefs.set_element(j, coeffs_.get_element(i));
	}
	return coefs;
}

void
SymmetryInternalCoordinate::printDetail(ostream& os) const
{
    print(os, true, true);
    os << endl;
    for (int i=0; i < simples_.size(); ++i)
    {
        if ( fabs(coeffs_[i]) > 1e-8 )
        {
            os << stream_printf("%8.4f", coeffs_.get_element(i));
            simples_[i]->print(os, false);
            os << endl;
        }
    }
}

double
SymmetryInternalCoordinate::getValue() const
{
    double value = 0.0;
    for (int i=0; i < simples_.size(); ++i)
    {
        double coeff = coeffs_[i];
        ConstSimpleInternalCoordinatePtr simple = simples_[i];
        double simple_value = simple->getValue();
        value += coeff * simple_value;
    }
    return value;
}

void
SymmetryInternalCoordinate::recompute()
{
    VectorPtr bvec;
    for (int n=0; n < mol_->natoms(); ++n)
    {
        ConstAtomPtr atom = mol_->getAtom(n);
        VectorPtr bvec(3);
        for (int i=0; i < simples_.size(); ++i)
        {
            bvec.accumulate(simples_[i]->getBVector(atom) * coeffs_[i]);
        }
        bmatrix_.assign_row(bvec, n);
    }
}

InternalCoordinatePtr
SymmetryInternalCoordinate::copy() const
{
    SymmetryInternalCoordinatePtr newsymm(new SymmetryInternalCoordinate(coeffs_, simples_, mol_));
    return newsymm;
}

InternalCoordinatePtr
SymmetryInternalCoordinate::copy(
    const ConstMoleculePtr& mol,
    const Set<ConstSimpleInternalCoordinatePtr>& simples
) const
{
    vector<ConstSimpleInternalCoordinatePtr> newsimples;

    //figure out which simples I have
    for (int i=0; i < simples_.size(); ++i)
    {
        for (int j=0; j < simples.size(); ++j)
        {
            if (!simples_[i]->matches(simples[j]))
                continue;

            newsimples.push_back(simples[j]);
        }
    }

    if (simples_.size() != newsimples.size())
    {
        printDetail(cerr);
        cerr << endl;
        except(
         stream_printf("Symmetry internal coordinate copy failed.  %d simple internals match, but I need %d", newsimples.size(), simples_.size())
     );
    }

    SymmetryInternalCoordinatePtr newsymm(new SymmetryInternalCoordinate(coeffs_, newsimples, mol));
    return newsymm;
}

SymmetryInternalCoordinatePtr
SymmetryInternalCoordinate::carbon_copy() const
{
    SymmetryInternalCoordinatePtr cpy(new SymmetryInternalCoordinate(coeffs_, simples_, mol_));
    return cpy;
}

SimpleInternalCoordinate::SimpleInternalCoordinate(
    const vector<int>& connect, 
    const ConstMoleculePtr& mol,
    gigstr name
)
    : InternalCoordinate(mol, name), 
      connectivity_(connect)
{
    resetAtoms();
}

SimpleInternalCoordinate::SimpleInternalCoordinate(const XMLArchivePtr& arch)
    : InternalCoordinate(arch)
{
    serial_load(connectivity);
    serial_load(atom_map);
}

void
SimpleInternalCoordinate::serialize(const XMLArchivePtr& arch) const
{
    InternalCoordinate::serialize(arch);

    serial_save(connectivity);
    serial_save(atom_map);
}

void
SimpleInternalCoordinate::intder_normalize()
{
    //do nothing
}

void
SimpleInternalCoordinate::printDetail(ostream& os) const
{
    print(os, true, true);
}

ConstAtomPtr
SimpleInternalCoordinate::getAtom(int idx) const
{
    //connectivity array is 1 based counting, but mol is 0 based
    return mol_->getAtom(connectivity_[idx] - 1);
}

string
SimpleInternalCoordinate::connectivityString() const
{
    std::stringstream sstr;
    for (int i=0; i < connectivity_.size(); ++i)
        sstr << stream_printf("%d ", connectivity_[i]);
    return sstr.str();
}

double
SimpleInternalCoordinate::getValueForMolecule(const ConstMoleculePtr& mol) const
{
    SimpleInternalCoordinatePtr coord = simple_copy(mol);
    double value = coord->getValue();
    return value;
}

void
SimpleInternalCoordinate::dummies(
    vector<int>& dummy_list
) const
{ //return nothing 
}

void
SimpleInternalCoordinate::connectivity(
    vector<int>& connect
) const
{
    connect.clear();
    connect = connectivity_;
}

void
SimpleInternalCoordinate::resetAtoms()
{
    atom_map_.clear();
    for (int num=0; num < connectivity_.size(); ++num)
    {
        int atomnum = connectivity_[num] - 1; //the connect array is 1 based counting... all else is 0
        ConstAtomPtr atom = mol_->getAtom(atomnum);
        atom_map_[num + 1] = atom; //use 1 based counting to simplify code (match Wilson,D,C...)
    }
}

int
SimpleInternalCoordinate::mol_number(int num) const
{
    return connectivity_[num-1] - 1;
}


InternalCoordinatePtr
SimpleInternalCoordinate::copy(
    const ConstMoleculePtr& mol,
    const Set<ConstSimpleInternalCoordinatePtr>& simples
) const
{
    return simple_copy(mol);
}

InternalCoordinatePtr
SimpleInternalCoordinate::copy() const
{
    return simple_copy(mol_);
}

BondLength::BondLength(
    const vector<int>& connect, 
    const ConstMoleculePtr& mol
) : SimpleInternalCoordinate(connect, mol, "STRE")
{
    SetRuntime(BondLength);
    init();
}

BondLength::BondLength(const XMLArchivePtr& arch)
    : SimpleInternalCoordinate(arch)
{
    SetRuntime(BondLength);
}

void
BondLength::serialize(const XMLArchivePtr& arch) const
{
    SimpleInternalCoordinate::serialize(arch);
}

SimpleInternalCoordinatePtr
BondLength::simple_copy(const ConstMoleculePtr& mol) const
{
    BondLengthPtr bond(new BondLength(connectivity_, mol));
    return bond;
}

double
BondLength::compute(const ConstAtomPtr& atom1, const ConstAtomPtr& atom2)
{
    VectorPtr r12vec = atom2->getXYZ() - atom1->getXYZ();
    return r12vec.norm();
}

double
BondLength::getValue() const
{
    ConstAtomPtr at1 = atom_map_.find(1)->second;
    ConstAtomPtr at2 = atom_map_.find(2)->second;
    double length = BondLength::compute(at1, at2);
    return length;
}

void
BondLength::recompute()
{
    ConstAtomPtr at1 = atom_map_[1];
    ConstAtomPtr at2 = atom_map_[2];
    int n1 = mol_number(1);
    int n2 = mol_number(2);

    if (at1.get() == NULL || at2.get() == NULL)
        except("BondLength got assigned null atoms");
    
    VectorPtr b1 = at1->getXYZ() - at2->getXYZ(); b1.normalize();
    VectorPtr b2 = -1.0 * b1;

    bmatrix_.assign_row(b1, n1);
    bmatrix_.assign_row(b2, n2);
}

BondAngle::BondAngle(
    const vector<int>& connect, 
    const ConstMoleculePtr& mol
) : SimpleInternalCoordinate(connect, mol, "BEND")
{
    SetRuntime(BondAngle);
    init();
}

BondAngle::BondAngle(const XMLArchivePtr& arch)
    : SimpleInternalCoordinate(arch)
{
    SetRuntime(BondAngle);
}

void
BondAngle::serialize(const XMLArchivePtr& arch) const
{
    SimpleInternalCoordinate::serialize(arch);
}

SimpleInternalCoordinatePtr
BondAngle::simple_copy(const ConstMoleculePtr& mol) const
{
    BondAnglePtr bend(new BondAngle(connectivity_, mol));
    return bend;
}

bool
BondAngle::isValidValue(double val) const
{
    return (0 <= val) && (val <= PI);
}

double
BondAngle::compute(
    const ConstAtomPtr& atom1, 
    const ConstAtomPtr& atom2, 
    const ConstAtomPtr& atom3
)
{
    //get the position vectors
    ConstVectorPtr at1_xyz = atom1->getXYZ();
    ConstVectorPtr at2_xyz = atom2->getXYZ();
    ConstVectorPtr at3_xyz = atom3->getXYZ();
    
    //get the bond vectors
    VectorPtr r21vec = at1_xyz - at2_xyz;
    VectorPtr r23vec = at3_xyz - at2_xyz;
    //and bond lengths
    double r21 = r21vec.norm();
    double r23 = r23vec.norm();
    double costh = r21vec.dot(r23vec) / r21 / r23;
    return acos(costh);
}

double
BondAngle::getConvertedValue() const
{
    double angle = getValue();
    angle *= RADIAN_TO_DEGREE;
    return angle;
}

double
BondAngle::getValue() const
{
    ConstAtomPtr at1 = atom_map_.find(1)->second;
    ConstAtomPtr at2 = atom_map_.find(2)->second;
    ConstAtomPtr at3 = atom_map_.find(3)->second;
    double angle = BondAngle::compute(at1, at2, at3);
    return angle;
}

void
BondAngle::recompute()
{
    ConstAtomPtr at1 = atom_map_[1];
    ConstAtomPtr at2 = atom_map_[2];
    ConstAtomPtr at3 = atom_map_[3];
    int n1 = mol_number(1);
    int n2 = mol_number(2);
    int n3 = mol_number(3);

    double phi = BondAngle::compute(at1, at2, at3);
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    
    VectorPtr e21 = at1->getXYZ() - at2->getXYZ(); e21.normalize();
    VectorPtr e23 = at3->getXYZ() - at2->getXYZ(); e23.normalize();
    double r21 = BondLength::compute(at2, at1);
    double r23 = BondLength::compute(at2, at3);

    VectorPtr b1 = cosphi * e21 - e23; b1.scale(1.0/r21/sinphi);
    VectorPtr b2 = (r21 - r23*cosphi) * e21 + (r23 - r21*cosphi) * e23; b2.scale(1.0/r21/r23/sinphi);
    VectorPtr b3 = cosphi * e23 - e21; b3.scale(1.0/r23/sinphi);

    bmatrix_.assign_row(b1, n1);
    bmatrix_.assign_row(b2, n2);
    bmatrix_.assign_row(b3, n3);
}

OutOfPlaneBend::OutOfPlaneBend(
    const vector<int>& connect, 
    const ConstMoleculePtr& mol
) : SimpleInternalCoordinate(connect, mol, "OUT")
{
    SetRuntime(OutOfPlaneBend);
    init();
}

OutOfPlaneBend::OutOfPlaneBend(const XMLArchivePtr& arch)
    : SimpleInternalCoordinate(arch)
{
    SetRuntime(OutOfPlaneBend);
}

void
OutOfPlaneBend::serialize(const XMLArchivePtr& arch) const
{
    SimpleInternalCoordinate::serialize(arch);
}

SimpleInternalCoordinatePtr
OutOfPlaneBend::simple_copy(const ConstMoleculePtr& mol) const
{
    OutOfPlaneBendPtr bend(new OutOfPlaneBend(connectivity_, mol));
    return bend;
}

bool
OutOfPlaneBend::isValidValue(double val) const
{
    if (val > -PI/2 or val < PI/2)
        return true;
    else
        return false;
}

double
OutOfPlaneBend::getValue() const
{
    ConstAtomPtr at1 = atom_map_.find(1)->second;
    ConstAtomPtr at2 = atom_map_.find(2)->second;
    ConstAtomPtr at3 = atom_map_.find(3)->second;
    ConstAtomPtr at4 = atom_map_.find(4)->second;

    double phi1 = BondAngle::compute(at3, at2, at4);
    VectorPtr e21 = at1->getXYZ() - at2->getXYZ(); e21.normalize();
    VectorPtr e23 = at3->getXYZ() - at2->getXYZ(); e23.normalize();
    VectorPtr e24 = at4->getXYZ() - at2->getXYZ(); e24.normalize();

    double angle = cross(e23, e24).dot(e21) / sin(phi1);

    return asin(angle);
}

void
OutOfPlaneBend::recompute()
{
    ConstAtomPtr at1 = atom_map_[1];
    ConstAtomPtr at2 = atom_map_[2];
    ConstAtomPtr at3 = atom_map_[3];
    ConstAtomPtr at4 = atom_map_[4];
    int n1 = mol_number(1);
    int n2 = mol_number(2);
    int n3 = mol_number(3);
    int n4 = mol_number(4);

    double theta = getValue();
    double costheta = cos(theta);
    double tantheta = tan(theta);
    double phi1 = BondAngle::compute(at3, at2, at4);
    double sinphi1 = sin(phi1);
    double cosphi1 = cos(phi1);
    double sin2phi1 = sinphi1 * sinphi1;

    double r21 = BondLength::compute(at1, at2);
    double r23 = BondLength::compute(at3, at2);
    double r24 = BondLength::compute(at4, at2);

    VectorPtr e21 = at1->getXYZ() - at2->getXYZ(); e21.normalize();
    VectorPtr e23 = at3->getXYZ() - at2->getXYZ(); e23.normalize();
    VectorPtr e24 = at4->getXYZ() - at2->getXYZ(); e24.normalize();

    VectorPtr b1 = cross(e23, e24) * (1.0/r21/costheta/sinphi1) - e21 * (tantheta/r21);
    VectorPtr b3 = cross(e24, e21) * (1.0/r23/costheta/sinphi1) - (e23 - cosphi1 * e24) * (tantheta/sin2phi1/r23);
    VectorPtr b4 = cross(e21, e23) * (1.0/r24/costheta/sinphi1) - (e24 - cosphi1 * e23) * (tantheta/sin2phi1/r24);
    VectorPtr b2 = -1.0 * b1 - b3 - b4;

    bmatrix_.assign_row(b1, n1);
    bmatrix_.assign_row(b2, n2);
    bmatrix_.assign_row(b3, n3);
    bmatrix_.assign_row(b4, n4);
}

PeriodicCoordinate::PeriodicCoordinate(
    const vector<int>& connect, 
    const ConstMoleculePtr& mol, 
    double discontinuity,
    double period,
    gigstr name
) : SimpleInternalCoordinate(connect, mol, name)
{
    period_ = period;
    discontinuity_ = discontinuity;
}

PeriodicCoordinate::PeriodicCoordinate(const XMLArchivePtr& arch)
    : SimpleInternalCoordinate(arch)
{
    serial_load(discontinuity);
}

void
PeriodicCoordinate::serialize(const XMLArchivePtr& arch) const
{
    SimpleInternalCoordinate::serialize(arch);
    serial_save(discontinuity);
}

double
PeriodicCoordinate::canonicalizeValue(double val) const
{
    double newval = val - discontinuity_;
    double nperiods = floor(newval / period_);
    double remainder = newval - nperiods * period_; //find the modulo part
    double finalval = remainder + discontinuity_;
    return finalval;

#if 0 
    if      ( val < discontinuity_ )            return 2*PI + val;
    else if ( val > (2*PI + discontinuity_))    return val - 2*PI;
    else                                        return val;
#endif
}

Torsion::Torsion(
    const vector<int>& connect, 
    const ConstMoleculePtr& mol, 
    double discontinuity
) : PeriodicCoordinate(connect, mol, discontinuity, 2*PI, "TORS")
{
    SetRuntime(Torsion);
    init();
}

Torsion::Torsion(const XMLArchivePtr& arch)
    : PeriodicCoordinate(arch)
{
    SetRuntime(Torsion);
}

void
Torsion::serialize(const XMLArchivePtr& arch) const
{
    PeriodicCoordinate::serialize(arch);
}

SimpleInternalCoordinatePtr
Torsion::simple_copy(const ConstMoleculePtr& mol) const
{
    TorsionPtr coord(new Torsion(connectivity_, mol, discontinuity_));
    return coord;
}

double
Torsion::getValue() const
{
    ConstAtomPtr at1 = atom_map_.find(1)->second;
    ConstAtomPtr at2 = atom_map_.find(2)->second;
    ConstAtomPtr at3 = atom_map_.find(3)->second;
    ConstAtomPtr at4 = atom_map_.find(4)->second;

    //get some unit vectors
    VectorPtr e12 = at2->getXYZ() - at1->getXYZ(); e12.normalize();
    VectorPtr e21 = -1.0 * e12;
    VectorPtr e23 = at3->getXYZ() - at2->getXYZ(); e23.normalize();
    VectorPtr e34 = at4->getXYZ() - at3->getXYZ(); e34.normalize();

    //get some angles
    double phi2 = BondAngle::compute(at1, at2, at3);
    double phi3 = BondAngle::compute(at2, at3, at4);
    double sinphi2 = sin(phi2); double sinphi3 = sin(phi3);

    //and compute
    double costau = cross(e12, e23).dot(cross(e23, e34)) / (sinphi2 * sinphi3);
    double tau = inverse_cosine(costau);

    //but we now have to determine the sign... from the cross products of e21 and e34
    //we can determine if we have a clockwise or counter clockwise rotation
    VectorPtr sign_cross = cross(e21, e34);
    double sign_val = sign_cross.dot(e23);
    if (sign_val < 0) tau *= -1.0;

    return canonicalizeValue(tau);
}

void
Torsion::recompute()
{
    ConstAtomPtr at1 = atom_map_[1];
    ConstAtomPtr at2 = atom_map_[2];
    ConstAtomPtr at3 = atom_map_[3];
    ConstAtomPtr at4 = atom_map_[4];
    int n1 = mol_number(1);
    int n2 = mol_number(2);
    int n3 = mol_number(3);
    int n4 = mol_number(4);

    //all of the intermediates that we will need
    double r12, r23, r32, r43, phi2, phi3, sinphi2, sinphi3, cosphi2, cosphi3, sin2phi2, sin2phi3;
    VectorPtr e12, e23, e32, e43;

    r12 = BondLength::compute(at1, at2);
    r23 = r32 = BondLength::compute(at2, at3);
    r43 = BondLength::compute(at3, at4);
    phi2 = BondAngle::compute(at1, at2, at3);
    phi3 = BondAngle::compute(at2, at3, at4);
    sinphi2 = sin(phi2);
    sinphi3 = sin(phi3);
    cosphi2 = cos(phi2);
    cosphi3 = cos(phi3);
    sin2phi2 = sinphi2*sinphi2;
    sin2phi3 = sinphi3*sinphi3;

    e12 = at2->getXYZ() - at1->getXYZ(); e12.normalize();
    e23 = at3->getXYZ() - at2->getXYZ(); e23.normalize();
    e32 = -1.0 * e23;
    e43 = at3->getXYZ() - at4->getXYZ(); e43.normalize();

    VectorPtr b1 = cross(e12, e23) * (-1.0 / (r12 * sin2phi2));
    VectorPtr b2 = ((r23 - r12 * cosphi2)/(r23*r12*sin2phi2)) * cross(e12, e23)
                     + ((cosphi3)/(r23*sin2phi3)) * cross(e43, e32);
    VectorPtr b3 = ((r32 - r43 * cosphi3)/(r32*r43*sin2phi3)) * cross(e43, e32)
                     + ((cosphi2)/(r32*sin2phi2)) * cross(e12, e23);
    VectorPtr b4 = cross(e43, e32) * (-1.0 / (r43 * sin2phi3));

    bmatrix_.assign_row(b1, n1);
    bmatrix_.assign_row(b2, n2);
    bmatrix_.assign_row(b3, n3);
    bmatrix_.assign_row(b4, n4);
}


double
Torsion::getConvertedValue() const
{
    double angle = getValue();
    angle *= RADIAN_TO_DEGREE;
    return angle;
}

LinX::LinX(
    const vector<int>& connect, 
    const ConstMoleculePtr& mol
) : SimpleInternalCoordinate(connect, mol, "LINX")
{
    SetRuntime(LinX);
    vector<int> angle_atoms;
    angle_atoms.push_back(connect[1]);
    angle_atoms.push_back(connect[2]);
    angle_atoms.push_back(connect[3]);

    angle_ = new BondAngle(angle_atoms, mol);
    torsion_ = new Torsion(connect, mol);
    init();
}

LinX::LinX(const XMLArchivePtr& arch)
    : SimpleInternalCoordinate(arch)
{
    SetRuntime(LinX);
}

void
LinX::serialize(const XMLArchivePtr& arch) const
{
    SimpleInternalCoordinate::serialize(arch);
}

SimpleInternalCoordinatePtr
LinX::simple_copy(const ConstMoleculePtr& mol) const
{
    LinXPtr coord(new LinX(connectivity_, mol));
    return coord;
}

double
LinX::getValue() const
{
    ConstAtomPtr at1 = atom_map_.find(1)->second;
    ConstAtomPtr at2 = atom_map_.find(2)->second;
    ConstAtomPtr at3 = atom_map_.find(3)->second;
    ConstAtomPtr at4 = atom_map_.find(4)->second;

    VectorPtr e21 = at1->getXYZ() - at2->getXYZ(); e21.normalize();
    VectorPtr e32 = at2->getXYZ() - at3->getXYZ(); e32.normalize();
    VectorPtr e23 = -1.0 * e32;
    VectorPtr e34 = at4->getXYZ() - at3->getXYZ(); e34.normalize();

    double sintheta = sin(BondAngle::compute(at1, at2, at3));
    double val = (1.0 / sintheta) * cross(e21, e23).dot(cross(e32, e34)) ;
    return val;
}

void
LinX::recompute()
{
    ConstAtomPtr at1 = atom_map_[1];
    ConstAtomPtr at2 = atom_map_[2];
    ConstAtomPtr at3 = atom_map_[3];
    ConstAtomPtr at4 = atom_map_[4];
    int n1 = mol_number(1);
    int n2 = mol_number(2);
    int n3 = mol_number(3);
    int n4 = mol_number(4);

    angle_->recompute();
    torsion_->recompute();

    double tau = torsion_->getValue();
    double theta = angle_->getValue();
    double sintau = sin(tau);
    double costau = cos(tau);
    double costheta = cos(theta);
    double sintheta = sin(theta);

    VectorPtr b1 = -sintau * sintheta * torsion_->getBVector(at1);
    VectorPtr b2 = -sintau * sintheta * torsion_->getBVector(at2) + costheta * costau * angle_->getBVector(at2);
    VectorPtr b3 = -sintau * sintheta * torsion_->getBVector(at3) + costheta * costau * angle_->getBVector(at3);
    VectorPtr b4 = -sintau * sintheta * torsion_->getBVector(at4) + costheta * costau * angle_->getBVector(at4);

    bmatrix_.assign_row(b1, n1);
    bmatrix_.assign_row(b2, n2);
    bmatrix_.assign_row(b3, n3);
    bmatrix_.assign_row(b4, n4);
}

LinY::LinY(
    const vector<int>& connect, 
    const ConstMoleculePtr& mol
) 
    : SimpleInternalCoordinate(connect, mol, "LINY")
{
    SetRuntime(LinY);
    vector<int> angle_atoms;
    angle_atoms.push_back(connect[1]);
    angle_atoms.push_back(connect[2]);
    angle_atoms.push_back(connect[3]);

    angle_ = new BondAngle(angle_atoms, mol);
    torsion_ = new Torsion(connect, mol);
    init();
}

LinY::LinY(const XMLArchivePtr& arch)
    : SimpleInternalCoordinate(arch)
{
    SetRuntime(LinY);
}

void
LinY::serialize(const XMLArchivePtr& arch) const
{
    SimpleInternalCoordinate::serialize(arch);
}

SimpleInternalCoordinatePtr
LinY::simple_copy(const ConstMoleculePtr& mol) const
{
    LinYPtr coord(new LinY(connectivity_, mol));
    return coord;
}

double
LinY::getValue() const
{
    ConstAtomPtr at1 = atom_map_.find(1)->second;
    ConstAtomPtr at2 = atom_map_.find(2)->second;
    ConstAtomPtr at3 = atom_map_.find(3)->second;
    ConstAtomPtr at4 = atom_map_.find(4)->second;

    VectorPtr e21 = at1->getXYZ() - at2->getXYZ(); e21.normalize();
    VectorPtr e32 = at2->getXYZ() - at3->getXYZ(); e32.normalize();
    VectorPtr e34 = at4->getXYZ() - at3->getXYZ(); e34.normalize();

    double sintheta = sin(BondAngle::compute(at1, at2, at3));
    double val = (1.0 / sintheta) * e21.dot(cross(e32, e34));
    return val;
}

void
LinY::recompute()
{
    ConstAtomPtr at1 = atom_map_[1];
    ConstAtomPtr at2 = atom_map_[2];
    ConstAtomPtr at3 = atom_map_[3];
    ConstAtomPtr at4 = atom_map_[4];
    int n1 = mol_number(1);
    int n2 = mol_number(2);
    int n3 = mol_number(3);
    int n4 = mol_number(4);

    angle_->recompute();
    torsion_->recompute();

    double tau = torsion_->getValue();
    double theta = angle_->getValue();
    double sintau = sin(tau);
    double costau = cos(tau);
    double costheta = cos(theta);
    double sintheta = sin(theta);

    VectorPtr b1 = costau * sintheta * torsion_->getBVector(at1);
    VectorPtr b2 = costau * sintheta * torsion_->getBVector(at2) + costheta * sintau * angle_->getBVector(at2);
    VectorPtr b3 = costau * sintheta * torsion_->getBVector(at3) + costheta * sintau * angle_->getBVector(at3);
    VectorPtr b4 = costau * sintheta * torsion_->getBVector(at4) + costheta * sintau * angle_->getBVector(at4);

    bmatrix_.assign_row(b1, n1);
    bmatrix_.assign_row(b2, n2);
    bmatrix_.assign_row(b3, n3);
    bmatrix_.assign_row(b4, n4);
}

Lin1::Lin1(
    const vector<int>& connect, 
    const ConstMoleculePtr& mol, 
    const vector<double>& vec
) 
    : SimpleInternalCoordinate(connect, mol, "LIN1"),
      ed_(XYZ_DIM)
{
    SetRuntime(Lin1);
    for (int i=0; i < XYZ_DIM; ++i) 
        ed_.set_element(i, vec[i]);
    init();

    //const cast the molecule for now... no other way to do this
    //at least we are very conscious here of ditching the const modifier
    MoleculePtr unconstmol = boost::const_pointer_cast<Molecule, const Molecule>(mol);
    dummy_ = unconstmol->registerDummyAtom(ed_);
}

Lin1::Lin1(const XMLArchivePtr& arch)
    : SimpleInternalCoordinate(arch)
{
    SetRuntime(Lin1);
    serial_load(ed);
}

void
Lin1::serialize(const XMLArchivePtr& arch) const
{
    SimpleInternalCoordinate::serialize(arch);
    serial_save(ed);
}

bool
Lin1::matches(const ConstSimpleInternalCoordinatePtr& coord) const
{
    bool check = SimpleInternalCoordinate::matches(coord);
    if (!check)
        return false;

    ConstLin1Ptr cast(boost::dynamic_pointer_cast<const Lin1,const SimpleInternalCoordinate>(coord));
    if (!cast)
        return false;

    VectorPtr diff = ed_ - cast->ed();
    return (diff.maxabs() < 1e-8);
}

SimpleInternalCoordinatePtr
Lin1::simple_copy(const ConstMoleculePtr& mol) const
{
    vector<double> v;
    for (int i=0; i < 3; ++i) 
        v.push_back(ed_[i]);
    Lin1Ptr coord(new Lin1(connectivity_, mol, v));
    return coord;
}

void
Lin1::dummies(
    vector<int>& dummy_list
) const
{ 
    dummy_list.push_back(dummy_);
}

ConstVectorPtr
Lin1::ed() const
{
    return ed_;
}

double
Lin1::getValue() const
{
    ConstAtomPtr at1 = atom_map_.find(1)->second;
    ConstAtomPtr at2 = atom_map_.find(2)->second;
    ConstAtomPtr at3 = atom_map_.find(3)->second;

    VectorPtr e21 = at1->getXYZ() - at2->getXYZ(); e21.normalize();
    VectorPtr e23 = at3->getXYZ() - at2->getXYZ(); e23.normalize();
    double sinval = ed_.dot(cross(e23, e21));
    double val = inverse_sine(sinval);
    return val;
}

void
Lin1::recompute()
{
    ConstAtomPtr at1 = atom_map_[1];
    ConstAtomPtr at2 = atom_map_[2];
    ConstAtomPtr at3 = atom_map_[3];
    int n1 = mol_number(1);
    int n2 = mol_number(2);
    int n3 = mol_number(3);

    VectorPtr r21v = at1->getXYZ() - at2->getXYZ(); 
    VectorPtr e21 = r21v.copy(); e21.normalize();
    VectorPtr r23v = at3->getXYZ() - at2->getXYZ(); 
    VectorPtr e23 = r23v.copy(); e23.normalize();
    double r21 = BondLength::compute(at1, at2);
    double oor21 = 1.0/r21;
    double r23 = BondLength::compute(at2, at3);
    double oor23 = 1.0/r23;

    double costheta = cos(getValue());

    //build the 
    VectorPtr ex = new Vector(3); ex.set_element(X, 1.0);
    VectorPtr ey = new Vector(3); ey.set_element(Y, 1.0);
    VectorPtr ez = new Vector(3); ez.set_element(Z, 1.0);

    VectorPtr ex21 = cross(ex, e21);
    VectorPtr ey21 = cross(ey, e21);
    VectorPtr ez21 = cross(ez, e21);
    VectorPtr ex23 = cross(ex, e23);
    VectorPtr ey23 = cross(ey, e23);
    VectorPtr ez23 = cross(ez, e23);
    VectorPtr e23_21 = cross(e23, e21);

    //do atom1
    VectorPtr d1x = ex23 + (oor21) * r21v[X] * e23_21;
    double bv1x = ed_.dot(d1x);
    VectorPtr d1y = ey23 + (oor21) * r21v[Y] * e23_21; 
    double bv1y = ed_.dot(d1y);
    VectorPtr d1z = ez23 + (oor21) * r21v[Z] * e23_21; 
    double bv1z = ed_.dot(d1z);
    VectorPtr b1 = new Vector(3);
    b1.set_element(X, bv1x);
    b1.set_element(Y, bv1y);
    b1.set_element(Z, bv1z);
    b1.scale(-1.0/costheta/r21);

    //do atom2
    double scalefac23 = (-1.0/costheta/r23);
    double scalefac12 = (1.0/costheta/r21);
    VectorPtr bv2x = scalefac23 * (ex21 - e23_21 * r23v[X] * oor23)
                     + scalefac12 * (ex23 + e23_21 * r21v[X] * oor21);
    VectorPtr bv2y = scalefac23 * (ey21 - e23_21 * r23v[Y] * oor23)
                     + scalefac12 * (ey23 + e23_21 * r21v[Y] * oor21);
    VectorPtr bv2z = scalefac23 * (ez21 - e23_21 * r23v[Z] * oor23)
                     + scalefac12 * (ez23 + e23_21 * r21v[Z] * oor21);
    VectorPtr b2 = new Vector(3);
    b2.set_element(X, ed_.dot(bv2x));
    b2.set_element(Y, ed_.dot(bv2y));
    b2.set_element(Z, ed_.dot(bv2z));

    //do atom3
    VectorPtr d3x = ex21 - (oor23) * r23v[X] * e23_21;
    double bv3x = ed_.dot(d3x);
    VectorPtr d3y = ey21 - (oor23) * r23v[Y] * e23_21; 
    double bv3y = ed_.dot(d3y);
    VectorPtr d3z = ez21 - (oor23) * r23v[Z] * e23_21; 
    double bv3z = ed_.dot(d3z);
    VectorPtr b3 = new Vector(3);
    b3.set_element(X, bv3x);
    b3.set_element(Y, bv3y);
    b3.set_element(Z, bv3z);
    b3.scale(1.0/costheta/r23);

    bmatrix_.assign_row(b1, n1);
    bmatrix_.assign_row(b2, n2);
    bmatrix_.assign_row(b3, n3);
}

ConstMoleculePtr
InternalCoordinate::mol() const
{
    return mol_;
}

int
InternalCoordinate::getCoordinateNumber(
    const ConstInternalCoordinatePtr& coord,
    const Set<ConstInternalCoordinatePtr >& coords
)
{
    for (int i=0; i < coords.size(); ++i)
    {
        if (coord == coords[i]) return i;
    }

    stringstream sstr;
    sstr << "Cannot find match for coordinate: "; coord->print(sstr);
    except(sstr.str());
}

RectMatrixPtr
InternalCoordinate::mapBMatrix(
    const ConstSymmetryOperationPtr& op
) const
{
    return op->getPermutationMatrix(mol_) * getBMatrix() * op->matrix();
}

RectMatrixPtr
InternalCoordinate::formBMatrix(
    const Set<ConstInternalCoordinatePtr>& coords
)
{
    if (coords.size() == 0)
        except("No coordinates available to form B matrix");

    ConstMoleculePtr mol = coords[0]->mol();
    int natoms = mol->natoms();
    RectMatrixPtr B(coords.size(), XYZ_DIM * natoms); 
    for (int coord_num=0; coord_num < coords.size(); coord_num++)
    {
        ConstInternalCoordinatePtr coord = coords[coord_num];
        for (int atom=0; atom < natoms; atom++)
        {
            ConstAtomPtr current_atom = mol->getAtom(atom); 
            ConstVectorPtr Bvec = coord->getBVector(current_atom);
            for (int x=0; x < XYZ_DIM; x++)
            {
                int row = coord_num; int col = XYZ_DIM*atom + x;
                B.set_element(row, col, Bvec.get_element(x));
            }
        }
    }
    return B;
}

void
InternalCoordinate::testBVectors(double tol)
{
    double dispsize = KeywordSet::getKeyword("bvecdisp")->getValueDouble();

    int natoms = mol_->natoms();
    print();
    for (int n=0; n < natoms; n++)
    {
        //get the numeric b vector
        double xderiv = getUnitVectorDerivative(n, X, dispsize);
        double yderiv = getUnitVectorDerivative(n, Y, dispsize);
        double zderiv = getUnitVectorDerivative(n, Z, dispsize);

        //get the analytic b vector
        ConstAtomPtr atom = mol_->getAtom(n);
        ConstVectorPtr analbvec = getBVector(atom);

        //accumulate the difference
        double diff = 0;
        diff += fabs(analbvec.get_element(X) - xderiv);
        diff += fabs(analbvec.get_element(Y) - yderiv);
        diff += fabs(analbvec.get_element(Z) - zderiv);

        if (diff > tol)
        {
            cerr << "B Vector Test failed" << endl;
            print(cerr);

            cerr << "Numerical B vector" << endl;
            cerr << stream_printf("<%12.8f, %12.8f, %12.8f", xderiv, yderiv, zderiv) << endl;

            analbvec.print("Analytic B Vector", cerr);
            abort();
        }
    }
}

double
InternalCoordinate::getUnitVectorDerivative(int n, int coord, double disp)
{
    VectorPtr vec3 = new Vector(3);
    VectorPtr xdisp = vec3.copy(); xdisp.set_element(coord, disp);
    
    MoleculePtr testmol(new Molecule(mol_));
    double deriv;
    //get the minus 2 displacement
    testmol -> displaceAtom(n, -2 * xdisp);
    double fm2 = getValue();
    //get the minus 1 displacement
    testmol -> displaceAtom(n, xdisp);
    double fm1 = getValue();
    //get the plus 1 displacement
    testmol -> displaceAtom(n, 2 * xdisp);
    double fp1 = getValue();
    //get the plus 2 displacement
    testmol -> displaceAtom(n, xdisp);
    double fp2 = getValue();
    //recenter molcule
    testmol -> displaceAtom(n, -2 * xdisp);
    deriv = (fm2 - 8 * fm1 + 8 * fp1 - fp2) / 12 / disp;
    return deriv;
}

VectorPtr
InternalCoordinate::addDelocalizedInternalCoordinates(
    const ConstMoleculePtr& mol,
    const Set<ConstSimpleInternalCoordinatePtr>& simples,
    vector<InternalCoordinatePtr>& coords,
    double tol
)
{
    vector<SymmetryInternalCoordinatePtr> symmcoords;
    VectorPtr evals = InternalCoordinate::addDelocalizedInternalCoordinates(mol, simples, symmcoords, tol);
    for (int i=0; i < symmcoords.size(); ++i)
        coords.push_back(symmcoords[i]);
    return evals;
}

VectorPtr
InternalCoordinate::addDelocalizedInternalCoordinates(
    const ConstMoleculePtr& mol,
    const Set<ConstSimpleInternalCoordinatePtr>& simples,
    vector<SymmetryInternalCoordinatePtr>& coords,
    double tol
)
{
    RectMatrixPtr Bmat = formBMatrix(simples);
    SymmMatrixPtr BB_T(Bmat.nrow());
    BB_T.accumulate_symmetric_product(Bmat);
    RectMatrixPtr evecs;
    VectorPtr evals;
    BB_T.eigen(evals, evecs);

	int debug = KeywordSet::getKeyword("coordinate debug")->getValueInteger();

    //figure out which evals are nonzero. this routine does not depend on
    //the order of eigenvalues
    vector<int> nonzero;
    for (int i=0; i < simples.size(); ++i)
    {
        double eval = fabs(evals[i]);
        if ( eval > tol )
            nonzero.push_back(i);
    }
    int nzero = simples.size() - nonzero.size();
	
	if (debug)
    	cout << "Throwing out " << nzero << " coordinates due to small eigenvalues" << endl;

    if (nonzero.size() == 0)
        except("Unable to build delocalized internal coordinates. No nonzero eigenvalues");
        

    VectorPtr nonzero_evals = new Vector(nonzero.size());

    vector<int>::iterator it;
    int iblock = 0;

    RectMatrixPtr nonzero_evecs(simples.size(), nonzero.size());
    for (it = nonzero.begin(); it != nonzero.end(); ++it, ++iblock)
    {
        int itot = *it;
        VectorPtr coeffs = evecs.get_column(itot);
        SymmetryInternalCoordinatePtr new_coord = new SymmetryInternalCoordinate(coeffs, simples, mol);
        coords.push_back(new_coord);
        nonzero_evals.set_element(iblock, evals.get_element(itot));
        nonzero_evecs.assign_column(coeffs, iblock);
    }

    //now orthogonalize the basis
    Bmat = formBMatrix(coords);
    BB_T = SymmMatrixPtr(nonzero.size());
    BB_T.accumulate_symmetric_product(Bmat);
    SymmMatrixPtr orthog = BB_T.invsqrt_matrix();
    RectMatrixPtr final_evecs = nonzero_evecs * orthog;

    coords.clear();
    for (int i=0; i < final_evecs.ncol(); ++i)
    {
        VectorPtr coeffs = final_evecs.get_column(i);
        SymmetryInternalCoordinatePtr new_coord = new SymmetryInternalCoordinate(coeffs, simples, mol);
        coords.push_back(new_coord);
    }

    return nonzero_evals;
}

CoordinateSubspace::CoordinateSubspace(
    const ConstMoleculePtr& mol
) : mol_(mol), madeabelian_(false)
{
    SetRuntime(CoordinateSubspace);
}

CoordinateSubspace::CoordinateSubspace(
    const XMLArchivePtr& arch
) : Serializable(arch)
{
    SetRuntime(CoordinateSubspace);
    serial_load(madeabelian);
    serial_load(coords);
    serial_load(mol);
}

void
CoordinateSubspace::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(madeabelian);
    serial_save(coords);
    serial_save(mol);
}

int
CoordinateSubspace::order() const
{
    return coords_.size();
}

void
CoordinateSubspace::addCoordinate(
    const SymmetryInternalCoordinatePtr& coord
)
{
    coords_.push_back(coord);
}

double
CoordinateSubspace::character(
    const ConstSymmetryOperationPtr& oper
) const
{
    if (coords_.size() == 0)
        except("CoordinateSubspace cannot compute character for subspace with not coordinates");

    RectMatrixPtr permutation = oper->getPermutationMatrix(mol_);

    double trace = 0;
    //we are essentially computing the diagonal elements of the tform matrix
    for (int i=0; i < coords_.size(); ++i)
    {
        ConstInternalCoordinatePtr coord = coords_[i];
        //transform the b vector
        RectMatrixPtr trans_bmat = permutation * coord->getBMatrix() * oper->matrix();
        VectorPtr trans_bvec = trans_bmat.toVector();
        trans_bvec.normalize();

        VectorPtr ref_bvec = coord->getBMatrix()->toVector();
        ref_bvec.normalize();
        
        //determine the projection onto the original basis vector
        double overlap = ref_bvec.dot(trans_bvec);
        
        trace += overlap;
    }

    return trace;
}

ConstSymmetryInternalCoordinatePtr
CoordinateSubspace::getCoordinate(int i) const
{
    if (i >= coords_.size())
        except(stream_printf("Invalid coordinate number %d for subspace of order %d", i, coords_.size()));

    return coords_[i];
}

SymmetryInternalCoordinatePtr
CoordinateSubspace::getCoordinate(int i)
{
    if (i >= coords_.size())
        except(stream_printf("Invalid coordinate number %d for subspace of order %d", i, coords_.size()));

    return coords_[i];
}

void
CoordinateSubspace::printCharacters(
    const Set<ConstSymmetryOperationPtr>& opers,
    ostream& os
) const
{
    vector<double> chars;
    Set<ConstSymmetryOperationPtr>::const_iterator it;
    for (it = opers.begin(); it != opers.end(); ++it)
        chars.push_back( character(*it) );

    os << stream_printf("subspace dim = %d ", coords_.size()); coords_[0]->printDetail();
    for (int col=0; col < opers.size(); ++col)
        cout << stream_printf("%8s", opers[col]->name().c_str());
    os << endl;
    for (int col=0; col < opers.size(); ++col)
        os << stream_printf("%8.3f", chars[col]);
}

VectorPtr
CoordinateSubspace::projection(
    const ConstVectorPtr& bvec
)
{
    VectorPtr proj = bvec.clone();
    for (int i=0; i < coords_.size(); ++i)
    {
        ConstInternalCoordinatePtr coord = coords_[i];
        VectorPtr ref_bvec = coord->getBMatrix()->toVector();

        double dot = ref_bvec.dot(bvec);
        proj.accumulate(dot * ref_bvec);
    }
    return proj;
}

void
CoordinateSubspace::print(ostream& os) const
{
    cout << stream_printf("Subspace of order %d", order()) << endl;
    for (int i=0; i < coords_.size(); ++i)
        coords_[i]->printDetail(os);
}

void
CoordinateSubspace::abelianify(const ConstPointGroupPtr& pg)
{
    //nothing to do
    if (coords_.size() == 1 || madeabelian_)
        return;

    vector<CoordinateSubspacePtr> subspaces;
    pg->formBasis(coords_, subspaces);

    //see if the number of subspaces generated is correct
    if (subspaces.size() != coords_.size())
    {
        madeabelian_ = false;
        return;
    }

    //reset the coordinates
    coords_.clear();
    for (int i=0; i < subspaces.size(); ++i)
    {
        coords_.push_back(subspaces[i]->getCoordinate(0));
    }

    //compute the characters - we might be able to do better
    vector<PointGroupClassPtr> classes;
    pg->getClasses(classes);
    vector<VectorPtr> charvec;
    for (int i=0; i < subspaces.size(); ++i)
    {
        VectorPtr v(classes.size());
        for (int cls=0; cls < classes.size(); ++cls)
        {
            ConstSymmetryOperationPtr op = classes[cls]->getRepresentative();
            double subchar = subspaces[i]->character(op);
            v.set_element(cls, subchar);
        }
        charvec.push_back(v);
    }

    //test the dot products to see if we can do better
    for (int i=0; i < charvec.size(); ++i)
    {
        for (int j=i+1; j < charvec.size(); ++j)
        {
            double dotprod = charvec[i].dot(charvec[j]);
            if ( fabs(dotprod) > 1e-6 )
            {
                //accept what we currently have, but try again maybe
                madeabelian_ = false;
                return;
            }
        }
    }

    //everything looks good
    madeabelian_ = true;
}

