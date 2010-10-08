#include <src/utilities.h>
#include <src/hessian.h>
#include <src/qff.h>
#include <src/molecule.h>
#include <src/symmetry.h>
#include <src/defines.h>
#include <src/exception.h>
#include <src/gigmatrix.h>
#include <sstream>

#define X 0
#define Y 1
#define Z 2


using namespace gigide;
using namespace std;
using namespace smartptr;

SerialDeclare(Molecule)
SerialDeclare(Atom)

/** Set static values to being */
bool Molecule::statics_done = false;
bool Molecule::needpg_ = false;
map<string, double> Molecule::MASSES;

//static smartptr::SerialRuntime Molecule_serialize("Molecule", &smartptr::serialize<Molecule>);
//static smartptr::SerialRuntime Atom_serialize("Atom", &smartptr::serialize<Atom>);

Molecule::Molecule(
    ConstRectMatrixPtr xyz, 
    const vector<string>& atomlist
)
{
    SetRuntime(Molecule);
    initialize_statics();

    if (atomlist.size() == 0)
        except("No atoms.  Cannot create molecule");

    xyz_ = xyz.copy();
    for (int i=0; i < atomlist.size(); ++i)
    {
        string symbol = atomlist[i];
        VectorPtr atom_xyz = xyz.get_row(i);
        AtomPtr atom(new Atom(atom_xyz, symbol));
        atoms_.push_back(atom);
    }

    init();

    incref();
    pg_ = new PointGroup(this);
    decref();
}

Molecule::Molecule(const ConstMoleculePtr& mol)
{
    SetRuntime(Molecule);
    initialize_statics();

    vector<string> atomlist;
    xyz_ = mol->getXYZ().copy();
    for (int i=0; i < mol->natoms(); ++i)
    {
        ConstAtomPtr atom = mol->getAtom(i);
        AtomPtr copy(new Atom(atom->getXYZ(), atom->symbol()));
        atoms_.push_back(copy);
    }

    init();

    incref();
    pg_ = new PointGroup(this);
    decref();
}

Molecule::Molecule(const XMLArchivePtr& arch)
    : Serializable(arch)
{
    SetRuntime(Molecule);
    initialize_statics();

    serial_load(xyz);
    serial_load(pg);
    serial_load(atoms);
    serial_load(dummies);

    init();
}

Molecule::~Molecule()
{
}

void
Molecule::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(xyz);
    serial_save(pg);
    serial_save(atoms);
    serial_save(dummies);
}

void
Molecule::init()
{
    for (int i=0; i < natoms(); ++i)
    {
        //get a map linking atom to number
        number_map_[atoms_[i].get()] = i;
    }

    //figure out whether it is linear
    if (atoms_.size() == 1) 
    {
        islinear_ = false;
    }
    else if (atoms_.size() == 2)
    {
        islinear_ = true;
    }
    else
    {
        VectorPtr e01 = atoms_[1]->getXYZ() - atoms_[0]->getXYZ(); e01.normalize();
        for (int i=2; i < atoms_.size(); ++i)
        {
            VectorPtr e0i = atoms_[i]->getXYZ() - atoms_[0]->getXYZ(); e0i.normalize();
            double dotprod = fabs(e01.dot(e0i));
            //if the dot product is 1, then the atoms are on the same line
            if ( fabs(dotprod - 1.0) > 1e-12 ) //not linear
            {
                islinear_ = false;
                break;
            }
            islinear_ = true; //all in a line
        }
    }
}

void
Molecule::print(ostream& os) const
{
    os << "Molecule: " << endl;
    os << stream_printf("natoms = %d", natoms()) << endl;
    os << "point group = " << pg_->name() << endl;
    for (int i=0; i < natoms(); ++i)
    {
        getAtom(i)->print(os);
        os << endl;
    }
}

ConstPointGroupPtr
Molecule::getPointGroup() const
{
    return pg_;

}

PointGroupPtr
Molecule::getPointGroup()
{
    return pg_;
}

int
Molecule::natoms() const
{
    return atoms_.size();
}

double
Molecule::getMass(const std::string& symbol)
{
    return MASSES[symbol];
}

ConstAtomPtr
Molecule::getAtom(int num) const
{
    if (num < 0 || num >= (int)atoms_.size()) {
        except(stream_printf("Invalid atom number %d natoms = %d\n", num, atoms_.size()));
    }
    return atoms_[num];
}

ConstAtomPtr
Molecule::getDummy(int num) const
{
    if (num < 0 || num >= dummies_.size()) {
        except(stream_printf("Invalid dummy atom number %d natoms = %d\n", num, dummies_.size()));
    }
    return dummies_[num];
}

int
Molecule::getNumber(const ConstAtomPtr& atom) const
{
    map<const Atom*, int>::const_iterator it = number_map_.find(atom.get());
    return it->second;
}

void 
Molecule::scaleXYZ(double scfactor)
{ 
    for (int i=0; i < atoms_.size(); i++)
        atoms_[i]->scaleXYZ(scfactor);
}

ConstRectMatrixPtr 
Molecule::getXYZ() const
{
    return xyz_;

}

int
Molecule::ndummies() const
{
    return dummies_.size();
}

void
Molecule::computeAbelianSymmetry()
{

    //add the matrix operations
    vector<SymmetryOperationPtr> symmops;
    symmops.push_back(SymmetryOperation::c2x_op());
    symmops.push_back(SymmetryOperation::c2y_op());
    symmops.push_back(SymmetryOperation::c2z_op());
    symmops.push_back(SymmetryOperation::sigmaxy_op());
    symmops.push_back(SymmetryOperation::sigmaxz_op());
    symmops.push_back(SymmetryOperation::sigmayz_op());
    symmops.push_back(SymmetryOperation::inversion_op());

    //determine the symmetry elements of the molecule
    for (int i=0; i < symmops.size(); ++i)
    {
        SymmetryOperationPtr op = symmops[i];
        if (hasSymmOp(op))
        {
            pg_->addOperation(op);
        }
    }
}

void
Molecule::recomputePointGroup(const ConstMoleculePtr& parent)
{
    pg_ = new PointGroup(this);

    if (!needpg_) //point group will never be used
        return;

    vector<ConstSymmetryOperationPtr> symmops; parent->getSymmetryElements(symmops);
    vector<ConstSymmetryOperationPtr>::const_iterator it;
    for (it = symmops.begin(); it != symmops.end(); ++it)
    {
        ConstSymmetryOperationPtr op = *it;
        if (hasSymmOp(op))
            pg_->addOperation(op);
    }
    
    if (parent->getPointGroup()->isAbelian())
    {
        pg_->formAbelianClasses();
        pg_->close();
    }
    else
    {
        pg_->formClosedGroup();
    }

}

void
Molecule::computePointGroup()
{
    string subgrp = KeywordSet::getKeyword("subgroup")->getValueString();
    if (subgrp != "c1")
        computeAbelianSymmetry();
    pg_->formClosedGroup();
}

void
Molecule::addSymmetryOperation(
    const ConstSymmetryOperationPtr& symmop
)
{
    int debug = KeywordSet::getKeyword("symmetry debug")->getValueInteger();
    if (hasSymmOp(symmop))
    {
        if (debug)
        {
            cout << "Adding operation" << endl;
            symmop->print();
        }
        pg_->addOperation(symmop);
    }
    else
    {
        stringstream sstr;
        sstr << "Symmetry operation not valid for this molecule!" << endl;
        symmop->print(sstr);
        symmop->matrix().print("symm op", sstr);
        except(sstr.str());
    }
}

RectMatrixPtr
buildSymmOpMatrix(
    SymmetryOperationPtr op,
    vector<SymmetryInternalCoordinatePtr>& coords
)
{
    RectMatrixPtr m(coords.size(), coords.size());
    for (int i=0; i < coords.size(); ++i)
    {
        VectorPtr Bi = coords[i]->mapBMatrix(op).toVector();  
        for (int j=0; j < coords.size(); ++j)
        {
            VectorPtr Bj = coords[j]->getBMatrix().toVector(); 
            double Bij = Bi.dot(Bj);
            m.set_element(i,j,Bij);
        }
    }
    return m;
}

RectMatrixPtr
coefficientMatrix(
    vector<InternalCoordinatePtr>& coords
)
{
    SymmetryInternalCoordinatePtr symmcoord = boost::dynamic_pointer_cast<SymmetryInternalCoordinate,InternalCoordinate>(coords[0].get());
    int nsimples = symmcoord->getCoefficients().n();
    RectMatrixPtr coefs(nsimples, coords.size());

    for (int i=1; i < coords.size(); ++i)
    {
        symmcoord = boost::dynamic_pointer_cast<SymmetryInternalCoordinate,InternalCoordinate>(coords[i].get());
        coefs.assign_column(symmcoord->getCoefficients(), i);
    }
    return coefs;
}

void
Molecule::getSymmetryElements(
    vector<ConstSymmetryOperationPtr >& symmops
) const
{
    pg_->getSymmetryElements(symmops);
}

bool
Molecule::isLinear() const
{
    return islinear_;
}

int
Molecule::ninternals() const
{
    int threen = 3 * natoms();
    if (atoms_.size() == 1)
        return 0;
    else if (isLinear())
        return threen - 5;
    else
        return threen - 6;
}

bool
Molecule::hasSymmOp(ConstSymmetryOperationPtr op) const
{
    ConstRectMatrixPtr matrix = op->matrix();
    for (int n=0; n < atoms_.size(); n++)
    {
        AtomPtr atom = atoms_[n];
        VectorPtr trans_xyz = matrix * atom->getXYZ();
        AtomPtr test_atom = new Atom(trans_xyz, atom->symbol());
        if (findAtomNumber(test_atom) == -1)
            return false;
    }
    //success on all atoms... this is a symm op
    return true;
}

int
Molecule::findAtomNumber(const ConstAtomPtr& atom) const
{
    double symm_tolerance = KeywordSet::getKeyword("symmetry tolerance")->getValueDouble();
    double tolerance = pow(10, -symm_tolerance);
    for (int n=0; n < natoms(); ++n)
    {
        AtomPtr check_atom = atoms_[n];
        bool symbols_equal = (atom->symbol() == check_atom->symbol());
        VectorPtr xyz_diff = atom->getXYZ() - check_atom->getXYZ();
        bool pos_equal = (xyz_diff.maxabs() < tolerance);
        if (symbols_equal && pos_equal)
            return n;
    }
    //found no matches
    return -1;
}

void
Molecule::displace(ConstRectMatrixPtr disp)
{
    xyz_.accumulate(disp);
    resetAtomXYZ();
}

void
Molecule::resetAtomXYZ()
{
    for (int num=0; num < atoms_.size(); num++)
    {
        atoms_[num]->setXYZ( xyz_.get_row(num) );
    }
}

void 
Molecule::initialize_statics()
{
    if (statics_done)
        return;

    //intialize the masses
    MASSES["H" ] =   1.007825032;
    MASSES["HE"] =   4.002603254;
    MASSES["LI"] =   7.016004548;
    MASSES["BE"] =   9.012182201;
    MASSES["B" ] =  11.009305406;
    MASSES["C" ] =  12.000000000;
    MASSES["N" ] =  14.003074005;
    MASSES["O" ] =  15.994914620;
    MASSES["F" ] =  18.998403224;
    MASSES["NE"] =  19.992440175;
    MASSES["NA"] =  22.989769281;
    MASSES["MG"] =  23.985041699;
    MASSES["AL"] =  26.981538627;
    MASSES["SI"] =  27.976926532;
    MASSES["P" ] =  30.973761629;
    MASSES["S" ] =  31.972070999;
    MASSES["CL"] =  34.968852682;
    MASSES["AR"] =  39.962383123;
    MASSES["K" ] =  38.963706679;
    MASSES["CA"] =  39.962590983;
    MASSES["SC"] =  44.955911909;
    MASSES["TI"] =  47.947946281;
    MASSES["V" ] =  50.943959507;
    MASSES["CR"] =  51.940507472;
    MASSES["MN"] =  54.938045141;
    MASSES["FE"] =  55.934937475;
    MASSES["CO"] =  58.933195048;
    MASSES["NI"] =  57.935342907;
    MASSES["CU"] =  62.929597474;
    MASSES["ZN"] =  63.929142222;
    MASSES["GA"] =  68.925573587;
    MASSES["GE"] =  73.921177767;
    MASSES["AS"] =  74.921596478;
    MASSES["SE"] =  79.916521271;
    MASSES["BR"] =  78.918337087;
    MASSES["KR"] =  85.910610729;
    MASSES["RB"] =  84.911789737;
    MASSES["SR"] =  87.905612124;
    MASSES["Y" ] =  88.905848295;
    MASSES["ZR"] =  89.904704416;
    MASSES["NB"] =  92.906378058;
    MASSES["MO"] =  97.905408169;
    MASSES["TC"] =  98.906254747;
    MASSES["RU"] = 101.904349312;
    MASSES["RH"] = 102.905504292;
    MASSES["PD"] = 105.903485715;
    MASSES["AG"] = 106.90509682 ;
    MASSES["CD"] = 113.90335854 ;
    MASSES["IN"] = 114.903878484;
    MASSES["SN"] = 119.902194676;
    MASSES["SB"] = 120.903815686;
    MASSES["TE"] = 129.906224399;
    MASSES["I" ] = 126.904472681;
    MASSES["XE"] = 131.904153457;
    MASSES["CS"] = 132.905451932;
    MASSES["BA"] = 137.905247237;
    MASSES["LA"] = 138.906353267;
    MASSES["CE"] = 139.905438706;
    MASSES["PR"] = 140.907652769;
    MASSES["ND"] = 144.912749023;
    MASSES["PM"] = 151.919732425;
    MASSES["SM"] = 152.921230339;
    MASSES["EU"] = 157.924103912;
    MASSES["GD"] = 158.925346757;
    MASSES["TB"] = 163.929174751;
    MASSES["DY"] = 164.93032207 ;
    MASSES["HO"] = 165.930293061;
    MASSES["ER"] = 168.93421325 ;
    MASSES["TM"] = 173.938862089;
    MASSES["YB"] = 174.940771819;
    MASSES["LU"] = 179.946549953;
    MASSES["HF"] = 180.947995763;
    MASSES["TA"] = 183.950931188;
    MASSES["W" ] = 186.955753109;
    MASSES["RE"] = 191.96148069 ;
    MASSES["OS"] = 192.96292643 ;
    MASSES["IR"] = 194.964791134;
    MASSES["PT"] = 196.966568662;
    MASSES["AU"] = 201.970643011;
    MASSES["HG"] = 204.974427541;
    MASSES["TL"] = 207.976652071;
    MASSES["PB"] = 208.980398734;
    MASSES["BI"] = 208.982430435;
    MASSES["PO"] = 210.987496271;
    MASSES["AT"] = 222.017577738;
    MASSES["RN"] = 222.01755173 ;
    MASSES["FR"] = 228.031070292;
    MASSES["RA"] = 227.027752127;
    MASSES["AC"] = 232.038055325;
    MASSES["TH"] = 231.03588399 ;
    MASSES["PA"] = 238.050788247;
    MASSES["U" ] = 237.048173444;
    MASSES["NP"] = 242.058742611;
    MASSES["PU"] = 243.06138108 ;
    MASSES["AM"] = 247.07035354 ;
    MASSES["CM"] = 247.07030708 ;
    MASSES["BK"] = 251.079586788;
    MASSES["CF"] = 252.082978512;
    MASSES["ES"] = 257.095104724;
    MASSES["FM"] = 258.098431319;
    MASSES["MD"] = 255.093241131;
    MASSES["NO"] = 260.105504   ;
    MASSES["LR"] = 263.112547   ;
    MASSES["RF"] = 255.107398   ;
    MASSES["DB"] = 259.114500   ;
    MASSES["SG"] = 262.122892   ;
    MASSES["BH"] = 263.128558   ;
    MASSES["HS"] = 265.136151   ;
    MASSES["MT"] = 281.162061   ;
    MASSES["DS"] = 272.153615   ;
    MASSES["RG"] = 283.171792   ;

    needpg_ = false;
    vector<string> moldetails; KeywordSet::getKeyword("mol details")->getValueVectorString(moldetails);
    vector<string>::iterator it;
    for (it = moldetails.begin(); it != moldetails.end(); ++it)
    {
        if (*it == "pg")
            needpg_ = true;
    }


    statics_done = true;
}

MoleculePtr
Molecule::copy() const
{
    RectMatrixPtr newxyz = xyz_.copy();
    vector<string> atom_symbols;
    for (int num=0; num < atoms_.size(); num++)
        atom_symbols.push_back(atoms_[num]->symbol());
    //vector<ConstSymmetryOperationPtr> symmops; getSymmetryElements(symmops);
    MoleculePtr newmol(new Molecule(newxyz, atom_symbols));
    return newmol;
}

void
Molecule::setXYZ(ConstRectMatrixPtr newxyz)
{
    xyz_.assign(newxyz);
    resetAtomXYZ();
}

string
Molecule::getXYZString(const std::string& units, int indent_level) const
{
    string indent = "";
    for (int i=0; i < indent_level; ++i) 
        indent += "\t";

    //set up the xyz coordinates string
    ostringstream xyz_stream;

    //do the first atom
    double unit_conversion = getConversion(units);
    double xcoord = xyz_.get_element(0, X); double ycoord = xyz_.get_element(0, Y); double zcoord = xyz_.get_element(0, Z);
    xyz_stream << stream_printf("%s%16.12f\t%16.12f\t%16.12f", 
                                indent.c_str(), 
                                xcoord * unit_conversion, 
                                ycoord * unit_conversion, 
                                zcoord * unit_conversion
                              );

    //now do the rest of the atoms, adding a newline before each atom
    for (int num=1; num < atoms_.size(); ++num)
    {
        xcoord = xyz_.get_element(num,X); ycoord = xyz_.get_element(num, Y); zcoord = xyz_.get_element(num, Z);
        xyz_stream << stream_printf("\n%s%16.12f\t%16.12f\t%16.12f", 
                                    indent.c_str(), 
                                    xcoord * unit_conversion, 
                                    ycoord * unit_conversion, 
                                    zcoord * unit_conversion
                                   );
    }
    return xyz_stream.str();
}

void
Molecule::displaceAtom(int n, ConstVectorPtr disp)
{
    AtomPtr atom = atoms_[n];
    atom->displaceXYZ(disp);
}

int
Molecule::registerDummyAtom(ConstVectorPtr xyz)
{
    AtomPtr atom = new Atom(xyz, "X");
    dummies_.push_back(atom);
    int atom_num = dummies_.size();
    return atom_num + natoms();
}

double
Molecule::getConversion(const std::string& units)
{
    string bondunits = KeywordSet::getKeyword("bond units")->getValueString();
    if (units ==  bondunits)
        return 1.0;
    else if (units == angstrom and bondunits == bohr)
        return BOHR_TO_ANGSTROM;
    else if (units == bohr and bondunits == angstrom)
        return ANGSTROM_TO_BOHR;
    else
        except(stream_printf("Invalid conversion from %s to %s", bondunits.c_str(), units.c_str()));
}


Atom::Atom(
    ConstVectorPtr xyz, 
    const std::string& symbol
) : 
    symbol_(symbol), 
    mass_(Molecule::getMass(symbol)),
    xyz_(xyz.copy())
{
    SetRuntime(Atom);
}

Atom::Atom(const XMLArchivePtr& arch)
  : Serializable(arch)
{
    SetRuntime(Atom);
    serial_load(xyz);
    serial_load(symbol);
    serial_load(mass);
}

void
Atom::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(xyz);
    serial_save(symbol);
    serial_save(mass);
}

string
Atom::symbol() const
{
    return symbol_;
}

double
Atom::mass() const
{
    return mass_;
}

void
Atom::displaceXYZ(ConstVectorPtr disp)
{
    xyz_.accumulate(disp);
}

ConstVectorPtr
Atom::getXYZ() const
{
    return xyz_;
}

void
Atom::setXYZ(ConstVectorPtr newxyz)
{
    xyz_.assign(newxyz);
}

void
Atom::scaleXYZ(double s)
{
    xyz_.scale(s);
}

void
Atom::print(ostream& os) const
{
    os << stream_printf("%3s %12.8f %12.8f %12.8f",
                        symbol_.c_str(),
                        xyz_[0],
                        xyz_[1],
                        xyz_[2]);
}

