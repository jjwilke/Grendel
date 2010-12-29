#include "gigide.h"

#define X 0
#define Y 1
#define Z 2

using namespace std;
using namespace gigide;
using namespace smartptr;

SerialDeclare(PointGroup);
SerialDeclare(PointGroupClass);
SerialDeclare(Rotation);
SerialDeclare(Reflection);
SerialDeclare(Inversion);
SerialDeclare(ImproperRotation);
SerialDeclare(IdentityElement);


VectorPtr
mapCoordinate(
    SymmetryOperationPtr op,
    InternalCoordinatePtr coord
)
{
    return coord->mapBMatrix(op).toVector();
}

SymmetryOperation::OperationType
PointGroupClass::type() const
{
    return type_;
}

SymmetryOperation::SymmetryOperation(
    OperationType type
) : type_(type)
{
}

SymmetryOperation::SymmetryOperation(
    const XMLArchivePtr& arch
) : Serializable(arch)
{
    serial_load_enum(type);
}

void
SymmetryOperation::serialize(
    const XMLArchivePtr& arch
) const
{
    Serializable::serialize(arch);
    serial_save_enum(type);
}

void
SymmetryOperation::print(ostream& os) const
{
    os << description();
}

bool
SymmetryOperation::isEquivalent(
    const ConstSymmetryOperationPtr& oper
) const
{
    double exp = KeywordSet::getKeyword("symmetry tolerance")->getValueDouble();
    double tol = pow(10, -exp);
    return equals(matrix_, oper->matrix(), tol);
}

SymmetryOperation::OperationType
SymmetryOperation::type() const
{
    return type_;
}

ConstRectMatrixPtr
SymmetryOperation::matrix() const
{
    return matrix_;
}

RectMatrixPtr
SymmetryOperation::getPermutationMatrix(const ConstMoleculePtr& mol) const
{
    int natoms = mol->natoms();
    RectMatrixPtr perm(natoms, natoms);
    for (int n=0; n < natoms; ++n)
    {
        ConstAtomPtr atom = mol->getAtom(n);
        VectorPtr trans_xyz = matrix_ * atom->getXYZ();
        AtomPtr test_atom = new Atom(trans_xyz, atom->symbol());
        int mapped_number = mol->findAtomNumber(test_atom);
        //sanity check
        //if this isn't really a symmetry operation, complain
        if (mapped_number == -1) 
        {
            atom->getXYZ().print("xyz", cerr);
            trans_xyz.print("trans xyz", cerr);
            matrix_.print("symmop", cerr);
            except(stream_printf("Symmetry operation %s is not valid for molecule:", name().c_str()));
        }
        
        //set the permutation elements
        perm.set_element(n,mapped_number, 1.0);
    }
    return perm;
}

string
SymmetryOperation::name() const
{
    return name_;
}

string
SymmetryOperation::opname(SymmetryOperation::OperationType type)
{
    switch (type)
    {
        case rotation:
            return "rotation";
        case reflection:
            return "reflection";
        case improper_rotation:
            return "improper rotation";
        case identity:
            return "identity";
        case inversion:
            return "inversion";
        default:
            cerr << "invalid symmetry operation type" << endl;
            abort();
    }
}

string
SymmetryOperation::axis_str(const VectorPtr& vec)
{
    stringstream sstr;
    sstr << stream_printf("[%15.12f %15.12f %15.12f]",
                     vec.get_element(X),
                     vec.get_element(Y),
                     vec.get_element(Z)
                    );
    return sstr.str();
}

SymmetryOperationPtr
SymmetryOperation::multiply(
    const SymmetryOperationPtr& r
) const
{
    RectMatrixPtr prod = matrix_ * r->matrix();
    SymmetryOperationPtr op = buildOperation(prod);

    if (!equals(op->matrix(), prod)) //error
    {
        stringstream sstr;
        sstr << "Product = " << endl;
        op->print(sstr);
        op->matrix().print("Wrong matrix", sstr);
        prod.print("Correct matrix", sstr);
        sstr << "Symmetry operation multiplication has failed.  L * R != Product" << endl;
        sstr << "L = " << endl;
        print(sstr);
        matrix_.print("L Matrix", sstr);
        sstr << "R = " << endl;
        r->print(sstr);
        r->matrix().print("R Matrix", sstr);
        except(sstr.str());
    }

    return op;
}

SymmetryOperationPtr
SymmetryOperation::buildOperation(
    ConstRectMatrixPtr m
)
{
    SymmetryOperationPtr op;
    RectMatrixPtr L, R;
    VectorPtr revals, ievals;
    m.eigen(revals, ievals, L, R);

    int nplus1 = 0;
    int nminus1 = 0;
    for (int i=0; i < 3; ++i)
    {
        double eval = revals[i];
        if      ( fabs(eval - 1.0) < 1e-8 )
            ++nplus1;
        else if ( fabs(eval + 1.0) < 1e-8 )
            ++nminus1;
    }

    //depending on the eigenvalue types, we have different types of operations
    if (nplus1 == 3) //identity
    {
        op = new IdentityElement();
    }
    else if (nminus1 == 3)
    {
        op = new Inversion();
    }
    else if (nminus1 == 1 && nplus1 == 2)
    {
        op = buildReflection(
                revals,
                R,
                m
             );
    }
    else if (nminus1 == 1 && nplus1 == 0)
    {
        op = buildImproperRotation(
                revals,
                R,
                m
             );
    }
    else if (nplus1 == 1)
    {
        op = buildRotation(
                revals,
                R,
                m 
             ); 
    }
    else
    {
        stringstream sstr;
        sstr << "Matrix does not correspond to symmetry operation" << endl;
        m.print("Matrix", sstr);
        revals.print("Eigenvalues", sstr);
        except(sstr.str());
    }

    if (!equals(op->matrix(), m)) //error
    {
        stringstream sstr;
        sstr << "Symmetry operation build has failed" << endl;
        op->print(sstr);
        op->matrix().print("Wrong matrix", sstr);
        m.print("Correct matrix", sstr);
        except(sstr.str());
    }

    return op;
}

VectorPtr
getAxis(
    ConstVectorPtr evals,
    ConstRectMatrixPtr axes,
    double axiseval,
    double& theta
)
{
    int axis;
    int offaxis;
    for (int i=0; i < 3; ++i)
    {
        if ( fabs(evals[i] - axiseval) < 1e-8 )
            axis = i;
        else
            offaxis = i;
    }

    if (axis == -1)
        except("Unable to find axis");
    if (offaxis == -1)
        except("Unable to find off axis vector");

    theta = inverse_cosine(evals[offaxis]);
    return axes.get_column(axis);
}

bool
getOrderAndExponent(
    double theta,
    int& order,
    int& exponent
)
{
    double tho2pi = theta / 2.0 / PI;
    bool check = getFraction(tho2pi, exponent, order);
    return check;
}

RotationPtr
SymmetryOperation::buildRotation(
    ConstVectorPtr evals,
    ConstRectMatrixPtr evecs,
    ConstRectMatrixPtr m
)
{
    double theta;
    int order, exponent;
    VectorPtr axis = getAxis(evals, evecs, 1.0, theta);
    bool check = getOrderAndExponent(theta, order, exponent);
    RotationPtr op = new Rotation(axis, order, exponent);
    if (!equals(op->matrix(), m, 1e-8)) //error, currently the sign can get effed up
    {
        op = new Rotation(axis, order, order - exponent);
    }
    return op;
}

ReflectionPtr
SymmetryOperation::buildReflection(
    ConstVectorPtr evals,
    ConstRectMatrixPtr evecs,
    ConstRectMatrixPtr m
)
{
    double theta;
    int order, exponent;
    ConstVectorPtr axis = getAxis(evals, evecs, -1.0, theta);
    ReflectionPtr op = new Reflection(axis);
    return op;
}

ImproperRotationPtr
SymmetryOperation::buildImproperRotation(
    ConstVectorPtr evals,
    ConstRectMatrixPtr evecs,
    ConstRectMatrixPtr m
)
{
    double theta;
    int order, exponent;
    ConstVectorPtr axis = getAxis(evals, evecs, -1.0, theta);
    bool check = getOrderAndExponent(theta, order, exponent);
    ImproperRotationPtr op = new ImproperRotation(axis, order, exponent);
    if (!equals(op->matrix(), m, 1e-8)) //error, currently the sign can get effed up
    {
        op = new ImproperRotation(axis, order, order - exponent);
    }
    return op;
}

string
SymmetryOperation::description() const
{
    return name_;
}

SymmetryOperationPtr
SymmetryOperation::conjugate(const ConstSymmetryOperationPtr& oper) const
{
    RectMatrixPtr temp = matrix_ * oper->matrix().t();
    RectMatrixPtr m = oper->matrix() * temp;
    SymmetryOperationPtr op = buildOperation(m);
    return op;
}

void
SymmetryOperation::init_statics()
{
    if (initdone_)
        return;

    //now build the vector of symm operations
    VectorPtr xaxis = new Vector(3); xaxis.set_element(X, 1.0);
    VectorPtr yaxis = new Vector(3); yaxis.set_element(Y, 1.0);
    VectorPtr zaxis = new Vector(3); zaxis.set_element(Z, 1.0);

    c2x_op_ = new Rotation(ConstVectorPtr(xaxis), 2);
    c2y_op_ = new Rotation(ConstVectorPtr(yaxis), 2);
    c2z_op_ = new Rotation(ConstVectorPtr(zaxis), 2);
    inversion_op_ = new Inversion();
    sigmaxy_op_ = new Reflection(ConstVectorPtr(zaxis));
    sigmaxz_op_ = new Reflection(ConstVectorPtr(yaxis));
    sigmayz_op_ = new Reflection(ConstVectorPtr(xaxis));
    identity_op_ = new IdentityElement();

    initdone_ = true;
}

SymmetryOperationPtr
SymmetryOperation::identity_op()
{
    if (!initdone_)
        except("Static initialization of symmetry operation not yet completed");

    return identity_op_;
}

SymmetryOperationPtr
SymmetryOperation::inversion_op()
{
    if (!initdone_)
        except("Static initialization of symmetry operation not yet completed");

    return inversion_op_;
}

SymmetryOperationPtr
SymmetryOperation::c2x_op()
{
    if (!initdone_)
        except("Static initialization of symmetry operation not yet completed");

    return c2x_op_;
}

SymmetryOperationPtr
SymmetryOperation::c2y_op()
{
    if (!initdone_)
        except("Static initialization of symmetry operation not yet completed");

    return c2y_op_;
}

SymmetryOperationPtr
SymmetryOperation::c2z_op()
{
    if (!initdone_)
        except("Static initialization of symmetry operation not yet completed");

    return c2z_op_;
}

SymmetryOperationPtr
SymmetryOperation::sigmaxy_op()
{
    if (!initdone_)
        except("Static initialization of symmetry operation not yet completed");

    return sigmaxy_op_;
}

SymmetryOperationPtr
SymmetryOperation::sigmaxz_op()
{
    if (!initdone_)
        except("Static initialization of symmetry operation not yet completed");

    return sigmaxz_op_;
}

SymmetryOperationPtr
SymmetryOperation::sigmayz_op()
{
    if (!initdone_)
        except("Static initialization of symmetry operation not yet completed");

    return sigmayz_op_;
}

Rotation::Rotation(
    ConstVectorPtr axis,
    int order,
    int exponent
) : SymmetryOperation(SymmetryOperation::rotation),
    order_(order),
    exponent_(exponent)
{
    SetRuntime(Rotation);
    axis_ = axis.copy(); axis_.normalize();
    init();
}

void
Rotation::init()
{
    double phi = (2.0 * PI / order_) * exponent_;
    double a = axis_.get_element(0); //alpha
    double b = axis_.get_element(1); //beta
    double g = axis_.get_element(2); //gamma
    double cosphi = cos(phi);
    double sinphi = sin(phi);

    if (exponent_ >= order_) //nope, not allowed
        except(stream_printf("Exponent %d is not valid for rotation order %d", exponent_, order_));

    double Rxx = a*a + (b*b + g*g) * cosphi;
    double Rxy = a*b*(1 - cosphi) + g*sinphi;
    double Rxz = a*g*(1 - cosphi) - b*sinphi;

    double Ryx = a*b*(1 - cosphi) - g*sinphi;
    double Ryy = b*b + (a*a + g*g) * cosphi;
    double Ryz = b*g*(1-cosphi) + a*sinphi;

    double Rzx = a*g*(1-cosphi) + b*sinphi;
    double Rzy = b*g*(1-cosphi) - a*sinphi;
    double Rzz = g*g + (a*a + b*b)*cosphi;

    matrix_ = RectMatrixPtr(3,3);
    matrix_.set_element(X,X,Rxx);
    matrix_.set_element(X,Y,Rxy);
    matrix_.set_element(X,Z,Rxz);
    matrix_.set_element(Y,X,Ryx);
    matrix_.set_element(Y,Y,Ryy);
    matrix_.set_element(Y,Z,Ryz);
    matrix_.set_element(Z,X,Rzx);
    matrix_.set_element(Z,Y,Rzy);
    matrix_.set_element(Z,Z,Rzz);

    stringstream sstr;
    sstr << "C" << order_;

    if      (fabs(axis_.get_element(X) - 1.0) < 1e-8) sstr << "X"; 
    else if (fabs(axis_.get_element(Y) - 1.0) < 1e-8) sstr << "Y";
    else if (fabs(axis_.get_element(Z) - 1.0) < 1e-8) sstr << "Z";
    else if (fabs(axis_.get_element(X) + 1.0) < 1e-8) sstr << "X";
    else if (fabs(axis_.get_element(Y) + 1.0) < 1e-8) sstr << "Y";
    else if (fabs(axis_.get_element(Z) + 1.0) < 1e-8) sstr << "Z";

    name_ = sstr.str();
}

Rotation::Rotation(const XMLArchivePtr& arch)
    : SymmetryOperation(arch)
{
    SetRuntime(Rotation);
    serial_load(order);
    serial_load(exponent);
    serial_load(axis);
    init();
}

void
Rotation::serialize(const XMLArchivePtr& arch) const
{
    SymmetryOperation::serialize(arch);
    serial_save(order);
    serial_save(exponent);
    serial_save(axis);
}

string
Rotation::description() const
{
    string descr = name() + " " + axis_str(axis_);
    return descr;
}

int
Rotation::order() const
{
    return order_;
}

int
Rotation::exponent() const
{
    return exponent_;
}

const VectorPtr&
Rotation::axis() const
{
    return axis_;
}

ImproperRotation::ImproperRotation(
    ConstVectorPtr axis,
    int order,
    int exponent
) : SymmetryOperation(SymmetryOperation::improper_rotation),
    order_(order),
    exponent_(exponent),
    axis_(axis.copy())
{
    SetRuntime(ImproperRotation);
    axis_.normalize();
    init();
}

void
ImproperRotation::init()
{
    double phi = (2.0 * PI / order_) * exponent_;
    double a = axis_.get_element(0); //alpha
    double b = axis_.get_element(1); //beta
    double g = axis_.get_element(2); //gamma
    double cosphi = cos(phi);
    double sinphi = sin(phi);

    double Rxx = -a*a + (b*b + g*g) * cosphi;
    double Rxy = -a*b*(1 + cosphi) + g*sinphi;
    double Rxz = -a*g*(1 + cosphi) - b*sinphi;

    double Ryx = -a*b*(1 + cosphi) - g*sinphi;
    double Ryy = -b*b + (a*a + g*g) * cosphi;
    double Ryz = -b*g*(1+cosphi) + a*sinphi;

    double Rzx = -a*g*(1+cosphi) + b*sinphi;
    double Rzy = -b*g*(1+cosphi) - a*sinphi;
    double Rzz = -g*g + (a*a + b*b)*cosphi;

    matrix_ = RectMatrixPtr(3,3);
    matrix_.set_element(X,X,Rxx);
    matrix_.set_element(X,Y,Rxy);
    matrix_.set_element(X,Z,Rxz);
    matrix_.set_element(Y,X,Ryx);
    matrix_.set_element(Y,Y,Ryy);
    matrix_.set_element(Y,Z,Ryz);
    matrix_.set_element(Z,X,Rzx);
    matrix_.set_element(Z,Y,Rzy);
    matrix_.set_element(Z,Z,Rzz);

    stringstream sstr;
    sstr << "S" << order_;
    name_ = sstr.str();
}

ImproperRotation::ImproperRotation(const XMLArchivePtr& arch)
    : SymmetryOperation(arch)
{
    SetRuntime(ImproperRotation);
    serial_load(order);
    serial_load(exponent);
    serial_load(axis);
    init();
}

void
ImproperRotation::serialize(const XMLArchivePtr& arch) const
{
    SymmetryOperation::serialize(arch);
    serial_save(order);
    serial_save(exponent);
    serial_save(axis);
}

ConstVectorPtr
ImproperRotation::axis() const
{
    return axis_;
}

string
ImproperRotation::description() const
{
    string descr = name() + " " + axis_str(axis_);
    return descr;
}

int
ImproperRotation::order() const
{
    return order_;
}

int
ImproperRotation::exponent() const
{
    return exponent_;
}

Reflection::Reflection(
    ConstVectorPtr axis
) : SymmetryOperation(SymmetryOperation::reflection)
{
    SetRuntime(Reflection);
    init(axis);
}

Reflection::Reflection(
   const XYZPointPtr& pt1,
   const XYZPointPtr& pt2,
   const XYZPointPtr& pt3
) : SymmetryOperation(SymmetryOperation::reflection)
{
    SetRuntime(Reflection);
    VectorPtr vec1 = pt3->getXYZ() - pt1->getXYZ();
    VectorPtr vec2 = pt2->getXYZ() - pt1->getXYZ();
    VectorPtr axis = cross(vec1,vec2);
    init(axis);
}

Reflection::Reflection(const XMLArchivePtr& arch)
    : SymmetryOperation(arch)
{
    SetRuntime(Reflection);

    boost::intrusive_ptr<Vector> tmp;
    serial_call_load(arch, tmp, "axis");
    VectorPtr axis(tmp.get()); 
    init(axis);
}

void
Reflection::serialize(const XMLArchivePtr& arch) const
{
    SymmetryOperation::serialize(arch);
    serial_save(axis);
}

const VectorPtr&
Reflection::axis() const
{
    return axis_;
}

void
Reflection::init(
    ConstVectorPtr axis
)
{
    axis_ = axis.copy(); axis_.normalize();

    double a = axis_.get_element(0); //alpha
    double b = axis_.get_element(1); //beta
    double g = axis_.get_element(2); //gamma

    double Rxx = 1 - 2*a*a;
    double Rxy = -2*a*b;
    double Rxz = -2*a*g;

    double Ryx = -2*a*b;
    double Ryy = 1 - 2*b*b;
    double Ryz = -2*b*g;

    double Rzx = -2*a*g;
    double Rzy = -2*b*g;
    double Rzz = 1 - 2*g*g;

    matrix_ = RectMatrixPtr(3,3);
    matrix_.set_element(X,X,Rxx);
    matrix_.set_element(X,Y,Rxy);
    matrix_.set_element(X,Z,Rxz);
    matrix_.set_element(Y,X,Ryx);
    matrix_.set_element(Y,Y,Ryy);
    matrix_.set_element(Y,Z,Ryz);
    matrix_.set_element(Z,X,Rzx);
    matrix_.set_element(Z,Y,Rzy);
    matrix_.set_element(Z,Z,Rzz);

    stringstream sstr;
    sstr << "SIGMA";

    if      (fabs(axis_.get_element(X) - 1.0) < 1e-8) sstr << "YZ"; 
    else if (fabs(axis_.get_element(Y) - 1.0) < 1e-8) sstr << "XZ";
    else if (fabs(axis_.get_element(Z) - 1.0) < 1e-8) sstr << "XY";
    else if (fabs(axis_.get_element(X) + 1.0) < 1e-8) sstr << "YZ";
    else if (fabs(axis_.get_element(Y) + 1.0) < 1e-8) sstr << "XZ";
    else if (fabs(axis_.get_element(Z) + 1.0) < 1e-8) sstr << "XY";

    name_ = sstr.str();
}

string
Reflection::description() const
{
    string descr = name() + " " + axis_str(axis_);
    return descr;
}

Inversion::Inversion() 
    : SymmetryOperation(SymmetryOperation::inversion)
{
    SetRuntime(Inversion);
    init();
}

void
Inversion::init()
{
    matrix_ = RectMatrixPtr(3,3);
    matrix_.set_element(X,X,-1.0);
    matrix_.set_element(Y,Y,-1.0);
    matrix_.set_element(Z,Z,-1.0);
    name_ = "I";
}

Inversion::Inversion(const XMLArchivePtr& arch)
    : SymmetryOperation(arch)
{
    SetRuntime(Inversion);
    init();
}

void
Inversion::serialize(const XMLArchivePtr& arch) const
{
    SymmetryOperation::serialize(arch);
}

IdentityElement::IdentityElement() 
    : SymmetryOperation(SymmetryOperation::identity)
{
    SetRuntime(IdentityElement);
    init();
}

void
IdentityElement::init()
{
    matrix_ = RectMatrixPtr(3,3);
    matrix_.set_element(X,X,1.0);
    matrix_.set_element(Y,Y,1.0);
    matrix_.set_element(Z,Z,1.0);
    name_ = "E";
}

IdentityElement::IdentityElement(const XMLArchivePtr& arch)
    : SymmetryOperation(arch)
{
    SetRuntime(IdentityElement);
    init();
}

void
IdentityElement::serialize(const XMLArchivePtr& arch) const
{
    SymmetryOperation::serialize(arch);
}

PointGroup::PointGroup(
    const ConstMoleculePtr& mol,
    string title
) : 
    mol_(mol),
    pgstr_(title)
{
    SetRuntime(PointGroup);
    valid_ = false;
    classescomputed_ = false;
    isabelian_ = false;
    init_statics();
    SymmetryOperation::init_statics();

    //add at least the identity element
    addOperation(SymmetryOperation::identity_op());

}

PointGroup::PointGroup(const XMLArchivePtr& arch)
    : Serializable(arch)
{
    SetRuntime(PointGroup);

    SymmetryOperation::init_statics();

    serial_load(valid);
    serial_load(classescomputed);
    serial_load(pgstr);
    serial_load(isabelian);
    serial_load(mol);
    serial_load(symmops);
    serial_load(subgroups);
    serial_load(subspaces);

    close();
}

PointGroup::~PointGroup()
{
}

void
PointGroup::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(valid);
    serial_save(classescomputed);
    serial_save(pgstr);
    serial_save(isabelian);
    serial_save(mol);
    serial_save(symmops);
    serial_save(subgroups);
    serial_save(subspaces);
}

int
PointGroup::order() const
{
    return symmops_.size();
}

int
PointGroup::nclasses() const
{
    return classes_.size();
}

void
PointGroup::printClasses(ostream& os) const
{
    vector< PointGroupClassPtr >::const_iterator it;
    for (it = classes_.begin(); it != classes_.end(); ++it)
        os << stream_printf("%8s", (*it)->getRepresentative()->name().c_str());
}

void
PointGroup::print(ostream& os) const
{
    os << "Point Group " << pgstr_ << " of order " << order();
    vector<ConstSymmetryOperationPtr> symmops; getSymmetryElements(symmops);
    vector<ConstSymmetryOperationPtr>::iterator it;
    for (it = symmops.begin(); it != symmops.end(); ++it)
    {
        os << endl;
        (*it)->print(os);
    }
}

bool
PointGroup::isTotallySymmetric(const VectorPtr& chars) const
{
    if (chars.n() != classes_.size())
    {
        cerr << stream_printf("Coordinate has characters for %d classes, but there are %d classes.",
                         chars.n(),
                         classes_.size()) << endl;
        abort();
    }

    double dotprod = 0;
    for (int i=0; i < classes_.size(); ++i)
    {
        dotprod += chars.get_element(i) * classes_[i]->order();
    }

    if ( fabs(dotprod) > 1e-5 ) //totally symmetric
        return true;
    else
        return false;
}

void
PointGroup::getSymmetryElements(
    vector<ConstSymmetryOperationPtr>& symmops
) const
{
    if (!valid_) except("Point group accessor method called before point group has been closed");

    vector<PointGroupClassPtr>::const_iterator it;
    vector<PointGroupClassPtr> classes; getClasses(classes);
    for (it = classes.begin(); it != classes.end(); ++it)
        (*it)->getSymmetryElements(symmops);
}

void
PointGroup::getClasses(vector<PointGroupClassPtr>& classvec) const
{
    if (!valid_) except("Point group accessor method called before point group has been closed");

    classvec = classes_;
}

void
PointGroup::init_statics()
{
    if (initdone_)
        return;

    symmplanes_["SIGMAXY"] = 1;
    symmplanes_["SIGMAYZ"] = 1;
    symmplanes_["SIGMAXZ"] = 1;
    symmaxes_["C2X"] = 1;
    symmaxes_["C2Y"] = 1;
    symmaxes_["C2Z"] = 1;

    initdone_ = true;
}

bool
PointGroup::hasOperation(ConstRectMatrixPtr op) const
{
    if (!valid_) except("Point group accessor method called before point group has been closed");

    double exp = KeywordSet::getKeyword("symmetry tolerance")->getValueDouble();
    double tol = pow(10, -exp);

    vector<ConstSymmetryOperationPtr>::const_iterator it;
    for (it = symmops_.begin(); it != symmops_.end(); ++it)
    {
        ConstSymmetryOperationPtr refop = *it;
        //matrices are equal
        if (equals(refop->matrix(), op, tol))
            return true;
    }
    return false;
}

bool
PointGroup::hasOperation(const ConstSymmetryOperationPtr& op) const
{
    if (!valid_) except("Point group accessor method called before point group has been closed");

    return hasOperation(op->matrix());
}

void
PointGroup::recompute()
{
}

PointGroupClass::PointGroupClass(
    const ConstMoleculePtr& mol,
    SymmetryOperation::OperationType type
) : mol_(mol), type_(type)
{
    SetRuntime(PointGroupClass);
}

PointGroupClass::PointGroupClass(
    const XMLArchivePtr& arch
) 
{
    SetRuntime(PointGroupClass);
    serial_load_enum(type);
    serial_load(members);
}

void
PointGroupClass::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save_enum(type);
    serial_save(members);
}

int
PointGroupClass::order() const
{
    return members_.size();
}

void
PointGroupClass::print(ostream& os) const
{
    os << stream_printf("Class of order %d for operation type %s", order(), SymmetryOperation::opname(type_).c_str());
    for (int i=0; i < order(); ++i)
    {
        os << endl;
        members_[i]->print(os);
    }
}

void
PointGroupClass::getSymmetryElements(
    vector<ConstSymmetryOperationPtr>& symmops
) const
{
    Set<ConstSymmetryOperationPtr>::const_iterator it;
    for (it = members_.begin(); it != members_.end(); ++it)
        symmops.push_back(*it);
}

ConstSymmetryOperationPtr
PointGroupClass::getRepresentative() const
{
    if (members_.size() == 0) 
        except("Point group class is empty and should not be");

    return members_[0];
}

bool
PointGroupClass::hasOperation(const ConstSymmetryOperationPtr& oper) const
{
    //check 
    Set<ConstSymmetryOperationPtr>::const_iterator it;
    for (it = members_.begin(); it != members_.end(); ++it)
    {
        ConstSymmetryOperationPtr member = *it;
        if (member->isEquivalent(oper))
            return true; //already have it
    }
    return false;
}

void
PointGroupClass::addOperation(const ConstSymmetryOperationPtr& oper)
{
    if (!hasOperation(oper))
    {
        members_.append(oper);
    }
}

bool
PointGroupClass::testMembership(
    const ConstSymmetryOperationPtr& op,
    const Set<ConstSymmetryOperationPtr>& symmops
) const
{
    for (int i=0; i < symmops.size(); ++i)
    {
        ConstSymmetryOperationPtr conjugating_op = symmops[i];
        SymmetryOperationPtr conj = op->conjugate(conjugating_op);
        if (hasOperation(conj))
            return true;
    }
    return false;
}

void
PointGroup::addOperation(const ConstSymmetryOperationPtr& op)
{
    //check 
    vector<ConstSymmetryOperationPtr>::const_iterator it;
    for (it = symmops_.begin(); it != symmops_.end(); ++it)
    {
        ConstSymmetryOperationPtr member = *it;
        if (member->isEquivalent(op))
            return; //already have it
    }

    symmops_.push_back(op);
}

Set<PointGroupPtr>
PointGroup::getSubgroups() const
{
    return subgroups_;
}

void
PointGroup::addSubgroup(
    PointGroupPtr pg
)
{
    pg->formClosedGroup();
    subgroups_.push_back(pg);
}

bool
PointGroup::findClass(const ConstSymmetryOperationPtr& op)
{
    vector<PointGroupClassPtr>::const_iterator itclass;
    for (itclass = classes_.begin(); itclass != classes_.end(); ++itclass)
    {
        bool found = (*itclass)->testMembership(op, symmops_);
        if (found) 
        {
            (*itclass)->addOperation(op);
            return true;
        }
    }
    return false; //nothing found
}

bool
PointGroup::isAbelian() const
{
    return isabelian_;
}

void
PointGroup::formClasses()
{
    if (classescomputed_) 
        return;

    vector<ConstSymmetryOperationPtr>::iterator it;
    for (it = symmops_.begin(); it != symmops_.end(); ++it)
    {
        bool found = findClass(*it);
        if (!found) //this operation needs to start a new class
        {
            addClass(*it);
        }
    }

    isabelian_ = true;
    for (int i=0; i < classes_.size(); ++i)
    {
        if (classes_[i]->order() != 1)
            isabelian_ = false;
    }

    classescomputed_ = true;
}

void
PointGroup::formAbelianClasses()
{
    if (classescomputed_) 
        return;

    vector<ConstSymmetryOperationPtr>::iterator it;
    for (it = symmops_.begin(); it != symmops_.end(); ++it)
    {
        addClass(*it);
    }

    isabelian_ = true;
    classescomputed_ = true;
}

void
PointGroup::addClass(const ConstSymmetryOperationPtr& op)
{
    PointGroupClassPtr cls(new PointGroupClass(mol_, op->type())); 
    cls->addOperation(op);
    classes_.push_back(cls);
}

void
PointGroup::close()
{
    valid_ = true; //accessor methods can now be called... this must be at the beginning!
    
    int nplanes = 0;
    int naxes = 0;
    bool hasinv = false;
    vector<ConstSymmetryOperationPtr>::iterator iti;
    for (iti = symmops_.begin(); iti != symmops_.end(); ++iti)
    {
        ConstSymmetryOperationPtr op(*iti);
        if (op->type() == SymmetryOperation::inversion) hasinv = true;
        else if (op->type() == SymmetryOperation::rotation) ++naxes;
        else if (op->type() == SymmetryOperation::reflection) ++nplanes;
    }

    if (pgstr_.size() != 0) //we already have a name
        return;

    pgstr_ = "C1";
    if (nplanes >= 3 && naxes >= 3 && hasinv)
        pgstr_ = "D2H";
    else if (nplanes >= 1 && naxes >= 1 && hasinv)
        pgstr_ = "C2H";
    else if (nplanes >= 0 && naxes >= 3 && !hasinv)
        pgstr_ = "D2";
    else if (nplanes >= 2 && naxes >= 1 && !hasinv)
        pgstr_ = "C2V";
    else if      (nplanes == 1 && naxes == 0 && !hasinv)
        pgstr_ = "CS";
    else if (nplanes >= 0 && naxes >= 0 && hasinv)
        pgstr_ = "CI";
    else if (nplanes >= 0 && naxes >= 1 && !hasinv)
        pgstr_ = "C2";
    else
        pgstr_ = "C1";
}

bool
PointGroup::computeAbelianness()
{
    for (int i=0; i < classes_.size(); ++i)
    {
        ConstSymmetryOperationPtr opi = classes_[i]->getRepresentative();
        for (int j=i+1; j < classes_.size(); ++j)
        {
            ConstSymmetryOperationPtr opj = classes_[j]->getRepresentative();
            
            RectMatrixPtr ab = opi->matrix() * opj->matrix();
            RectMatrixPtr ba = opj->matrix() * opi->matrix();
            if (!equals(ab,ba))
                return false; //not abelian
        }
    }
    return true;
}

void
PointGroup::formClosedGroup()
{
    valid_ = true;

    //iterate conjugation until the group is closed under multiplication
    int lastorder = 0;
    int currentorder = order();
    while (lastorder != currentorder)
    {
        lastorder = currentorder;

        //first attempt to add inverses of everything 
        int size = symmops_.size();
        for (int i=0; i < size; ++i)
        {
            ConstRectMatrixPtr inv = symmops_[i]->matrix().t();
            if (!hasOperation(inv))
            {
                SymmetryOperationPtr newop = SymmetryOperation::buildOperation(inv);
                symmops_.push_back(newop);
            }
        }

        //multiply all elements by each other
        for (int i=0; i < size; ++i)
        {
            for (int j=0; j < size; ++j)
            {
                ConstRectMatrixPtr product = symmops_[i]->matrix() * symmops_[j]->matrix();
                if (!hasOperation(product))
                {
                    SymmetryOperationPtr newop = SymmetryOperation::buildOperation(product);
                    symmops_.push_back(newop);
                }
            } 
        }
        currentorder = order();
    }

    formClasses();

    isabelian_ = computeAbelianness();

    close();
}

void
PointGroup::formBasis(
    const Set<ConstSymmetryInternalCoordinatePtr>& input_internals,
    vector<CoordinateSubspacePtr>& subspaces
) const
{
    Set<ConstSimpleInternalCoordinatePtr> simples = input_internals[0]->getSimples();

    vector<SymmetryInternalCoordinatePtr> internals;

    //create a new set of non-const internals
    for (int i=0; i < input_internals.size(); ++i)
    {
        internals.push_back(input_internals[i]->carbon_copy());
    }

    double exp = KeywordSet::getKeyword("symmetry tolerance")->getValueDouble();
    double tol = pow(10, -exp);

    //decouple the coordinate as much as possible for all of the symmetry operations
    for (int n=0; n < symmops_.size(); ++n)
    {
        ConstSymmetryOperationPtr op = symmops_[n];
        for (int i=0; i < internals.size(); ++i)
        {
            for (int j=0; j < internals.size(); ++j)
            {
                if (i==j)
                    continue;

                VectorPtr Bi_row = internals[i]->mapBMatrix(op).toVector();
                VectorPtr Bi_col = internals[i]->getBMatrix().toVector();
                VectorPtr Bj_row = internals[j]->mapBMatrix(op).toVector();
                VectorPtr Bj_col = internals[j]->getBMatrix().toVector();

                double a12 = Bi_row.dot(Bj_col);

                if ( fabs(a12) < tol ) //already decoupled
                    continue;

                double a22 = Bj_row.dot(Bj_col);
                double a11 = Bi_row.dot(Bi_col);

                SymmMatrixPtr m(2);
                m.set_element(0,0,a11); 
                m.set_element(0,1,a12);
                m.set_element(1,1,a22);

                VectorPtr evals;
                RectMatrixPtr evecs;
                m.eigen(evals, evecs);

                //these coordinates can be decoupled
                ConstVectorPtr Ci = internals[i]->getCoefficients();
                ConstVectorPtr Cj = internals[j]->getCoefficients();

                VectorPtr newCi = evecs.get_element(0,0) * Ci + evecs.get_element(1,0) * Cj;
                VectorPtr newCj = evecs.get_element(0,1) * Ci + evecs.get_element(1,1) * Cj;

                SymmetryInternalCoordinatePtr newVi = new SymmetryInternalCoordinate(newCi, simples, mol_);
                SymmetryInternalCoordinatePtr newVj = new SymmetryInternalCoordinate(newCj, simples, mol_);

                internals[i] = newVi;
                internals[j] = newVj;
            }
        }
    }

    RectMatrixPtr overlap(internals.size(), internals.size());
    for (int n=0; n < symmops_.size(); ++n)
    {
        ConstSymmetryOperationPtr op = symmops_[n];
        for (int i=0; i < internals.size(); ++i)
        {
            for (int j=0; j < internals.size(); ++j)
            {
                VectorPtr Bi = internals[i]->mapBMatrix(op).toVector();
                VectorPtr Bj = internals[j]->getBMatrix().toVector();
                double Sij = fabs(Bi.dot(Bj));
                overlap.accumulate_element(i,j,Sij);
            }
        }
    }

    map<int,int> popped; //keep track of which internal coordinates we have already used
    for (int i=0; i < internals.size(); ++i)
    {
        if (popped[i]) //this coordinate has already been used
            continue;


        CoordinateSubspacePtr subspace(new CoordinateSubspace(mol_));
        subspaces.push_back(subspace);
        subspace->addCoordinate(internals[i]);
        for (int j=i+1; j < internals.size(); ++j)
        {
            double Sij = overlap.get_element(i,j);
            if (Sij > 1e-6) //these coordinates are coupled
            {
                popped[j] = 1;
                subspace->addCoordinate(internals[j]);
            }
        }
    }

    //now loop the coordinate subspaces and attempt to adapt them to subgroup symmetry
    for (int space=0; space < subspaces.size(); ++space)
    {
        CoordinateSubspacePtr subspace = subspaces[space];
        for (int sg=0; sg < subgroups_.size(); ++sg)
        {
            PointGroupPtr pg = subgroups_[sg]; 
            subspace->abelianify(pg);
        }
    }
}

void
PointGroup::formBasis(
    const Set<ConstSimpleInternalCoordinatePtr>& simples
)
{
    if (subspaces_.size() != 0) //already done
        return;

    //get rid of redundant coords by passing delocalized
    vector<SymmetryInternalCoordinatePtr> internals;
    VectorPtr evals = InternalCoordinate::addDelocalizedInternalCoordinates(
        mol_,
        simples,
        internals
    );
    formBasis(internals, subspaces_);
        
    //form bases on the subspaces
    for (int i=0; i < subgroups_.size(); ++i)
        subgroups_[i]->formBasis(simples);
}

void
PointGroup::getSubspaces(
    vector<CoordinateSubspacePtr>& subspaces
) const
{
    subspaces = subspaces_;
}

bool PointGroup::initdone_ = false;
map<string, int> PointGroup::symmplanes_;
map<string, int> PointGroup::symmaxes_;
SymmetryOperationPtr SymmetryOperation::identity_op_;
SymmetryOperationPtr SymmetryOperation::c2x_op_;
SymmetryOperationPtr SymmetryOperation::c2y_op_;
SymmetryOperationPtr SymmetryOperation::c2z_op_;
SymmetryOperationPtr SymmetryOperation::inversion_op_;
SymmetryOperationPtr SymmetryOperation::sigmaxy_op_;
SymmetryOperationPtr SymmetryOperation::sigmaxz_op_;
SymmetryOperationPtr SymmetryOperation::sigmayz_op_;
bool SymmetryOperation::initdone_ = false;
