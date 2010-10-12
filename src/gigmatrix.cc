#include <src/gigmatrix.h>
#include <cstdarg>
#include <src/printstream.h>
#include <algorithm>
#include <src/timer.h>
#include <cstring>

#undef heisenbug
#define heisenbug cout << "Heisenbug: " << __FILE__ << " " << __LINE__ << endl

#define VALGRIND 0

extern "C" {

extern void dgemm(const char*, const char*, const int*,
  const int*, const int*, const double*, const double*, const int*,
  const double*, const int*, const double*, double*, const int*);

extern void dgeev(const char*, const char*, const int*, const double*,
  const int*, const double*, const double*, const double*, const int*,
  const double*, const int*, const double*, const int*, const int*);

extern void dsyev(const char*, const char*, const int*, const double*,
  const int*, const double*, const double*, const int*, const int*);

}

using namespace gigide;
using namespace std;

SerialDeclare(Matrix);
SerialDeclare(Vector);

#define X 0
#define Y 1
#define Z 2

void
gigide::error(
    const Matrix* l,
    const Matrix* r,
    const std::string& op
)
{
    cerr << stream_printf("Matrices of size %dx%d and %dx%d not aligned for %s",
                  l->nrow(), l->ncol(),
                  r->nrow(), r->ncol(),
                  op.c_str()) << endl;
    abort();
}

void
gigide::error(
    const Vector* l,
    const Vector* r,
    const std::string& op
)
{
    cerr << stream_printf("Vector of size %d and %d not aligned for %s",
                  l->n(),
                  r->n(),
                  op.c_str()) << endl;
    abort();
}

void
gigide::scale_array(
    double* v,
    double s,
    mdim_t n
)
{
    double* vptr = v;
    for (mdim_t i=0; i < n; ++i, ++vptr)
        (*vptr) *= s;
}

void
gigide::assign_array(
    double* v,
    double d,
    mdim_t n
)
{
    double* vptr = v;
    for (mdim_t i=0; i < n; ++i, ++vptr)
        (*vptr) = d;
}

double 
gigide::max_array(
    const double* v,
    mdim_t n
)
{
    double max = 0;
    for (mdim_t i=0; i < n; ++i)
    {
        double newval = fabs(v[i]);
        if (newval > max)
            max = newval;
    }
    return max;
}

void
gigide::transpose(
    double* vals,
    mdim_t nrow,
    mdim_t ncol
)
{
    double* scratch = new double[nrow*ncol];
    double* scratchptr = scratch;
    double* ptr;
    for (mdim_t i=0; i < ncol; ++i)
    {
        ptr = vals + i;
        for (mdim_t j=0; j < nrow; ++j, ptr += ncol, ++scratchptr)
        {
            *scratchptr = *ptr;
        }
    }
    memcpy(vals, scratch, nrow*ncol*sizeof(double));

    delete[] scratch;
}

void
gigide::eigenvalues(
    const Matrix* matrix,
    Vector* evals,
    Matrix* evecs
)
{

    mdim_t n = matrix->nrow();
    mdim_t nblock = n * n;

    //evecs end up here
    double* m = new double[nblock];
    memcpy(m, matrix->data(), nblock * sizeof(double));
    double* vals = new double[n];
    int worksize = 8*n;
    double* work = new double[worksize];


#if VALGRIND
    double* matptr = m;
    for (mdim_t i=0; i < n; ++i)
    {
        vals[i] = 1;
        for (mdim_t j=0; j < n; ++j, ++matptr)
        {
            if (i==j)
            {
                *matptr = 1;
            }
            else
            {
                *matptr = 0;
            }
        }
    }
#else
    //loss of precision... but what can you do
    int n_ = (int) n;
    int info = 0;
    dsyev("V", "U", &n_, m, &n_, vals, work, &worksize, &info);
    if (info != 0)
    {
        cerr << "Eigenvalue routined failed" << endl;
        matrix->print("matrix", cerr);
        abort();
    }
#endif


    transpose(m, n, n);
    evecs->assign(m);
    evals->assign(vals);

    delete[] work;
    delete[] m;
    delete[] vals;
}

void
gigide::eigenvalues(
    const Matrix* matrix,
    Vector* rvals,
    Vector* ivals,
    Matrix* Levecs,
    Matrix* Revecs
)
{
    if (matrix->nrow() != matrix->ncol())
    {
        cerr << "Attempting to compute eigenvalues of non-square matrix" << endl;
        abort();
    }

    mdim_t n = matrix->nrow();
    mdim_t nblock = n * n;
    //must memcpy since dgeev overwrites array
    double* vals = new double[nblock]; memcpy(vals, matrix->data(), nblock * sizeof(double));
    double* levecs = new double[nblock];
    double* revecs = new double[nblock];

    double* revals = new double[n];
    double* ievals = new double[n];
    int worksize = 8*n;
    double* work = new double[worksize];


#if VALGRIND
    double* revecptr = revecs;
    double* levecptr = levecs;
    double* matptr = vals;
    for (mdim_t i=0; i < n; ++i)
    {
        revals[i] = 1;
        ievals[i] = 0;
        for (mdim_t j=0; j < n; ++j, ++revecptr, ++levecptr, ++matptr)
        {
            if (i==j)
            {
                *revecptr = 1;
                *levecptr = 1;
            }
            else
            {
                *revecptr = 0;
                *levecptr = 0;
            }
        }
    }
#else
    //loss of precision... but what can you do
    int n_ = (int) n;
    int info = 0;
    dgeev("V", "V", &n_, vals, &n_, revals, ievals, levecs, &n_, revecs, &n_, work, &worksize, &info);

    if (info != 0)
    {
        cerr << "Eigenvalue routined failed" << endl;
        matrix->print("matrix", cerr);
        abort();
    }
#endif


    transpose(revecs, n, n);
    transpose(levecs, n, n);

    //these are transposed
    Levecs->assign(revecs);
    Revecs->assign(levecs);
    rvals->assign(revals);
    ivals->assign(ievals);

    delete[] revals;
    delete[] ievals;
    delete[] levecs;
    delete[] revecs;
    delete[] vals;
    delete[] work;
}

double*
gigide::subtract_arrays(
    const double* l,
    const double* r,
    mdim_t n
)
{
    const double* lptr = l;
    const double* rptr = r;
    double* vals = new double[n];
    double* valptr = vals;
    
    for (mdim_t i=0; i < n; ++i, ++lptr, ++rptr, ++valptr)
        (*valptr) = (*lptr) - (*rptr);

    return vals;
}

double*
gigide::add_arrays(
    const double* l,
    const double* r,
    mdim_t n
)
{
    const double* lptr = l;
    const double* rptr = r;
    double* vals = new double[n];
    double* valptr = vals;
    
    for (mdim_t i=0; i < n; ++i, ++lptr, ++rptr, ++valptr)
        (*valptr) = (*lptr) + (*rptr);

    return vals;
}

void
gigide::accumulate_array(
    double* target,
    const double* src,
    mdim_t n
)
{
    double* tptr = target;
    const double* sptr = src;
    for (mdim_t i=0; i < n; ++i, ++tptr, ++sptr)
        (*tptr) += (*sptr);
}

Matrix*
gigide::multiply(
    const Matrix* l,
    const Matrix* r,
    bool transpose_l,
    bool transpose_r
)
{
    const double* ldata = l->data();
    const double* rdata = r->data();

    const char* opl;
    const char* opr;
    mdim_t nrow, ncol, nlink_l, nlink_r;

    //compute l * r
    if (transpose_l)
    {
        opl = "N";
        nrow = l->ncol();
        nlink_l = l->nrow();
    }
    else
    {
        opl = "T";
        nrow = l->nrow();
        nlink_l = l->ncol();
    }

    if (transpose_r)
    {
        opr = "N";
        ncol = r->nrow();
        nlink_r = r->ncol();
    }
    else
    {
        opr = "T";
        ncol = r->ncol();
        nlink_r = r->nrow();
    }

    if (nlink_l != nlink_r)
    {
        stringstream sstr;
        sstr << "multiplication(" << opl << "," << opr << ")";
        error(l, r, sstr.str());
    }

    mdim_t nlink = nlink_r = nlink_l;
    mdim_t ldl = l->ncol();
    mdim_t ldr = r->ncol();

    mdim_t nblock = nrow * ncol;

    double* prod_data = new double[nblock];
    memset(prod_data, 0, nblock * sizeof(double));
    double alpha = 1.0;
    double beta = 0.0;

#if VALGRIND
    //cout << "Multiplying" << endl;
    memset(prod_data, 0, ncol * nrow * sizeof(double));
    double* prodptr = prod_data;
    for (mdim_t i=0; i < nrow; ++i)
    {
        for (mdim_t j=0; j < ncol; ++j, ++prodptr)
        {
            double Xij = 0;
            const double* jptr = rdata + j * ncol;
            const double* iptr = ldata + i * nrow;
            for (mdim_t k=0; k < nlink; ++k, ++iptr, ++jptr)
            {
                Xij += (*iptr) * (*jptr);
                //cout << stream_printf("L[%d][%d] * R[%d][%d] = %12.8f * %12.8f -> %12.8f", i, k, k, j, (*iptr), (*jptr), Xij) << endl;
            }
            *prodptr = Xij;
        }
    }
#else
    int nlink_ = (int) nlink;
    int nrow_ = (int) nrow;
    int ncol_ = (int) ncol;
    int ldl_ = (int) ldl;
    int ldr_ = (int) ldr;
    dgemm(opl, opr, &nrow_, &ncol_, &nlink_, &alpha, ldata, &ldl_, rdata, &ldr_, &beta, prod_data, &nrow_);
#endif

    //this is all transposed
    transpose(prod_data, ncol, nrow);

    Matrix* prod = new Matrix(nrow, ncol);
    prod->assign(prod_data);

    delete[] prod_data;

    return prod;
}

Matrix::Matrix(mdim_t nrow, mdim_t ncol)
    : nrow_(nrow), ncol_(ncol)
{
    SetRuntime(Matrix);

    data_ = new double[nrow * ncol];
    memset(data_, 0, nrow * ncol * sizeof(double));
}

Matrix::Matrix(double* d, mdim_t nrow, mdim_t ncol)
    : nrow_(nrow), ncol_(ncol), data_(d)
{
    SetRuntime(Matrix);
}

Matrix::Matrix(Vector* v, mdim_t nrow, mdim_t ncol)
    : nrow_(nrow), ncol_(ncol)
{
    SetRuntime(Matrix);

    data_ = new double[nrow * ncol];
    assign(v->data());
}

Matrix::Matrix(const XMLArchivePtr& arch)
{
    SetRuntime(Matrix);

    mdim_t size;
    arch->getBinary<double>(data_, size, "data");
    serial_load(nrow);
    serial_load(ncol);

    if (nrow_ * ncol_ != size)
    {
        cerr << "Matrix data not aligned with nrow x ncol" << endl;
        abort();
    }
}

void
Matrix::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    mdim_t size = nrow_ * ncol_;
    arch->setBinary<double>(data_, size, "data");
    serial_save(nrow);
    serial_save(ncol);
}

Matrix::~Matrix()
{
    delete[] data_;
}

Vector*
Matrix::toVector() const
{
    mdim_t nblock = nrow_ * ncol_;
    Vector* v = new Vector(nblock);
    v->assign(data());
    return v;
}

Matrix*
Matrix::add(const Matrix* r) const
{
    if (nrow() != r->nrow() || ncol() != r->ncol())
        error(this, r, "addition");

    const double* ldata = data();
    const double* rdata = r->data();

    double* vals = add_arrays(ldata, rdata, nrow() * ncol());
    Matrix* newptr = new Matrix(nrow(), ncol());
    newptr->assign(vals);

    delete[] vals;

    return newptr;
}


Matrix*
Matrix::subtract(const Matrix* r) const
{
    if (nrow() != r->nrow() || ncol() != r->ncol())
        error(this, r, "subtraction");

    const double* ldata = data();
    const double* rdata = r->data();

    double* vals = subtract_arrays(ldata, rdata, nrow_ * ncol_);
    Matrix* newptr = new Matrix(nrow(), ncol());
    newptr->assign(vals);

    delete[] vals;

    return newptr;
}


Matrix*
Matrix::mult(const Matrix* r) const
{
    Matrix* newptr = multiply(this, r);
    return newptr;
}

Matrix*
Matrix::mult(double d) const
{
    const double* vals = data();
    const double* ptr = vals;
    mdim_t n = nrow() * ncol();

    double* newvals = new double[n];
    double* newptr = newvals;
    for (mdim_t i=0; i < n; ++i, ++ptr, ++newptr)
        (*newptr) = d * (*ptr);
    Matrix* m = new Matrix(nrow(), ncol());
    m->assign(newvals);
    
    delete[] newvals;

    return m;
}

Vector*
Matrix::mult(const Vector* v) const
{
    Matrix* m = v->toMatrix(v->n(), 1);
    Matrix* t = mult(m);
    Vector* newv = t->toVector();
    delete m;
    delete t;
    return newv;
}

Matrix*
Matrix::t() const 
{
    mdim_t nblock = nrow_ * ncol_;
    double* vals = new double[nblock];
    memcpy(vals, data_, nblock * sizeof(double));
    transpose(vals, nrow_, ncol_);
    Matrix* m = new Matrix(vals, ncol_, nrow_);
    return m;
}

const double*
Matrix::data() const
{
    return data_;
}

void
Matrix::accumulate(const Matrix* r)
{
    accumulate_array(data_, r->data(), nrow_ * ncol_);
}

mdim_t
Matrix::nrow() const
{
    return nrow_;
}

mdim_t
Matrix::ncol() const
{
    return ncol_;
}

void
Matrix::scale(double s)
{
    scale_array(data_, s, nrow_ * ncol_);
}

void
Matrix::assign(const double* newvals)
{
    memcpy(data_, newvals, ncol_ * nrow_ * sizeof(double));
}

void
Matrix::assign(double d)
{
    assign_array(data_, d, ncol_ * nrow_);
}

void
Matrix::set_element(mdim_t i, mdim_t j, double val)
{
    mdim_t index = i * ncol_ + j;
    data_[index] = val;
}

double
Matrix::get_element(mdim_t i, mdim_t j) const
{
    mdim_t index = i * ncol_ + j;
    return data_[index];
}

void
Matrix::accumulate_element(mdim_t i, mdim_t j, double val)
{
    mdim_t index = i * ncol_ + j;
    data_[index] += val;
}

void
Matrix::printValues(std::ostream& os) const
{
    for (mdim_t row=0; row < nrow_; ++row)
    {
        os << "  ";
        for (mdim_t col=0; col < ncol_; ++col)
        {
            os << stream_printf("  %16.10f", get_element(row, col));
        }
        os << endl;
    }
}

void
Matrix::print(const std::string& title, std::ostream& os) const
{
    os << stream_printf("%d x %d Matrix %s", nrow_, ncol_, title.data()) << endl;
    printValues(os);
}

Vector::Vector(mdim_t n) :
    n_(n)
{
    SetRuntime(Vector);

    data_ = new double[n];
    memset(data_, 0, n * sizeof(double));
}

Vector::Vector(double* d, mdim_t n) :
    n_(n), data_(d)
{
    SetRuntime(Vector);
}

Vector::Vector(const XMLArchivePtr& arch)
{
    SetRuntime(Vector);

    mdim_t size;
    arch->getBinary(data_, size, "data");
    serial_load(n);

    if (size != n_)
    {
        cerr << stream_printf("Data size %ld does not match dim %ld", size, n_) << endl;
        abort();
    }
}

void
Vector::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);

    arch->setBinary(data_, n_, "data");
    serial_save(n);
}

Matrix*
Vector::toMatrix(mdim_t nrow, mdim_t ncol) const
{
    if (nrow * ncol != n())
    {
        cerr << stream_printf("Vector of size %d cannot be converted to matrix of size %d x %d",
                              n(), nrow, ncol) << endl;
        abort();
    }
        
    Matrix* m = new Matrix(nrow, ncol);
    m->assign(data());
    return m;
}

Vector::~Vector()
{
    delete[] data_;
}

Vector*
Vector::add(const Vector* v) const
{
    if (n() != v->n())
        error(this, v, "addition");

    double* vals = add_arrays(data_, v->data(), n());

    Vector* newptr = new Vector(n());
    newptr->assign(vals);

    delete[] vals;

    return newptr;
}

Vector*
Vector::subtract(const Vector* v) const
{
    if (n() != v->n())
        error(this, v, "subtraction");

    double* newvals = subtract_arrays(data(), v->data(), n()); 
    Vector* vec = new Vector(n());
    vec->assign(newvals);

    delete[] newvals;

    return vec;
}

Vector*
Vector::mult(double d) const
{
    double* vals = new double[n()];
    for (mdim_t i=0; i < n(); ++i)
        vals[i] = data_[i] * d;

    Vector* v = new Vector(n());
    v->assign(vals);

    delete[] vals;

    return v;
}

const double*
Vector::data() const
{
    return data_;
}

void
Vector::accumulate(const Vector* v)
{
    accumulate_array(data_, v->data(), n());
}

void
Vector::scale(double s)
{
    scale_array(data_, s, n_);
}

void
Vector::assign(const double* vals)
{
    memcpy(data_, vals, n() * sizeof(double));
}

void
Vector::assign(double d)
{
    assign_array(data_, d, n_);
}

double
Vector::dot(const Vector* v) const
{
    if (n() != v->n())
        error(this, v, "dot product");

    const double* ldata = data();
    const double* rdata = v->data();
    const double* lptr = ldata;
    const double* rptr = rdata;
    
    double sum = 0;
    for (mdim_t i=0; i < n(); ++i, ++lptr, ++rptr)
        sum += (*lptr) * (*rptr);
        
    return sum;
}

double
Vector::get_element(mdim_t i) const
{
    return data_[i];
}

void
Vector::set_element(mdim_t i, double val)
{
    data_[i] = val;
}

void
Vector::accumulate_element(mdim_t i, double val)
{
    data_[i] += val;
}

void
Vector::sort()
{
    std::sort(data_, data_ + n_);
}

void
Vector::printValues(std::ostream& os) const
{
    for (mdim_t i=0; i < n_; ++i)
    {
        os << stream_printf("  %16.10f", get_element(i));
        os << endl;
    }
}

void
Vector::print(const std::string& title, std::ostream& os) const
{
    os << stream_printf("Vector %s of size %d", title.data(), n_) << endl;
    printValues(os);
}

mdim_t
Vector::n() const
{
    return n_;
}

Vector1::Vector1(double a0)
    : VectorPtr(1)
{
    set_element(0, a0);
}

Vector2::Vector2(double a0, double a1)
    : VectorPtr(2)
{
    set_element(0, a0);
    set_element(1, a1);
}

Vector3::Vector3(double a0, double a1, double a2)
    : VectorPtr(3)
{
    set_element(0, a0);
    set_element(1, a1);
    set_element(2, a2);
}

Vector4::Vector4(double a0, double a1, double a2, double a3)
    : VectorPtr(4)
{
    set_element(0, a0);
    set_element(1, a1);
    set_element(2, a2);
    set_element(3, a3);
}

Vector8::Vector8(double a0, double a1, double a2, double a3,
                 double a4, double a5, double a6, double a7)
    : VectorPtr(8)
{
    set_element(0, a0);
    set_element(1, a1);
    set_element(2, a2);
    set_element(3, a3);
    set_element(4, a4);
    set_element(5, a5);
    set_element(6, a6);
    set_element(7, a7);
}


VectorPtr
gigide::cross(const ConstVectorPtr& v, const ConstVectorPtr& w)
{
    v.nullcheck();
    w.nullcheck();

    VectorPtr cp = v.clone();
    double newx = v[Y]*w[Z] - v[Z]*w[Y];
    double newy = v[Z]*w[X] - v[X]*w[Z];
    double newz = v[X]*w[Y] - v[Y]*w[X];
    cp.set_element(X, newx);
    cp.set_element(Y, newy);
    cp.set_element(Z, newz);
    return cp;
}

bool
gigide::equals(const ConstRectMatrixPtr& l, const ConstRectMatrixPtr& r, double tol)
{
    l.nullcheck();
    r.nullcheck();

    RectMatrixPtr diff = l - r;
    return diff.maxabs() < tol;
}

bool
gigide::equals(const ConstVectorPtr& l, const ConstVectorPtr& r, double tol)
{
    l.nullcheck();
    r.nullcheck();

    VectorPtr diff = l - r;
    return diff.maxabs() < tol;
}

RectMatrixPtr gigide::operator*(const ConstRectMatrixPtr& l, const ConstRectMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Matrix* m = l.get()->mult(r.get());    

    return m;
}

RectMatrixPtr gigide::operator*(const ConstRectMatrixPtr& l, const ConstSymmMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();
    Matrix* m = l.get()->mult(r.get());    

    return m;
}

RectMatrixPtr gigide::operator*(const ConstSymmMatrixPtr& l, const ConstRectMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();
    Matrix* m = l.get()->mult(r.get());    

    return m;
}

RectMatrixPtr gigide::operator*(const ConstSymmMatrixPtr& l, const ConstSymmMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();
    Matrix* m = l.get()->mult(r.get());    

    return m;
}

RectMatrixPtr gigide::operator*(const ConstRectMatrixPtr& m, double d)
{
    m.nullcheck();

    Matrix* p = m.get()->mult(d);

    return p;
}

RectMatrixPtr gigide::operator*(double d, const ConstRectMatrixPtr& m)
{
    m.nullcheck();

    Matrix* p = m.get()->mult(d);

    return p;
}

RectMatrixPtr gigide::operator+(const ConstRectMatrixPtr& l, const ConstRectMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Matrix* m = l.get()->subtract(r.get());

    return m;
}

RectMatrixPtr gigide::operator-(const ConstRectMatrixPtr& l, const ConstRectMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Matrix* m = l.get()->subtract(r.get());

    return m;
}

SymmMatrixPtr gigide::operator+(const ConstSymmMatrixPtr& l, const ConstSymmMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Matrix* m = l.get()->add(r.get());

    return m;
}

SymmMatrixPtr gigide::operator-(const ConstSymmMatrixPtr& l, const ConstSymmMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Matrix* m = l.get()->subtract(r.get());

    return m;
}

SymmMatrixPtr gigide::operator*(const ConstSymmMatrixPtr& m, double d)
{
    m.nullcheck();

    Matrix* p = m.get()->mult(d);

    return p;
}

SymmMatrixPtr gigide::operator*(double d, const ConstSymmMatrixPtr& m)
{
    m.nullcheck();

    Matrix* p = m.get()->mult(d);

    return p;
}

VectorPtr gigide::operator*(const ConstRectMatrixPtr& m, const ConstVectorPtr& v)
{
    m.nullcheck();
    v.nullcheck();
    Vector* p =  m.get()->mult(v.get());

    return p;
}

VectorPtr gigide::operator*(const ConstSymmMatrixPtr& m, const ConstVectorPtr& v)
{
    m.nullcheck();
    v.nullcheck();
    Vector* p =  m.get()->mult(v.get());

    return p;
}

VectorPtr gigide::operator*(double d, const ConstVectorPtr& v)
{
    v.nullcheck();
    Vector* p = v.get()->mult(d);

    return p;
}

VectorPtr gigide::operator*(const ConstVectorPtr& v, double d)
{
    v.nullcheck();

    Vector* p = v.get()->mult(d);

    return p;
}

VectorPtr gigide::operator+(const ConstVectorPtr& l, const ConstVectorPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Vector* v = l.get()->add(r.get());

    return v;
}

VectorPtr gigide::operator-(const ConstVectorPtr& l, const ConstVectorPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Vector* v = l.get()->subtract(r.get());

    return v;
}

SymmMatrixPtr::SymmMatrixPtr()
    : Parent()
{
}

SymmMatrixPtr::SymmMatrixPtr(Matrix* m)
    : Parent(m)
{
}

SymmMatrixPtr::SymmMatrixPtr(mdim_t n)
    : Parent(n)
{
}


SymmMatrixPtr::SymmMatrixPtr(const SymmMatrixPtr& m)
    : Parent(m)
{
}

ConstSymmMatrixPtr::ConstSymmMatrixPtr(const Matrix* m)
    : Parent(m)
{
}

ConstSymmMatrixPtr::ConstSymmMatrixPtr(const SymmMatrixPtr& m)
    : Parent(m.get())
{
}

ConstSymmMatrixPtr::ConstSymmMatrixPtr(const ConstSymmMatrixPtr& m)
    : Parent(m.get())
{
}

RectMatrixPtr::RectMatrixPtr()
    : Parent()
{
}

RectMatrixPtr::RectMatrixPtr(Matrix* m)
    : Parent(m)
{
}

RectMatrixPtr::RectMatrixPtr(mdim_t nrow, mdim_t ncol)
    : Parent(nrow, ncol)
{
}

RectMatrixPtr::RectMatrixPtr(const RectMatrixPtr& m)
    : Parent(m)
{
}

ConstRectMatrixPtr::ConstRectMatrixPtr()
    : Parent()
{
}

ConstRectMatrixPtr::ConstRectMatrixPtr(const Matrix* m)
    : Parent(m)
{
}

ConstRectMatrixPtr::ConstRectMatrixPtr(const RectMatrixPtr& m)
    : Parent(m.get())
{
}

ConstRectMatrixPtr::ConstRectMatrixPtr(const ConstRectMatrixPtr& m)
    : Parent(m.get())
{
}

VectorPtr::VectorPtr()
    : Parent()
{
}

VectorPtr::VectorPtr(mdim_t n)
    : Parent(n)
{
}

VectorPtr::VectorPtr(Vector* v)
    : Parent(v)
{
}

VectorPtr::VectorPtr(const VectorPtr& v)
    : Parent(v.get())
{
}

ConstVectorPtr::ConstVectorPtr()
    : Parent()
{
}

ConstVectorPtr::ConstVectorPtr(mdim_t n)
    : Parent(n)
{
}

ConstVectorPtr::ConstVectorPtr(const Vector* v)
    : Parent(v)
{
}

ConstVectorPtr::ConstVectorPtr(const VectorPtr& v)
    : Parent(v.get())
{
}

ConstVectorPtr::ConstVectorPtr(const ConstVectorPtr& v)
    : Parent(v.get())
{
}

void
gigide::instantiate()
{
    RectMatrixPtr m;
    VectorPtr v(3);
    SymmMatrixPtr s;
    const double* vals;

    RectMatrixTemplate<Matrix>* mptr = new RectMatrixTemplate<Matrix>;
    mptr->print("test");
    delete mptr;

    m.t();
    m.get_subblock(0,0,0,0);
    m.copy();
    m.clone();
    m.set_element(0,0,0);
    m.accumulate_element(0,0,0);
    m.assign_subblock(m,0,0);
    m.assign(vals);
    m.assign(m);
    m.assign_row(v,0);
    m.assign_column(v,0);
    m.accumulate(m);
    m.eigen(v,v,m,m);
    m.print("test");
    m.printValues();
    m.scale(0.0);
    m.toVector();
    m.nonnull();
    m.null();
    m.get_row(0);
    m.get_column(0);
    m.zero();
    m.nrow();
    m.ncol();
    vals = m.data();

    ConstRectMatrixPtr cm(m);
    cm.t();
    cm.get_subblock(0,0,0,0);
    cm.copy();
    cm.clone();
    cm.eigen(v,v,m,m);
    cm.print("test");
    cm.printValues();
    cm.toVector();
    cm.nonnull();
    cm.null();
    cm.get_row(0);
    cm.get_column(0);
    m.nrow();
    m.ncol();
    vals = m.data();


    s.i();
    s.sqrt_matrix();
    s.invsqrt_matrix();
    s.n();
    s.set_element(0,0,0);
    s.accumulate_element(0,0,0);
    s.assign(s);
    s.eigen(v,m);
    s.accumulate_symmetric_product(m);
    s.accumulate_transform(m,s);
    s.accumulate_transform(s,s);
    s.clone();
    s.copy();

    ConstSymmMatrixPtr cs(s);
    cs.i();
    cs.sqrt_matrix();
    cs.invsqrt_matrix();
    cs.n();
    cs.eigen(v,m);
    cs.clone();
    cs.copy();

    v[0];
    v.dot(v);
    v.get_element(0);
    v.set_element(0,0);
    v.accumulate_element(0,0);
    v.assign(vals);
    v.assign(v);
    v.printValues();
    v.print("test");
    v.null();
    v.nonnull();
    v.n();
    v.normalize();
    v.sort();
    v.accumulate(v);
    v.scale(0.5);
    v.assign(0.0);
    v.clone();
    v.copy();
    v.zero();
    v.maxabs();
    v.norm();
    v.symmetric_outer_product();
    v.data();
    v.nullcheck();
    v.sort();
    
    ConstVectorPtr cv(v);
    cv[0];
    cv.dot(v);
    cv.get_element(0);
    cv.printValues();
    cv.print("test");
    cv.null();
    cv.nonnull();
    cv.n();
    cv.clone();
    cv.copy();
    cv.maxabs();
    cv.norm();
    cv.symmetric_outer_product();
    cv.data();
    cv.nullcheck();

}


