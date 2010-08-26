#include <src/gigmatrix.h>
#include <cstdarg>
#include <math.h>
#include <src/printstream.h>
#include <algorithm>

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
error(
    const Matrix* l,
    const Matrix* r,
    string op
)
{
    cerr << stream_printf("Matrices of size %dx%d and %dx%d not aligned for %s",
                  l->nrow(), l->ncol(),
                  r->nrow(), r->ncol(),
                  op.c_str()) << endl;
    abort();
}

void
error(
    const Vector* l,
    const Vector* r,
    string op
)
{
    cerr << stream_printf("Vector of size %d and %d not aligned for %s",
                  l->n(),
                  r->n(),
                  op.c_str()) << endl;
    abort();
}

void
scale_array(
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
assign_array(
    double* v,
    double d,
    mdim_t n
)
{
    double* vptr = v;
    for (mdim_t i=0; i < n; ++i, ++vptr)
        (*vptr) = d;
}

double max_array(
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
transpose(
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
eigenvalues(
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
eigenvalues(
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
subtract_arrays(
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
add_arrays(
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
accumulate_array(
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
multiply(
    const Matrix* l,
    const Matrix* r,
    bool transpose_l = false,
    bool transpose_r = false
);

Matrix*
multiply(
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

    if (nrow < 0 || ncol < 0)
        abort();

    data_ = new double[nrow * ncol];
    memset(data_, 0, nrow * ncol * sizeof(double));
}

Matrix::Matrix(double* d, mdim_t nrow, mdim_t ncol)
    : nrow_(nrow), ncol_(ncol), data_(d)
{
    SetRuntime(Matrix);

    if (nrow < 0 || ncol < 0)
        abort();
}

Matrix::Matrix(Vector* v, mdim_t nrow, mdim_t ncol)
    : nrow_(nrow), ncol_(ncol)
{
    SetRuntime(Matrix);

    if (nrow < 0 || ncol < 0)
        abort();

    data_ = new double[nrow * ncol];
    assign(v->data());
}

Matrix::Matrix(const ArchivePtr& arch)
{
    SetRuntime(Matrix);

    mdim_t size;
    void* ptr;
    ConstArchivePtr node = arch->loadBinary(&ptr, size, "data");
    data_ = (double *) ptr;
    node->getAttribute(nrow_, "nrow");
    node->getAttribute(ncol_, "ncol");

    if (node.get() == NULL)
    {
        cerr << "No data on matrix node" << endl;
        abort();
    }

    if (size != nrow_ * ncol_ * sizeof(double))
    {
        cerr << stream_printf("Data size %ld does not match dims %ld x %ld", size, nrow_, ncol_) << endl;
        abort();
    }
}

void
Matrix::serialize(const ArchivePtr& arch) const
{
    Serializable::serialize(arch);
    mdim_t size = nrow_ * ncol_ * sizeof(double);
    ArchivePtr node = arch->writeBinary(data_, size, "data");
    node->setAttribute(nrow_, "nrow");
    node->setAttribute(ncol_, "ncol");
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
Matrix::print(std::string title, std::ostream& os) const
{
    os << stream_printf("%d x %d Matrix %s", nrow_, ncol_, title.c_str()) << endl;
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

    //data_ = new double[n];
    //memcpy(data_, d, n * sizeof(double));
}

Vector::Vector(const ConstArchivePtr& arch)
{
    SetRuntime(Vector);

    mdim_t size;
    void* ptr;
    ConstArchivePtr node = arch->loadBinary(&ptr, size, "data");
    data_ = (double *) ptr;
    node->getAttribute(n_, "n");

    if (node.get() == NULL)
    {
        cerr << "No data on matrix node" << endl;
        abort();
    }

    if (size != n_ * sizeof(double))
    {
        cerr << stream_printf("Data size %ld does not match dim %ld", size, n_) << endl;
        abort();
    }
}

void
Vector::serialize(const ArchivePtr& arch) const
{
    Serializable::serialize(arch);

    ArchivePtr node = arch->writeBinary(data_, n_ * sizeof(double), "data");
    node->setAttribute(n_, "n");

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
Vector::print(std::string title, std::ostream& os) const
{
    os << stream_printf("Vector %s of size %d", title.c_str(), n_) << endl;
    printValues(os);
}

mdim_t
Vector::n() const
{
    return n_;
}

template <class T> 
MatrixTemplate<T>::MatrixTemplate()
    : Parent(NULL)
{
}

template <class T> 
MatrixTemplate<T>::MatrixTemplate(T* m)
    : Parent(m)
{
}

template <class T> 
bool
MatrixTemplate<T>::null() const
{
    return get() == NULL;
}

template <class T> 
bool
MatrixTemplate<T>::nonnull() const
{
    return get() != NULL;
}

template <class T> 
void
MatrixTemplate<T>::scale(double s)
{
    get()->scale(s);
}

template <class T> 
void
MatrixTemplate<T>::zero()
{
    get()->assign(0.0);
}

template <class T> 
double
MatrixTemplate<T>::get_element(mdim_t i, mdim_t j) const
{
    return get()->get_element(i,j);
}

template <class T> 
const VectorPtr
MatrixTemplate<T>::toConstVector() const
{
    mdim_t size = get()->nrow() * get()->ncol();
    Vector* v = new Vector(const_cast<double*>(data()), size);
    const VectorPtr vptr(v);
    return vptr;
}

template <class T> 
VectorPtr
MatrixTemplate<T>::toVector() const
{
    Vector* v = get()->toVector();
    return v;
}

template <class T> 
const double*
MatrixTemplate<T>::data() const
{
    return get()->data();
}

template <class T> 
void
MatrixTemplate<T>::printValues(ostream& os) const
{
    nullcheck();
    get()->printValues(os);
}

template <class T> 
void
MatrixTemplate<T>::print(string title, ostream& os) const
{
    nullcheck();
    get()->print(title, os);
}

template <class T> 
void
MatrixTemplate<T>::nullcheck() const
{
    if (null())
    {
        cerr << "Called method on null matrix" << endl;
        abort();
    }
}

template <class T> 
double
MatrixTemplate<T>::maxabs() const
{
    nullcheck();
    mdim_t nrow = get()->nrow();
    mdim_t ncol = get()->ncol();
    return max_array(data(), nrow * ncol);
}

template <class T> 
VectorPtr
MatrixTemplate<T>::get_row(mdim_t row) const
{
    nullcheck();
    mdim_t ncol = get()->ncol();
    Vector* v = new Vector(ncol);
    for (mdim_t col=0; col < ncol; ++col)
        v->set_element(col, get_element(row, col));
    return v;
}

template <class T> 
VectorPtr
MatrixTemplate<T>::get_column(mdim_t col) const
{
    nullcheck();
    mdim_t nrow = get()->nrow();
    Vector* v = new Vector(nrow);
    for (mdim_t row=0; row < nrow; ++row)
        v->set_element(row, get_element(row, col));
    return v;
}


template <class T> 
mdim_t
MatrixTemplate<T>::nrow() const
{
    nullcheck();
    return get()->nrow();
}


template <class T> 
mdim_t
MatrixTemplate<T>::ncol() const
{
    nullcheck();
    return get()->ncol();
}

template <class T> 
RectMatrixTemplate<T>::RectMatrixTemplate()
    : Parent(NULL)
{
}

template <class T> 
RectMatrixTemplate<T>::~RectMatrixTemplate()
{
}

template <class T> 
RectMatrixTemplate<T>::RectMatrixTemplate(T* m)
    : Parent(m)
{
}

template <class T> 
RectMatrixTemplate<T>::RectMatrixTemplate(mdim_t nrow, mdim_t ncol)
    : Parent(new Matrix(nrow, ncol))
{
}

template <class T> 
RectMatrixTemplate<T>::RectMatrixTemplate(const RectMatrixPtr& m)
    : Parent(m.get())
{
}

template <class T> 
void
RectMatrixTemplate<T>::accumulate_element(mdim_t i, mdim_t j, double val)
{
    nullcheck();
    get()->accumulate_element(i,j,val);
}

template <class T>
void
RectMatrixTemplate<T>::set_element(mdim_t i, mdim_t j, double val)
{
    nullcheck();
    get()->set_element(i,j,val);
}

template <class T> 
void
RectMatrixTemplate<T>::accumulate(const ConstRectMatrixPtr& m)
{
    nullcheck();
    get()->accumulate(m.get());
}

template <class T> 
void
RectMatrixTemplate<T>::assign(const ConstRectMatrixPtr& m)
{
    nullcheck();
    get()->assign(m.data());
}

template <class T> 
void
RectMatrixTemplate<T>::assign(const double* vals)
{
    nullcheck();
    get()->assign(vals);
}


template <class T> 
RectMatrixPtr
RectMatrixTemplate<T>::t() const
{
    nullcheck();
    return get()->t();
}

template <class T> 
void
RectMatrixTemplate<T>::assign_row(const ConstVectorPtr& v, mdim_t row)
{
    nullcheck();
    for (mdim_t col=0; col < Parent::ncol(); ++col)
        set_element(row, col, v[col]);
}

template <class T> 
void
RectMatrixTemplate<T>::assign_column(const ConstVectorPtr& v, mdim_t col)
{
    nullcheck();
    for (mdim_t row=0; row < Parent::nrow(); ++row)
        set_element(row, col, v[row]);
}

template <class T> 
RectMatrixPtr
RectMatrixTemplate<T>::clone() const
{
    nullcheck();
    return new Matrix(Parent::nrow(), Parent::ncol());
}

template <class T> 
RectMatrixPtr
RectMatrixTemplate<T>::copy() const
{
    nullcheck();
    Matrix* m = new Matrix(Parent::nrow(), Parent::ncol());
    m->assign(Parent::data());
    return m;
}

template <class T> 
RectMatrixPtr
RectMatrixTemplate<T>::get_subblock(mdim_t rowstart, mdim_t rowstop, mdim_t colstart, mdim_t colstop) const
{
    nullcheck();
    mdim_t nr = rowstop - rowstart + 1;
    mdim_t nc = colstop - colstart + 1;
    Matrix* ptr = new Matrix(nr, nc);
    mdim_t r = 0;
    for (mdim_t row=rowstart; row <= rowstop; ++row, ++r)
    {
        mdim_t c = 0;
        for (mdim_t col=colstart; col <= colstop; ++col, ++c)
        {
            ptr->set_element(r, c, Parent::get_element(row, col)); 
        }
    }
    return ptr;
}

template <class T> 
void
RectMatrixTemplate<T>::assign_subblock(const ConstRectMatrixPtr& block, mdim_t rowstart, mdim_t colstart)
{
    nullcheck();
    mdim_t nr = block.nrow();
    mdim_t nc = block.ncol();
    mdim_t row, r, col, c;
    for (r=0, row=rowstart; r < nr; ++r, ++row)
    {
        for (c=0, col=colstart; c < nc; ++c, ++col)
        {
            set_element(row, col, block.get_element(r,c));
        }
    }
}

template <class T> 
void
RectMatrixTemplate<T>::eigen(VectorPtr& revals, VectorPtr& ievals, RectMatrixPtr& Levecs, RectMatrixPtr& Revecs) const
{
    nullcheck();

    if (revals.null()) revals = new Vector(Parent::nrow());
    if (ievals.null()) ievals = new Vector(Parent::nrow());
    if (Levecs.null()) Levecs = new Matrix(Parent::nrow(), Parent::nrow());
    if (Revecs.null()) Revecs = new Matrix(Parent::nrow(), Parent::nrow());

    eigenvalues(get(), revals.get(), ievals.get(), Levecs.get(), Revecs.get());
}

template <class T>
VectorTemplate<T>::VectorTemplate(mdim_t n)
    : Parent(new Vector(n))
{
}

template <class T>
VectorTemplate<T>::VectorTemplate()
    : Parent(NULL)
{
}

template <class T>
VectorTemplate<T>::VectorTemplate(T* v)
    : Parent(v)
{
}

template <class T>
VectorTemplate<T>::VectorTemplate(const ConstVectorPtr& v)
    : Parent(v.get())
{
}

template <class T>
mdim_t
VectorTemplate<T>::n() const
{
    nullcheck();
    return get()->n();
}

template <class T>
double
VectorTemplate<T>::dot(const ConstVectorPtr& v) const
{
    nullcheck();
    return get()->dot(v.get());
}

template <class T>
void
VectorTemplate<T>::nullcheck() const
{
    if (null())
    {
        cerr << "Called method on null matrix" << endl;
        abort();
    }
}

template <class T>
void
VectorTemplate<T>::printValues(ostream& os) const
{
    nullcheck();
    get()->printValues(os);
}

template <class T>
void
VectorTemplate<T>::print(string title, ostream& os) const
{
    nullcheck();
    get()->print(title, os);
}

template <class T>
void
VectorTemplate<T>::assign(const ConstVectorPtr& v)
{
    nullcheck();
    get()->assign(v.data());
}

template <class T>
void
VectorTemplate<T>::assign(const double* vals)
{
    nullcheck();
    get()->assign(vals);
}

template <class T>
const double*
VectorTemplate<T>::data() const
{
    nullcheck();
    return get()->data();
}

template <class T>
double
VectorTemplate<T>::operator[](mdim_t i) const
{
    nullcheck();
    return get()->get_element(i);
}

template <class T>
double
VectorTemplate<T>::get_element(mdim_t i) const
{
    nullcheck();
    return get()->get_element(i);
}

template <class T>
void
VectorTemplate<T>::set_element(mdim_t i, double val)
{
    nullcheck();
    get()->set_element(i, val);
}

template <class T>
void
VectorTemplate<T>::accumulate_element(mdim_t i, double val)
{
    nullcheck();
    get()->accumulate_element(i, val);
}

template <class T>
void
VectorTemplate<T>::sort()
{
    get()->sort();
}

template <class T>
bool
VectorTemplate<T>::null() const
{
    return get() == NULL;
}

template <class T>
bool
VectorTemplate<T>::nonnull() const
{
    return get() != NULL;
}

template <class T>
void
VectorTemplate<T>::zero()
{
    nullcheck();
    get()->assign(0.0);
}

template <class T>
double
VectorTemplate<T>::norm() const
{
    nullcheck();
    double n2 = get()->dot(get());
    return sqrt(n2);
}

template <class T>
void
VectorTemplate<T>::accumulate(const ConstVectorPtr& v)
{
    nullcheck();
    v.nullcheck();
    get()->accumulate(v.get());
}

template <class T>
void
VectorTemplate<T>::scale(double d)
{
    nullcheck();
    get()->scale(d);
}

template <class T>
void
VectorTemplate<T>::assign(double d)
{
    nullcheck();
    get()->assign(d);
}

template <class T>
void
VectorTemplate<T>::normalize()
{
    nullcheck();
    double normsq = get()->dot(get()); 
    double oonorm = 1.0/sqrt(normsq);
    scale(oonorm);
}

template <class T>
VectorPtr
VectorTemplate<T>::copy() const
{
    nullcheck();
    Vector* v = new Vector(n());
    v->assign(data());
    return v;
}

template <class T>
VectorPtr
VectorTemplate<T>::clone() const
{
    nullcheck();
    Vector* v = new Vector(n());
    return v;
}

template <class T>
double
VectorTemplate<T>::maxabs() const
{
    nullcheck();
    return max_array(data(), n());
}

template <class T>
SymmMatrixPtr
VectorTemplate<T>::symmetric_outer_product() const
{
    nullcheck();
    const double* v = data();
    Matrix* m = new Matrix(n(), n());
    for (mdim_t i=0; i < n(); ++i)
    {
        for (mdim_t j=0; j < n(); ++j)
        {
            m->set_element(i, j, v[i] * v[j]);
        }
    }
    return m;
}

template <class T> 
SymmMatrixTemplate<T>::SymmMatrixTemplate()
    : Parent(NULL)
{
}

template <class T> 
SymmMatrixTemplate<T>::SymmMatrixTemplate(mdim_t n)
    : Parent(new Matrix(n,n))
{
}

template <class T> 
SymmMatrixTemplate<T>::SymmMatrixTemplate(T* m)
    : Parent(m)
{
}

template <class T> 
SymmMatrixTemplate<T>::SymmMatrixTemplate(const SymmMatrixPtr& m)
    : Parent(m.get())
{
    m.nullcheck();
}

template <class T> 
mdim_t
SymmMatrixTemplate<T>::n() const
{
    nullcheck();
    return get()->nrow();
}

template <class T> 
void
SymmMatrixTemplate<T>::set_element(mdim_t i, mdim_t j, double val)
{
    nullcheck();
    get()->set_element(i,j,val);
    get()->set_element(j,i,val);
}

template <class T> 
void
SymmMatrixTemplate<T>::accumulate_element(mdim_t i, mdim_t j, double val)
{
    nullcheck();
    get()->accumulate_element(i, j, val);
}

template <class T> 
void
SymmMatrixTemplate<T>::accumulate(const SymmMatrixPtr& m)
{
    nullcheck();
    get()->accumulate(m.get());
}

template <class T> 
SymmMatrixPtr
SymmMatrixTemplate<T>::clone() const
{
    nullcheck();
    return new Matrix(n(), n());
}

template <class T> 
SymmMatrixPtr
SymmMatrixTemplate<T>::copy() const
{
    nullcheck();
    Matrix* m = new Matrix(n(), n());
    m->assign(Parent::data());
    return m;
}

template <class T> 
void
SymmMatrixTemplate<T>::assign(const SymmMatrixPtr& m)
{
    get()->assign(m->data());
}

template <class T> 
void
SymmMatrixTemplate<T>::eigen(VectorPtr& evals, RectMatrixPtr& evecs) const
{
    nullcheck();

    if (evals.null()) evals = new Vector(n());
    if (evecs.null()) evecs = new Matrix(n(), n());

    eigenvalues(get(), evals.get(), evecs.get());
}

template <class T> 
SymmMatrixPtr
SymmMatrixTemplate<T>::sqrt_matrix(double tol) const
{
    nullcheck();
    RectMatrixPtr evecs;
    VectorPtr evals;
    SymmMatrixPtr epsilon(n());

    eigen(evals, evecs);
    for (mdim_t i=0; i < n(); ++i)
    {
        double val = evals[i];
        if ( fabs(val) > tol)
            epsilon.set_element(i, i, sqrt(val));
    }

    SymmMatrixPtr sqrtmat = clone();
    sqrtmat.accumulate_transform(evecs, epsilon);

    return sqrtmat;
}

template <class T> 
SymmMatrixPtr
SymmMatrixTemplate<T>::invsqrt_matrix(double tol) const
{
    nullcheck();
    RectMatrixPtr evecs;
    VectorPtr evals;
    SymmMatrixPtr epsilon(n());

    eigen(evals, evecs);
    for (mdim_t i=0; i < n(); ++i)
    {
        double val = evals[i];
        if ( fabs(val) > tol)
            epsilon.set_element(i, i, 1.0/sqrt(val));
    }

    SymmMatrixPtr sqrtmat = clone();
    sqrtmat.accumulate_transform(evecs, epsilon);

    return sqrtmat;
}

template <class T> 
SymmMatrixPtr
SymmMatrixTemplate<T>::i(double tol) const
{
    nullcheck();
    RectMatrixPtr evecs;
    VectorPtr evals;
    SymmMatrixPtr epsilon(n());

    eigen(evals, evecs);
    for (mdim_t i=0; i < n(); ++i)
    {
        double val = evals[i];
        if ( fabs(val) > tol )
            epsilon.set_element(i, i, 1.0/val);
    }

    SymmMatrixPtr inv = clone();
    inv.accumulate_transform(evecs, epsilon);

    return inv;
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

template <class T> 
void
SymmMatrixTemplate<T>::accumulate_transform(const ConstRectMatrixPtr& r, const ConstSymmMatrixPtr& s)
{
    nullcheck();
    r.nullcheck();
    s.nullcheck();

    Matrix* temp = multiply(r.get(), s.get());
    Matrix* final = multiply(temp, r.get(), false, true);
    get()->accumulate(final);
    delete temp;
    delete final;
}

template <class T> 
void
SymmMatrixTemplate<T>::accumulate_transform(const ConstSymmMatrixPtr& r, const ConstSymmMatrixPtr& s)
{
    nullcheck();
    r.nullcheck();
    s.nullcheck();
    Matrix* temp = multiply(r.get(), s.get());

    Matrix* final = multiply(temp, r.get());
    get()->accumulate(final);
    delete temp;
    delete final;
}

template <class T> 
void
SymmMatrixTemplate<T>::accumulate_symmetric_product(const ConstRectMatrixPtr& m)
{
    nullcheck();
    m.nullcheck();

    Matrix* prod = multiply(m.get(), m.get(), false, true); //multiply by transpose
    get()->accumulate(prod);
    delete prod;
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
instantiate()
{
    RectMatrixPtr m;
    VectorPtr v;
    SymmMatrixPtr s;
    const double* vals;

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


