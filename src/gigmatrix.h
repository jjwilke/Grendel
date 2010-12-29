#ifndef gigide_matrix_h
#define gigide_matrix_h

#include <src/smartptr/src/serialize.h>
#include <iostream>
#include <math.h>

namespace gigide {

class Vector;
class VectorPtr;
class ConstVectorPtr;
class RectMatrixPtr;
class ConstRectMatrixPtr;
class SymmMatrixPtr;
class ConstSymmMatrixPtr;



typedef size_t mdim_t;

class Matrix : public smartptr::Serializable {

    protected:
        double* data_;
        mdim_t nrow_;
        mdim_t ncol_;
    
    public:
        Matrix(mdim_t nrow, mdim_t ncol);

        Matrix(double* vals, mdim_t nrow, mdim_t ncol);

        Matrix(Vector* v, mdim_t nrow, mdim_t ncol);

        Matrix(const XMLArchivePtr& arch);

        void serialize(const XMLArchivePtr& arch) const;

        ~Matrix();

        double get_element(mdim_t i, mdim_t j) const;

        void set_element(mdim_t i, mdim_t j, double val);

        void accumulate_element(mdim_t i, mdim_t j, double val);

        void scale(double s);

        void assign(const double* vals);

        void assign(double s);

        const double* data() const;

        void printValues(std::ostream& os = std::cout) const;

        void print(const std::string& title, std::ostream& os = std::cout) const;

        Matrix* add(const Matrix* r) const;

        Matrix* subtract(const Matrix* r) const;

        Matrix* mult(double d) const;

        Matrix* mult(const Matrix* m) const;

        Vector* mult(const Vector* m) const;

        mdim_t nrow() const;

        mdim_t ncol() const;

        void accumulate(const Matrix* m);

        Matrix* t() const;

        Vector* toVector() const;

};

class Vector : public smartptr::Serializable {

    protected:
        double* data_;
        mdim_t n_;

    public:
        Vector(double* d, mdim_t n);

        Vector(mdim_t n);

        Vector(const XMLArchivePtr& arch);

        void serialize(const XMLArchivePtr& arch) const;

        ~Vector();

        Vector* add(const Vector* r) const;

        Vector* subtract(const Vector* r) const;

        Vector* mult(double d) const;

        mdim_t n() const;

        double dot(const Vector* v) const;

        void accumulate(const Vector* v);

        void sort();

        void scale(double s);

        void assign(const double* vals);

        void assign(double s);

        void print(const std::string& title, std::ostream& os = std::cout) const;

        void printValues(std::ostream& os = std::cout) const;

        double get_element(mdim_t i) const;

        void set_element(mdim_t i, double val);

        void accumulate_element(mdim_t i, double val);

        const double* data() const;

        Matrix* toMatrix(mdim_t nrow, mdim_t ncol) const;

};

class VectorPtr;

template <class T>
class MatrixTemplate : public boost::intrusive_ptr<T> {
    
    protected:
        typedef boost::intrusive_ptr<T> Parent;

    public:
        using Parent::get;
    
    public:
        MatrixTemplate();

        MatrixTemplate(T* m);

        bool null() const;

        bool nonnull() const;

        void scale(double s);

        double get_element(mdim_t i, mdim_t j) const;

        void printValues(std::ostream& os = std::cout) const;

        void print(const std::string& title, std::ostream& os = std::cout) const;

        void zero();

        mdim_t nrow() const;

        mdim_t ncol() const;

        const double* data() const;

        VectorPtr get_column(mdim_t col) const;

        VectorPtr get_row(mdim_t row) const;

        double maxabs() const;

        VectorPtr toVector() const;

        const VectorPtr toConstVector() const;

        void nullcheck() const;

};

class SymmMatrixPtr;
class RectMatrixPtr;

//non-const version
template <class T>
class RectMatrixTemplate : public MatrixTemplate<T> {
    
    private:
        typedef MatrixTemplate<T> Parent;

    public:
        using Parent::get;
        using Parent::nullcheck;

        RectMatrixTemplate();

        ~RectMatrixTemplate();

        RectMatrixTemplate(T* m);

        RectMatrixTemplate(mdim_t nrow, mdim_t ncol);

        RectMatrixTemplate(mdim_t nrow, mdim_t ncol, Matrix* m);

        RectMatrixTemplate(const RectMatrixPtr& m);

        RectMatrixPtr t() const;

        SymmMatrixPtr symmetrize() const;

        RectMatrixPtr get_subblock(mdim_t rowstart, mdim_t rowstop, mdim_t colstart, mdim_t colstop) const;

        RectMatrixPtr copy() const;

        RectMatrixPtr clone() const;

        void set_element(mdim_t i, mdim_t j, double val);

        void accumulate_element(mdim_t i, mdim_t j, double val);

        void assign_subblock(const ConstRectMatrixPtr& block, mdim_t rowstart, mdim_t colstart);

        void assign(const double* vals);

        void assign(const ConstRectMatrixPtr& m);

        void assign_row(const ConstVectorPtr& v, mdim_t row);

        void assign_column(const ConstVectorPtr& v, mdim_t col);

        void accumulate(const ConstRectMatrixPtr&);

        void eigen(VectorPtr& revals, VectorPtr& ievals, RectMatrixPtr& Levecs, RectMatrixPtr& Revecs) const;

};

template <class T>
class SymmMatrixTemplate : public MatrixTemplate<T> {
    
    private:
        typedef MatrixTemplate<T> Parent;
    
    public:
        using Parent::get;
        using Parent::nullcheck;

        SymmMatrixTemplate();

        SymmMatrixTemplate(mdim_t n);

        SymmMatrixTemplate(T* m);

        SymmMatrixTemplate(const SymmMatrixPtr& m);

        SymmMatrixPtr i(double tol = 1e-8) const;

        SymmMatrixPtr sqrt_matrix(double tol = 1e-8) const;

        SymmMatrixPtr invsqrt_matrix(double tol = 1e-8) const;

        mdim_t n() const;

        void set_element(mdim_t i, mdim_t j, double val);

        void accumulate_element(mdim_t i, mdim_t j, double val);

        void accumulate(const SymmMatrixPtr&);

        void assign(const SymmMatrixPtr& m);

        void eigen(VectorPtr& evals, RectMatrixPtr& evecs) const;

        void accumulate_symmetric_product(const ConstRectMatrixPtr& m);

        void accumulate_transform(const ConstRectMatrixPtr& t, const ConstSymmMatrixPtr& m);

        void accumulate_transform(const ConstSymmMatrixPtr& t, const ConstSymmMatrixPtr& m);

        SymmMatrixPtr clone() const;

        SymmMatrixPtr copy() const;

};

template <class T>
class VectorTemplate : public boost::intrusive_ptr<T> {

    protected:
        typedef boost::intrusive_ptr<T> Parent;

    public:
        using Parent::get;

        VectorTemplate(); 

        VectorTemplate(mdim_t n);

        VectorTemplate(T* v);

        VectorTemplate(const ConstVectorPtr& v);

        double operator[](mdim_t i) const;

        double dot(const ConstVectorPtr& v) const;

        double get_element(mdim_t i) const;

        void set_element(mdim_t i, double val);

        void accumulate_element(mdim_t i, double val);

        void assign(const double* vals);

        void assign(const ConstVectorPtr& v);

        void print(const std::string& title, std::ostream& os = std::cout) const;

        void printValues(std::ostream& os = std::cout) const;

        bool null() const;

        bool nonnull() const;

        mdim_t n() const;

        void normalize();

        void sort();

        void accumulate(const ConstVectorPtr& v);

        void scale(double s);

        void assign(double d);

        VectorPtr clone() const;

        VectorPtr copy() const;

        void zero();

        double maxabs() const;

        double norm() const;

        SymmMatrixPtr symmetric_outer_product() const;

        const double* data() const;

        void nullcheck() const;

};

class RectMatrixPtr : public RectMatrixTemplate<Matrix> {

    private:
        typedef RectMatrixTemplate<Matrix> Parent;

    public:
        using Parent::get;

        RectMatrixPtr();

        RectMatrixPtr(Matrix* m);

        RectMatrixPtr(mdim_t nrow, mdim_t ncol);

        RectMatrixPtr(const RectMatrixPtr& m);

        RectMatrixPtr(const boost::intrusive_ptr<Matrix>& m);

};

class ConstRectMatrixPtr : public RectMatrixTemplate<const Matrix> {

    private:
        typedef RectMatrixTemplate<const Matrix> Parent;

    public:
        using Parent::get;

        ConstRectMatrixPtr();

        ConstRectMatrixPtr(const Matrix* m);

        ConstRectMatrixPtr(const RectMatrixPtr& m);

        ConstRectMatrixPtr(const ConstRectMatrixPtr& m);
    

};

class SymmMatrixPtr : public SymmMatrixTemplate<Matrix> {
    
    private:
        typedef SymmMatrixTemplate<Matrix> Parent;

    public:
        using Parent::get;
        using Parent::nullcheck;

        SymmMatrixPtr();

        SymmMatrixPtr(Matrix* m);

        SymmMatrixPtr(mdim_t n);

        SymmMatrixPtr(const SymmMatrixPtr& m);

        SymmMatrixPtr(const boost::intrusive_ptr<Matrix>& m);

};

class ConstSymmMatrixPtr : public SymmMatrixTemplate<const Matrix> {

    private:
        typedef SymmMatrixTemplate<const Matrix> Parent;

    public:
        using Parent::get;

        ConstSymmMatrixPtr(const Matrix* m);

        ConstSymmMatrixPtr(const SymmMatrixPtr& m);

        ConstSymmMatrixPtr(const ConstSymmMatrixPtr& m);
    
};

class VectorPtr : public VectorTemplate<Vector> {

    private:
        typedef VectorTemplate<Vector> Parent;

    public:
        VectorPtr(); 

        VectorPtr(mdim_t n);

        VectorPtr(Vector* v);

        VectorPtr(const VectorPtr& v);

        VectorPtr(const boost::intrusive_ptr<Vector>& v);

};

class ConstVectorPtr : public VectorTemplate<const Vector> {

    private:
        typedef VectorTemplate<const Vector> Parent;

    public:
        ConstVectorPtr(); 

        ConstVectorPtr(mdim_t n);

        ConstVectorPtr(const Vector* v);

        ConstVectorPtr(const VectorPtr& v);

        ConstVectorPtr(const ConstVectorPtr& v);

};

class Vector1 : public VectorPtr {

    public:
        Vector1(double a0);

};

class Vector2 : public VectorPtr {

    public:
        Vector2(double a0, double a1);

};

class Vector3 : public VectorPtr {

    public:
        Vector3(double a0, double a1, double a2);

};

class Vector4 : public VectorPtr {

    public:
        Vector4(double a0, double a1, double a2, double a3);

};

class Vector8 : public VectorPtr {

    public:
        Vector8(double a0, double a1, double a2, double a3, 
                double a4, double a5, double a6, double a7);

};


RectMatrixPtr operator*(const ConstRectMatrixPtr& l, const ConstRectMatrixPtr& r);

RectMatrixPtr operator*(const ConstRectMatrixPtr& l, const ConstSymmMatrixPtr& r);

RectMatrixPtr operator*(const ConstSymmMatrixPtr& l, const ConstRectMatrixPtr& r);

RectMatrixPtr operator*(const ConstSymmMatrixPtr& l, const ConstSymmMatrixPtr& r);

RectMatrixPtr operator*(const ConstRectMatrixPtr&, double d);

RectMatrixPtr operator*(double d, const ConstRectMatrixPtr&);

RectMatrixPtr operator+(const ConstRectMatrixPtr& l, const ConstRectMatrixPtr& r);

RectMatrixPtr operator-(const ConstRectMatrixPtr& l, const ConstRectMatrixPtr& r);

SymmMatrixPtr operator+(const ConstSymmMatrixPtr& l, const ConstSymmMatrixPtr& r);

SymmMatrixPtr operator-(const ConstSymmMatrixPtr& l, const ConstSymmMatrixPtr& r);

SymmMatrixPtr operator*(const ConstSymmMatrixPtr&, double d);

SymmMatrixPtr operator*(double d, const ConstSymmMatrixPtr&);

VectorPtr operator*(const ConstRectMatrixPtr& m, const ConstVectorPtr& v);

VectorPtr operator*(const ConstSymmMatrixPtr& m, const ConstVectorPtr& v);

VectorPtr operator*(double d, const ConstVectorPtr& v);

VectorPtr operator*(const ConstVectorPtr&, double d);

VectorPtr operator+(const ConstVectorPtr& l, const ConstVectorPtr& r);

VectorPtr operator-(const ConstVectorPtr& l, const ConstVectorPtr& r);

VectorPtr cross(const ConstVectorPtr& l, const ConstVectorPtr& r);

bool equals(const ConstRectMatrixPtr& l, const ConstRectMatrixPtr& r, double tol = 1e-8);

bool equals(const ConstVectorPtr& l, const ConstVectorPtr& r, double tol = 1e-8);

void instantiate();

void
error(
    const Matrix* l,
    const Matrix* r,
    const std::string& op
);

void
error(
    const Vector* l,
    const Vector* r,
    const std::string& op
);

void
scale_array(
    double* v,
    double s,
    mdim_t n
);

void
assign_array(
    double* v,
    double d,
    mdim_t n
);

double 
max_array(
    const double* v,
    mdim_t n
);

void
transpose(
    double* vals,
    mdim_t nrow,
    mdim_t ncol
);

void
eigenvalues(
    const Matrix* matrix,
    Vector* evals,
    Matrix* evecs
);

void
eigenvalues(
    const Matrix* matrix,
    Vector* rvals,
    Vector* ivals,
    Matrix* Levecs,
    Matrix* Revecs
);

double*
subtract_arrays(
    const double* l,
    const double* r,
    mdim_t n
);

double*
add_arrays(
    const double* l,
    const double* r,
    mdim_t n
);

void
accumulate_array(
    double* target,
    const double* src,
    mdim_t n
);

Matrix*
multiply(
    const Matrix* l,
    const Matrix* r,
    bool transpose_l = false,
    bool transpose_r = false
);

template <class T> 
MatrixTemplate<T>::MatrixTemplate()
    : Parent(0)
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
    return !get();
}

template <class T> 
bool
MatrixTemplate<T>::nonnull() const
{
    return get();
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
MatrixTemplate<T>::printValues(std::ostream& os) const
{
    nullcheck();
    get()->printValues(os);
}

template <class T> 
void
MatrixTemplate<T>::print(const std::string& title, std::ostream& os) const
{
    if (get() == 0)
        os << "null matrix" << std::endl;
    else
        get()->print(title, os);
}

template <class T> 
void
MatrixTemplate<T>::nullcheck() const
{
    if (null())
    {
        std::cerr << "Called method on null matrix" << std::endl;
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
    : Parent(0)
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
    : Parent(0)
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
        std::cerr << "Called method on null matrix" << std::endl;
        abort();
    }
}

template <class T>
void
VectorTemplate<T>::printValues(std::ostream& os) const
{
    nullcheck();
    get()->printValues(os);
}

template <class T>
void
VectorTemplate<T>::print(const std::string& title, std::ostream& os) const
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
    return !get();
}

template <class T>
bool
VectorTemplate<T>::nonnull() const
{
    return get();
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
    : Parent(0)
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


}


namespace smartptr {
    SerialDecideConstSubptr(gigide::ConstSymmMatrixPtr);
    SerialDecideConstSubptr(gigide::ConstRectMatrixPtr);
    SerialDecideConstSubptr(gigide::ConstVectorPtr);
    SerialDecideSubptr(gigide::RectMatrixPtr);
    SerialDecideSubptr(gigide::SymmMatrixPtr);
    SerialDecideSubptr(gigide::VectorPtr);
}

#endif
