#ifndef gigide_matrix_h
#define gigide_matrix_h

#include <iostream>
#include <src/smartptr/src/serialize.h>

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

        void print(std::string title, std::ostream& os = std::cout) const;

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

        void print(std::string title, std::ostream& os = std::cout) const;

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

        void print(std::string title, std::ostream& os = std::cout) const;

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

        void print(std::string title, std::ostream& os = std::cout) const;

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

}

namespace smartptr {
    SerialDecideSubptr(gigide::RectMatrixPtr);
    SerialDecideSubptr(gigide::SymmMatrixPtr);
    SerialDecideSubptr(gigide::VectorPtr);
}

#endif
