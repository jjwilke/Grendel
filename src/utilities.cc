#include <src/defines.h>
#include <src/utilities.h>
#include <src/printstream.h>
#include <src/exception.h>
#include <sstream>
#include <cmath>

using namespace gigide;
using namespace std;

int 
gigide::factorial(int x)
{
    int fact = 1;
    for (int i=2; i <= x; i++)
        fact *= i;
    return fact;
}

int
gigide::roundNumber(double x)
{
    if (x < 0)
    {
        return -1 * lrint(-x);
    }
    else
        return lrint(x);
}

bool
gigide::getFraction(
    double val,
    int& numerator,
    int& denom
)
{
    double testval, diff;
    for (int n=1; n <= 9999; ++n)
    {
        testval = val * n;
        int num = round(testval);
        diff = fabs(testval - num); //see if this multiplication produced an integer
        if (diff < 1e-6)
        {
            numerator = num;
            denom = n;
            return true;
        }
    }
    return false;
}

void
gigide::print_vector(vector<int>& v, bool newline, const char* format)
{
    cout << "<";
    for (int i=0; i < v.size(); i++)
    {
        cout << stream_printf(format, v[i]) << ",";
    }
    cout << ">";
    if (newline) cout << endl;
}

void
gigide::print_vector(vector<double>& v, bool newline, const char* format)
{
    cout << "<";
    for (int i=0; i < v.size(); i++)
    {
        cout << stream_printf(format, v[i]) << ",";
    }
    cout << ">";
    if (newline) cout << endl;
}

double
gigide::inverse_cosine(double value)
{
    if      ( fabs(value - 1) < 1e-12 ) return 0.0;
    else if ( fabs(value + 1) < 1e-12 ) return PI;
    else return acos(value);
}

double
gigide::inverse_sine(double value)
{
    if      ( fabs(value - 1) < 1e-12 ) return PI/2.0;
    else if ( fabs(value + 1) < 1e-12 ) return 3*PI/2.0;
    else return asin(value);
}

void
assignValues(string text, Matrix* matrix)
{
    int nrow = matrix->nrow();
    int ncol = matrix->ncol();

    stringstream sstr(text);
    for (int i=0; i < nrow; ++i)
    {
        for (int j=0; j < ncol; ++j)
        {
            double val; sstr >> val;
            matrix->set_element(i,j,val);
        }
    }
}

#if 0
ArchivePtr
saveXMLMatrix(const Matrix* matrix, const ArchivePtr& arch, string tagname)
{
    ArchivePtr node = arch->writeBinary(matrix->data(), sizeof(double) * matrix->nrow() * matrix->ncol(), tagname);
    node->attribute(matrix->nrow(), "nrow");
    node->attribute(matrix->ncol(), "ncol");
    return node;
}

ArchivePtr
gigide::saveMatrix(const ConstSymmMatrixPtr& matrix, const ArchivePtr& arch, string tagname)
{
    return saveXMLMatrix(matrix.get(), arch, tagname);
}

ArchivePtr
gigide::saveMatrix(const ConstRectMatrixPtr& matrix, const ArchivePtr& arch, string tagname)
{
    return saveXMLMatrix(matrix.get(), arch, tagname);
}

ArchivePtr
gigide::saveVector(const ConstVectorPtr& v, const ArchivePtr& arch, std::string tagname)
{
    ArchivePtr node = arch->writeBinary(v->data(), sizeof(double) * v->n(), tagname);
    node->setAttribute(v.n(), "n");
    return node;
}

ConstArchivePtr
gigide::loadMatrix(RectMatrixPtr& m, const ConstArchivePtr& arch, std::string tagname)
{
    void* vals;
    size_t size;
    ConstArchivePtr node = arch->loadBinary(&vals, size, tagname);
    if (node.get() == NULL)
        return NULL;

    
    int nrow; node->attribute(nrow, "nrow");
    int ncol; node->attribute(ncol, "ncol");
    
    if (nrow != ncol)
        except(stream_printf("Matrix on tagname %s is not symmetric", tagname.c_str()));

    m = new Matrix((double *) vals, nrow, ncol);
    return node;
}

ConstArchivePtr
gigide::loadMatrix(SymmMatrixPtr& s, const ConstArchivePtr& arch, std::string tagname)
{
    void* vals;
    size_t size;
    ConstArchivePtr node = arch->loadBinary(&vals, size, tagname);
    if (node.get() == NULL)
        return NULL;

    
    int nrow; node->attribute(nrow, "nrow");
    int ncol; node->attribute(ncol, "ncol");
    
    if (nrow != ncol)
        except(stream_printf("Matrix on tagname %s is not symmetric", tagname.c_str()));

    s = new Matrix((double *) vals, nrow, ncol);
    return node;
}

ConstArchivePtr
gigide::loadVector(VectorPtr& v, const ConstArchivePtr& arch, std::string tagname)
{
    void* vals;
    size_t size;
    ConstArchivePtr node = arch->loadBinary(&vals, size, tagname);
    if (node.get() == NULL)
        return NULL;

    int n; node->getAttribute(n, "n");
    
    v = new Vector((double *) vals, n);
    return node;
}
#endif


