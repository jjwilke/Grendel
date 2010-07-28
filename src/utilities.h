
#ifndef _gigide_utils_h_
#define _gigide_utils_h_

#include <vector>
#include <src/gigmatrix.h>
#include <src/archive.h>

namespace gigide {

/**
    @param
    @return x!
*/
int factorial(int x);

/**
    @param x
    @return Nearest integer
*/
int roundNumber(double x);

/**
    @param value
    @return The inverse cosine between 0 and \f$\pi\f$
*/
double inverse_cosine(double value);

/**
    @param value
    @return The inverse sin between \f$ \frac{\pi}{2} \f$ and \f$ \frac{3\pi}{2} \f$
*/
double inverse_sine(double value);


/**
    @param val The number to convert to rational form
    @param num Will hold the numerator upon return
    @param denom Will hold the denominator upon return
    @return Whether the build succeeded with a max denominator of 10,000
*/
bool
getFraction(
    double val,
    int& num,
    int& denom
);

/**
    @param v The vector of values to print
    @param newline Whether to include a trailing newline
    @param format A formatting string for the integers
*/
void 
print_vector(std::vector<int>& v, bool newline = true, const char* format = "%2d");

/**
    @param v The vector of values to print
    @param newline Whether to include a trailing newline
    @param format A formatting string for the doubles
*/
void 
print_vector(std::vector<double>& v, bool newline = true, const char* format = "%7.4f");

#if 0
ArchivePtr saveMatrix(const ConstSymmMatrixPtr& matrix, const ArchivePtr& xml, std::string tagname);

ArchivePtr saveMatrix(const ConstRectMatrixPtr& matrix, const ArchivePtr& xml, std::string tagname);

ArchivePtr saveVector(const ConstVectorPtr& v, const ArchivePtr& xml, std::string tagname);

ConstArchivePtr loadMatrix(SymmMatrixPtr& s, const ConstArchivePtr& xml, std::string tagname);

ConstArchivePtr loadMatrix(RectMatrixPtr& m, const ConstArchivePtr& xml, std::string tagname);

ConstArchivePtr loadVector(VectorPtr& v, const ConstArchivePtr& xml, std::string tagname);
#endif

}

#endif
