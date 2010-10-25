#include <src/defines.h>
#include <src/utilities.h>
#include <src/exception.h>
#include <sstream>
#include <cmath>

#include <src/smartptr/src/printstream.h>

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



