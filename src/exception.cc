#include <src/exception.h>
#include <sstream>
#include <iostream>
#include <cstring>

using namespace gigide;
using namespace std;

GigideException::GigideException(
    const string& msg,
    const string& file,
    int line
) throw() : runtime_error(msg), 
    msg_(msg), 
    file_(file), 
    line_(line)
{
}

void
GigideException::rewrite_msg(const string& msg) throw()
{
    msg_ = msg;
}
    
const char*
GigideException::what() const throw()
{
    stringstream sstr;
    sstr << msg_ << "\n";
    //sstr << location();
    sstr << "file: " << file_ << "\n";
    sstr << "line: " << line_;
    sstr << "      " << endl;

    string str = sstr.str();
    char* charstr = new char[str.size()];
    memcpy(charstr, str.data(), str.size() * sizeof(char));
    return charstr;
}

string 
GigideException::file() const throw()
{
    return file_;
}

int 
GigideException::line() const throw()
{
    return line_;
}

string
GigideException::location() const throw()
{
    stringstream sstr;
    sstr << "file: " << file_ << "\n";
    sstr << "line: " << line_;

    const string& str = sstr.str();
    char* charstr = new char[str.size()];
    memcpy(charstr, str.c_str(), str.size() * sizeof(char));
    return charstr;
}

GigideException::~GigideException() throw()
{
}

SanityCheckError::SanityCheckError(
    const string& message,
    const string& file,
    int line
    ) throw() 
  : GigideException(message, file, line)
{
    stringstream sstr;
    sstr << "sanity check failed! " << message;
    rewrite_msg(sstr.str());
}

SanityCheckError::~SanityCheckError() throw() {}

InputException::InputException(
    const string& msg,
    const string& param_name,
    int value,
    const string& file,
    int line
) throw() : GigideException(msg, file, line)
{
    write_input_msg<int>(msg, param_name, value);
}

InputException::InputException(
    const string& msg,
    const string& param_name,
    const string& value,
    const string& file,
    int line
) throw() : GigideException(msg, file, line)
{
    write_input_msg<string>(msg, param_name, value);
}

InputException::InputException(
    const string& msg,
    const string& param_name,
    double value,
    const string& file,
    int line
) throw() : GigideException(msg, file, line)
{
    write_input_msg<double>(msg, param_name, value);
}

InputException::InputException(
    const string& msg,
    const string& param_name,
    const string& file,
    int line
) throw() : GigideException(msg, file, line)
{
    write_input_msg<string>(msg, param_name, "in input");
}

template <
    class T
> void
InputException::write_input_msg(
    const string& msg,
    const string& param_name,
    T value
) throw()
{
    stringstream sstr;
    sstr << msg << "\n";
    sstr << "value " << value << " is incorrect" << "\n";
    sstr << "please change " << param_name << " in input";
    rewrite_msg(sstr.str());
}

StepSizeError::StepSizeError(
    const string& value_name,
    double max,
    double actual,
    const string& file,
    int line) throw()
    : LimitExceeded<double>(value_name + " step size", max, actual, file, line)
{
}

StepSizeError::~StepSizeError() throw() {}

MaxIterationsExceeded::MaxIterationsExceeded(
    const string& routine_name,
    int max,
    const string& file,
    int line)  throw()
    : LimitExceeded<int>(routine_name + " iterations", max, max, file, line)
{
}

MaxIterationsExceeded::~MaxIterationsExceeded() throw() {}

ConvergenceError::ConvergenceError(
    const string& routine_name,
    int max,
    double desired_accuracy,
    double actual_accuracy,
    const string& file,
    int line) throw()
    : MaxIterationsExceeded(routine_name + " iterations", max, file, line), desired_acc_(desired_accuracy), actual_acc_(actual_accuracy)
{
    stringstream sstr;
    sstr << "could not converge " << routine_name << ".  desired " << desired_accuracy << " but got " << actual_accuracy << "\n";
    sstr << description();
    rewrite_msg(sstr.str());
}

ConvergenceError::~ConvergenceError() throw() {}

double 
ConvergenceError::desired_accuracy() const throw() {return desired_acc_;}

double 
ConvergenceError::actual_accuracy() const throw() {return actual_acc_;}

ResourceAllocationError::ResourceAllocationError(
    const string& resource_name,
    size_t max,
    size_t actual,
    const string& file,
    int line)  throw()
    : LimitExceeded<size_t>(resource_name, max, actual, file, line)
{
}

FeatureNotImplemented::FeatureNotImplemented(
    const string& module_name,
    const string& feature_name,
    const string& file,
    int line
) throw() 
 : GigideException("psi exception", file, line)
{
    stringstream sstr;
    sstr << feature_name << " not implemented in " << module_name;
    rewrite_msg(sstr.str());
}

FeatureNotImplemented::~FeatureNotImplemented() throw()
{
}

ResourceAllocationError::~ResourceAllocationError() throw() {}
