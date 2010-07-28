#ifndef pyregexp_h
#define pyregexp_h

#include <Python.h>
#include <string>
#include <vector>

namespace pyregexp {

enum PyFlags { IgnoreCase=1, UpperCase=2, FindAll=4, LowerCase=8, IncludeLineBreaks=16, StripWhitespace=32 };
enum Shifts { IgnoreCaseShift=0, UpperCaseShift=1, FindAllShift=2, LowerCaseShift=3, IncludeLineBreaksShift=4, StripWhitespaceShift=5 };

PyObject* get_regexp_value(std::string regexp, std::string text, std::string method_name, int flags = 0);

int get_regexp_int(std::string regexp, std::string text, int flags = 0);
 
double get_regexp_double(std::string regexp, std::string text, int flags = 0);

std::string get_regexp_string(std::string regexp, std::string text, int flags = 0);

double* get_regexp_double_array(std::string regexp, std::string text, size_t& length, int flags = 0);

int* get_regexp_int_array(std::string regexp, std::string text, size_t& length, int flags = 0);

void
findmatch(std::vector<std::string>& matches, std::string regexp, std::string text, int numgroups, int flags = 0);

void
findall(std::vector<std::vector<std::string> >& matches, std::string regexp, std::string text, int numgroups, int flags = 0);

void
group_matches(std::vector< std::vector<std::string> >& raw_matches, std::vector<std::string>& matches);

bool has_regexp_match(std::string regexp, std::string text, int flags = 0);

int get_regexp_num_matches(std::string regexp, std::string text, int flags = 0);

}

#endif
