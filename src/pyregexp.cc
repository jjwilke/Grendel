#include <stdio.h>
#include <pyregexp.h>
#include <string>
#include <iostream>
#include <sstream>
#include <ostream>
#include <fstream>

using namespace std;
using namespace pyregexp;

#undef heisenbug
#define heisenbug cout << "Heisenbug: " << __FILE__ << " " << __LINE__ << endl

PyObject* 
pyregexp::get_regexp_value(string regexp, string text, string method_name, int flags)
{
    PyObject *preturn, *pmod, *pfunc, *pargs;

    pmod = PyImport_ImportModule("PyRegExp");
	if (!pmod)
	{
		printf("Could not find module PyRegExp\n");
        abort();
	}
    Py_INCREF(pmod);
    pfunc = PyObject_GetAttrString(pmod, const_cast<char*>(method_name.c_str())); 
    Py_INCREF(pfunc);
    pargs = Py_BuildValue("ssi", regexp.c_str(), text.c_str(), flags);
    Py_INCREF(pargs);
    preturn = PyEval_CallObject(pfunc, pargs);
    Py_INCREF(preturn);
    //Py_DECREF(pargs);
    //Py_DECREF(pfunc);
    //Py_DECREF(pmod);
    return preturn;
}

bool 
pyregexp::has_regexp_match(string regexp, string text, int flags)
{
    int ret_val;
    PyObject* preturn = get_regexp_value(regexp, text, "has_match", flags);
    PyArg_Parse(preturn, "i", &ret_val);
    //Py_DECREF(preturn);
    return ret_val;
}

int 
pyregexp::get_regexp_int(string regexp, string text, int flags)
{
    int ret_val;
    PyObject* preturn = get_regexp_value(regexp, text, "get_int_from_regexp", flags);
    PyArg_Parse(preturn, "i", &ret_val);
    //Py_DECREF(preturn);
    return ret_val;
}

int 
pyregexp::get_regexp_num_matches(string regexp, string text, int flags)
{
    int ret_val;
    PyObject* preturn = get_regexp_value(regexp, text, "get_num_matches_from_regexp", flags);
    PyArg_Parse(preturn, "i", &ret_val);
    //Py_DECREF(preturn);
    return ret_val;
}

double 
pyregexp::get_regexp_double(string regexp, string text, int flags)
{
    double ret_val;
    PyObject* preturn = get_regexp_value(regexp, text, "get_double_from_regexp", flags);
    PyArg_Parse(preturn, "d", &ret_val);
    //Py_DECREF(preturn);
    return ret_val;
}

string 
pyregexp::get_regexp_string(string regexp, string text, int flags)
{
    PyObject* preturn = get_regexp_value(regexp, text, "get_string_from_regexp", flags);
    char* ret_val;
    PyArg_Parse(preturn, "s", &ret_val);
    string test = ret_val; //must copy out
    //Py_DECREF(preturn); //destroys the original char, so must copy out
    return test;
}

double* 
pyregexp::get_regexp_double_array(string regexp, string text, size_t& length, int flags)
{
    PyObject* preturn = get_regexp_value(regexp, text, "get_double_array_from_regexp", flags);

    length = PyList_Size(preturn); 
    double* vals = new double[length];
    for (size_t i=0; i < length; ++i)
    {
        PyObject* obj = PyList_GetItem(preturn, i);
        double val;
        PyArg_Parse(obj, "d", &val);
        vals[i] = val;
    }

    //Py_DECREF(preturn);
    return vals;
}

int* 
pyregexp::get_regexp_int_array(string regexp, string text, size_t& length, int flags)
{
    PyObject* preturn = get_regexp_value(regexp, text, "get_int_array_from_regexp", flags);
    
    length = PyList_Size(preturn); 
    int* vals = new int[length];
    for (size_t i=0; i < length; ++i)
    {
        PyObject* obj = PyList_GetItem(preturn, i);
        int val;
        PyArg_Parse(obj, "i", &val);
        vals[i] = val;
    }

    //Py_DECREF(preturn);
    return vals;
}

void
pyregexp::findmatch(vector<string>& matches, string regexp, string text, int numgroups, int flags)
{
    if ( (flags >> FindAllShift) % 2 )
    {
        vector< vector<string> > raw_matches; 
        findall(raw_matches, regexp, text, numgroups, flags);
        group_matches(raw_matches, matches);
    }
    else
    {
        PyObject *pdict, *pfunc, *pmod, *pargs, *pval;
        PyGILState_STATE gstate;
        gstate = PyGILState_Ensure();
        pdict = PyDict_New();
        Py_INCREF(pdict);
        pmod = PyImport_ImportModule("PyRegExp");
        if (!pmod)
        {
            printf("Could not find module PyRegExp\n");
        }
        Py_INCREF(pmod);
        pfunc = PyObject_GetAttrString(pmod, "build_match_from_regexp");
        Py_INCREF(pfunc);
        pargs = Py_BuildValue("ssOi", regexp.c_str(), text.c_str(), pdict, flags);
        Py_INCREF(pargs);
        PyEval_CallObject(pfunc, pargs);

        for (int i=0; i < numgroups; i++)
        {
            //allow for 10,000 matches
            char* key = new char[5]; 
            char* newval;
            sprintf(key, "%d", i);
            pval = PyDict_GetItemString(pdict, key);
            Py_INCREF(pval);
            PyArg_Parse(pval, "s", &newval);
            matches.push_back(newval);
            //Py_DECREF(pval);
            delete[] key;
        }
        //Py_DECREF(pdict);
        //Py_DECREF(pargs);
        //Py_DECREF(pfunc);
        //Py_DECREF(pmod);
        PyGILState_Release(gstate);
    }
}

void
pyregexp::findall(vector<vector<string> >& matches, string regexp, string text, int numgroups, int flags)
{
    PyObject *pdict, *pfunc, *pmod, *pargs, *pval, *preturn;
    int nummatches = 0;

    pdict = PyDict_New();
    Py_INCREF(pdict);
    pmod = PyImport_ImportModule("PyRegExp");
	if (!pmod)
	{
		printf("Could not find module PyRegExp\n");
	}
    Py_INCREF(pmod);
    pfunc = PyObject_GetAttrString(pmod, "build_findall_from_regexp");
    Py_INCREF(pfunc);
    pargs = Py_BuildValue("ssOi", regexp.c_str(), text.c_str(), pdict, flags);
    Py_INCREF(pargs);
    preturn = PyEval_CallObject(pfunc, pargs);
    PyArg_Parse(preturn, "i", &nummatches);

    for (int match=0; match < nummatches; match++)
    {
        vector<string> nextmatch;
        for (int group=0; group < numgroups; group++)
        {
            //allow for 10,000 matches
            char* key = new char[5]; 
            char* newval;
            sprintf(key, "%d,%d", match, group);
            pval = PyDict_GetItemString(pdict, key);
            Py_INCREF(pval);
            PyArg_Parse(pval, "s", &newval);
            nextmatch.push_back(newval); //this copies the char* so there are no memory leak issues
            //Py_DECREF(pval);
            delete[] key;
        }
        matches.push_back(nextmatch);
    }
    //Py_DECREF(pdict);
    //Py_DECREF(pargs);
    //Py_DECREF(pfunc);
    //Py_DECREF(preturn);
}

void
pyregexp::group_matches(vector< vector<string> >& raw_matches, vector<string>& matches)
{
    for (int match=0; match < raw_matches.size(); ++match)
    {
        for (int group=0; group < raw_matches[match].size(); ++group)
        {
            matches.push_back( raw_matches[match][group] );
        }
    }
}


