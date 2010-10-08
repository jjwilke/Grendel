#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <config.h>
#include <stdio.h>
#include <src/printstream.h>
#include <src/utilities.h>
#include <src/fileio.h>
#include <src/qff.h>
#include <src/coordinates.h>
#include <src/keyword.h>
#include <src/defines.h>
#include <src/input.h>
#include <src/symmetry.h>
#include <src/exception.h>
#include <src/displacement.h>
#include <src/derivative.h>
#include <src/deftypes.h>
#include <src/gigmatrix.h>
#include <src/runtime.h>
#include <src/timer.h>

#include <src/smartptr/src/regexp.h>
#include <src/smartptr/src/set.h>

using namespace std;
using namespace gigide;

#define USE_PYTHON 0

void
run_final_intder()
{
    stringstream sstr;
    sstr << INTDER_PATH << " < intder.inp > intder.out";
    int status = system(sstr.str().c_str());
    if (status != 0)
        except("intder failed to run");
}

void
run_final_anharm()
{
    rename("file15", "file15.dat");
    rename("file20", "file20.dat");
    rename("file24", "file24.dat");

    stringstream sstr;
    sstr << ANHARM_PATH;
    remove("output.dat");
    system(sstr.str().c_str());
}

int main(int argc, const char* argv[])
{
    cout << "GRENDEL++ " << VERSION << endl;
    cout << "A code for general energy derivatives for electronic structure methods." << endl;
    cout << " OR " << endl;
    cout << "GIGIDE.  Gigide is a general internal coordinate derivative engine." << endl;
    cout << endl;

#if USE_PYTHON
    Py_Initialize();

    //configure the python path so it knows where to find the python libraries
    PyObject *pmod, *ppath, *pstring;
    pmod = PyImport_ImportModule("sys"); //I could validate, but this should never fail
    ppath = PyObject_GetAttrString(pmod, "path");
    pstring = PyString_FromString(GIGDATADIR);
    PyList_Append(ppath, pstring);

    Py_DECREF(pmod);
    Py_DECREF(ppath);
    Py_DECREF(pstring);
#endif

    gigtimer::Timer::start("gigide");
    try {
        //check the command line options
        map<string, string> argvmap;
        argvmap["--dispcart"] = "Generate dispcart file for given force field";
        argvmap["--calc"] = "Given a data.xml force field file, computes derivatives";
        argvmap["--generate"] = "Given internal coordinate values in keyword 'values' and guess XYZ, generate XYZ coordinates";
        argvmap["--bvectest"] = "Verify analytic B vector equations against numerical differentiation";
        argvmap["--displace"] = "Given XYZ coordinates, generates new XYZ with coordinates displacements given by keyword 'disps' ";
        argvmap["--checkfit"] = "Validate the fitting matrix";
        argvmap["--check"] = "Validate coordinates (characters and non-zero BB^T evals)";
        argvmap["--xmlprint"] = "Load gigide xml files and print details";

        //I am assuming that is never compiled with XLC
        if (argc < 2)
        {
            cerr << "You need to pass a program option so " << PACKAGE_NAME << " knows what to run. Allowed options are: " << endl << endl;
            stringstream sstr;
            print_map(argvmap, sstr);
            except(sstr.str());
        }
        else if (argc > 2)
        {
            cerr << PACKAGE_NAME << " only takes one option. Allowed options are: " << endl << endl;
            stringstream sstr;
            print_map(argvmap, sstr);
            except(sstr.str());
        }

        string option = argv[1];
        if (option == "--help")
        {
            stringstream sstr;
            print_map(argvmap, sstr);
            cout << sstr.str();
            return 0;
        }

        string descr = argvmap[option];
        if (descr.size() == 0) //not a valid option
        {
            argvmap.erase(option);
            cerr << "Invalid option passed to " << PACKAGE << ". Allowed options are: " << endl << endl;
            stringstream sstr;
            print_map(argvmap, sstr);
            except(sstr.str());
        }


        GigideRuntime::init_input();
        GigideRuntime::init_keymap();
        GigideRuntime::init_molecule();
        if (GigideRuntime::molecule()->natoms() == 1) //monoatomic... nothing further to do
        {
            cout << "There is only one atom, and I have nothing left to do.  Goodbye!" << endl;
            return 0;
        }
        GigideRuntime::init_coordinates();
        GigideRuntime::init_symmetry();
        //check precompute options
        if      (option == "--check")    return 0;
        else if (option == "--xmlprint") {GigideRuntime::run_xmlprint();}

        if      (option == "--dispcart")    GigideRuntime::run_dispcart();
        else if (option == "--calc")        {GigideRuntime::run_calc(); }//GigideRuntime::xml_commit();}
        else if (option == "--checkfit")    GigideRuntime::run_checkfit();
        else if (option == "--generate")    GigideRuntime::run_generate();
        else if (option == "--bvectest")    GigideRuntime::run_bvectest();
        else if (option == "--displace")    GigideRuntime::run_displace();

        //commit the xml
        
    }
    catch (GigideException e)
    {
        throw;
    }

    gigtimer::Timer::stop("gigide");
    //commit the xml

#if USE_PYTHON
    Py_Finalize();
#endif

    gigtimer::Timer::print();

    return 0;
}


