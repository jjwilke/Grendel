#include <src/runtime.h>
#include <src/timer.h>


using namespace std;
using namespace gigide;
using namespace smartptr;

void
GigideRuntime::run_intder()
{
    stringstream sstr;
    sstr << INTDER_PATH << " < intder.inp > intder.out";
    int status = system(sstr.str().c_str());
    if (status != 0)
        except("intder failed to run");
}

void
GigideRuntime::run_anharm()
{
    rename("file15", "file15.dat");
    rename("file20", "file20.dat");
    rename("file24", "file24.dat");

    stringstream sstr;
    sstr << ANHARM_PATH;
    remove("output.dat");
    system(sstr.str().c_str());
}

void
GigideRuntime::xml_commit()
{
    gigtimer::Timer::start("xml commit");
    //reset the writer
    archive_ = new smartptr::XMLArchive();

    serial_call_save(archive_, mol_, "molecule");
    serial_call_save(archive_, simples_, "simples");
    serial_call_save(archive_, coords_, "coordinates");
    serial_call_save(archive_, qff_, "qff");
    archive_->toFile("gigide.xml");
    gigtimer::Timer::stop("xml commit");
}

void
GigideRuntime::init_molecule()
{
    //build the molecule
    vector<string> atomlist;
    RectMatrixPtr xyz = input_->getGeometry(atomlist);
    stringstream sstr;
    string units = KeywordSet::getKeyword("bond units")->getValueString();
    sstr << "Coordinates (" << units << ")";
    xyz.print(sstr.str().c_str());
    mol_ = new Molecule(xyz, atomlist);
    input_->readSymmetryOperations(mol_);
    mol_->computePointGroup();
    cout << "Abelian point group: " << mol_->getPointGroup()->name() << endl;
    mol_->getPointGroup()->print();
    cout << endl;
}

void
GigideRuntime::init_coordinates()
{
    //setup the coordinates
    simples_.clear(); coords_.clear();
    input_->getSimpleCoordinates(mol_, simples_);
    if (simples_.size() == 0) //no coordinates!
        except("No simple internal coordinates given!");
}

void
GigideRuntime::init_symmetry()
{
    mol_->getPointGroup()->formBasis(simples_);

    {
        vector<SimpleInternalCoordinatePtr>::iterator it;
        for (it = simples_.begin(); it != simples_.end(); ++it)
        {
            SimpleInternalCoordinatePtr coord = *it;
            coord->computeCharacters();
            coord->print();
            cout << endl;
        }
    }

    input_->getSymmetryCoordinates(mol_, simples_, coords_);
    VectorPtr bvals = InternalCoordinate::checkCompleteness(coords_, mol_);
    bvals.sort();
    bvals.print("BB^T eigenvalues");

    cout << stream_printf("%53s", ""); 
    mol_->getPointGroup()->printClasses();
    cout << endl;

    {
        vector<InternalCoordinatePtr>::iterator it;
        int i=0;
        for (it = coords_.begin(); it != coords_.end(); ++it, ++i)
        {
            InternalCoordinatePtr coord = *it;
            //JJW 04/13/2010 for computing characters, the coordinate needs to be XYZ normalized
            //intder requires the coefficients to be normalized
            //for future versions this will need to change
            coord->intder_normalize();
            coord->computeCharacters();
            coord->printDetail();
        }
    }
}

void
GigideRuntime::init_qff()
{
    //and build the quartic force field
    qff_ = new ForceField(mol_, coords_, simples_);
}

void
GigideRuntime::init_input()
{
    //validate the input file to make sure all of the
    //sections are input properly
    string filetext = getFileText(DEFAULT_INPUT_FILE, "#"); //use number sign as comment
    input_ = new GigideInputFile(filetext);
    input_->validateInputFile();
}

void
GigideRuntime::init_keymap()
{
    //set up the options map
    keymap_ = input_->getOptions();
    cout << "Input options: " << endl;
    keymap_->print();
} 

ConstMoleculePtr
GigideRuntime::molecule()
{
    return mol_;
}

void
GigideRuntime::run_xmlprint()
{
    gigtimer::Timer::start("xml");
    XMLArchivePtr arch(new smartptr::XMLArchive("gigide.xml"));

    MoleculePtr testmol;
    serial_call_load(arch, testmol, "molecule");
    testmol->print();

    vector<SimpleInternalCoordinatePtr> simples;
    serial_call_load(arch, simples, "simples");
    for (int i=0; i < simples.size(); ++i)
    {
        simples[i]->print();
        cout << endl;
    }

    vector<InternalCoordinatePtr> coords;
    serial_call_load(arch, coords, "coordinates");
    for (int i=0; i < coords.size(); ++i)
    {
        coords[i]->print();
        cout << endl;
    }

    testmol->getPointGroup()->print();

    
    ForceFieldPtr qff;
    serial_call_load(arch, qff, "qff");
    if (!qff)
        except("Could not find force field on tag qff");
    gigtimer::Timer::stop("xml");

    qff->print();
    gigtimer::Timer::start("compute");
    qff->compute();
#if HAS_INTDER
    run_intder();
#endif
#if HAS_INTDER && HAS_ANHARM
    run_anharm();
#endif
    gigtimer::Timer::stop("compute");

}

void
GigideRuntime::run_dispcart()
{
    if (!qff_)
        init_qff();

    qff_->generateDisplacements();
    qff_->writeDisplacementsToFile("dispcart");
}

void
GigideRuntime::run_calc()
{
    if (!qff_)
        init_qff();

    gigtimer::Timer::start("read xml data");
    qff_->readXMLData();
    gigtimer::Timer::stop("read xml data");
    gigtimer::Timer::start("compute");
    qff_->compute();
#if HAS_INTDER
    run_intder();
#endif
#if HAS_INTDER && HAS_ANHARM
    run_anharm();
#endif
    gigtimer::Timer::stop("compute");

    xml_commit();
}

void
GigideRuntime::run_checkfit()
{
    qff_->validateFit();
}

void
GigideRuntime::run_generate()
{
    vector<double> values = KeywordSet::getKeyword("values")->getValueVectorDouble();
    MoleculePtr dispmol = Displacement::generateGeometry(values, mol_, coords_, simples_);
    dispmol->getXYZ().print("new geometry");
}

void
GigideRuntime::run_bvectest()
{
    for (int i=0; i < coords_.size(); ++i)
    {
        coords_[i]->testBVectors();
    }
}

void
GigideRuntime::run_displace()
{
    vector<double> values = KeywordSet::getKeyword("disps")->getValueVectorDouble();
    MoleculePtr mol = Displacement::displaceGeometry(values, mol_, coords_, simples_);
    mol->getXYZ().print("new geometry");
}

XMLArchivePtr GigideRuntime::archive_ = 0;
MoleculePtr GigideRuntime::mol_ = 0;
std::vector<SimpleInternalCoordinatePtr> GigideRuntime::simples_;
std::vector<InternalCoordinatePtr> GigideRuntime::coords_;
ForceFieldPtr GigideRuntime::qff_ = 0;
GigideInputFilePtr GigideRuntime::input_ = 0;
KeywordSetPtr GigideRuntime::keymap_ = 0;

