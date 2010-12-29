#include <config.h>
#include <iostream>
#include <algorithm>

#include "gigide.h"

using namespace std;
using namespace gigide;
using namespace smartptr;


string gigide::getFileText(const std::string& filename, const std::string& comment)
{
    ifstream infile(filename.c_str());

    if (!infile.good())
    {
        except(stream_printf("%s is not a valid file", filename.c_str()));
    }

    string filetext;
    string comment_regexp = "\\s*["; comment_regexp += comment; comment_regexp += "]";
    string uncomment_regexp = "\\s*(.*?)["; uncomment_regexp += comment; uncomment_regexp += "]";
    char* nextline = new char[MAX_LINE_SIZE + 5];
    //read the first line without adding trailing newline
    if (!infile.eof())
    {
        infile.getline(nextline, MAX_LINE_SIZE);
        string line = nextline;
        //see if it is a comment line
        if (comment.size() && has_regexp_match(comment_regexp, line))
            line = get_regexp_string(uncomment_regexp, line);

        if (line.size())
            filetext += line;
    }
    //read the rest of the file inserting new line
    //char at beginning
    while (!infile.eof())
    {
        filetext += "\n";
        infile.getline(nextline, MAX_LINE_SIZE);
        string line = nextline;
        //see if it is a comment line
        if (comment.size() && has_regexp_match(comment_regexp, line)) //only do if we have been given a comment
            line = get_regexp_string(uncomment_regexp, line);

        if (line.size())
            filetext += line;
    }

    delete[] nextline;

    return filetext;
}

FileWriter::FileWriter()
{
}

void
FileWriter::addEndl()
{
    contents_ << endl;
}

void
FileWriter::commit(const std::string& filename)
{
    ofstream outfile;
    outfile.open(filename.c_str());
    outfile << contents_.str() << endl;
    outfile.close();
}

void
FileWriter::add(const std::string& newtext)
{
    contents_ << newtext;
}

InputFileWriter::InputFileWriter(const ConstMoleculePtr& mol)
    : mol_(mol)
{
}

void
InputFileWriter::makeFile11(const ConstMoleculePtr& mol, ConstRectMatrixPtr grads)
{
    ofstream outfile;
    outfile.open("file11");
    outfile << "Transformation" << endl;
    outfile << stream_printf("%5d%15.7f", mol->natoms(), 0.0) << endl;

    double unit_conversion = Molecule::getConversion("bohr"); //must be in bohr
    
    //print geometry
    ConstRectMatrixPtr xyz = mol->getXYZ(); 
    for (int n=0; n < mol->natoms(); n++)
    {
        outfile << stream_printf("%20.10f", mol->getAtom(n)->mass());
        for (int x=0; x < XYZ_DIM; x++)
            outfile << stream_printf("%20.10f", xyz.get_element(n,x) * unit_conversion);
        outfile << endl;
    }

    //print gradients - this is only ever hartree/bohr
    for (int n=0; n < mol->natoms(); n++)
    {
        outfile << stream_printf("%20s","");
        for (int x=0; x < XYZ_DIM; x++)
            outfile << stream_printf("%20.10f", grads.get_element(n,x));
        outfile << endl;
    }
    outfile.close();
}

void
InputFileWriter::writeHeader()
{
}

void
InputFileWriter::writeFooter()
{
    addEndl();
}

void
InputFileWriter::closeSection()
{
}

void
InputFileWriter::makeFile12(
    const ConstMoleculePtr& mol, 
    const Set<ConstInternalCoordinatePtr>& coords,
    ConstVectorPtr grads
)
{
    ofstream outfile;
    outfile.open("file12");
    outfile << stream_printf("%5d%22.10f", mol->natoms(), 0.0) << endl;

    
    //print coordinates and values
    //double unit_conversion = Molecule::getConversion("bohr");
    //ConstRectMatrixPtr xyz = mol->getXYZ(); //must be in bohr
    outfile << stream_printf("%20.10f%20.10f", coords[0]->getValueForMolecule(mol), grads.get_element(0));
    for (int n=1; n < coords.size(); ++n)
    {
        outfile << endl;
        outfile << stream_printf("%20.10f%20.10f", coords[n]->getValueForMolecule(mol), grads.get_element(n));
    }
    outfile << endl; //end of file error on certain intders, absolutely necessary
    outfile.close();
}

void
InputFileWriter::makeFCFile(
    const ConstMoleculePtr& mol, 
    ConstRectMatrixPtr fc, 
    const std::string& filename
)
{
    ofstream outfile;
    outfile.open(filename.c_str());
    outfile << stream_printf("%5d%5d", mol->natoms(), mol->natoms() * 6) << endl;

    mdim_t size = fc.nrow() * fc.ncol();
    const double* vals = fc.data(); 
    const double* valptr = vals;

    for (mdim_t row=0; row < size / 3; ++row)
    {
        for (mdim_t col=0; col < 3; ++col, ++valptr)
        {
            outfile << stream_printf("%20.10f", *valptr);
        }
        outfile << endl;
    }
    outfile.close();
}

void
InputFileWriter::makeFile15(
    const ConstMoleculePtr& mol, 
    ConstRectMatrixPtr fc
)
{
    makeFCFile(mol, fc, "file15");
}

void
InputFileWriter::makeFile16(
    const ConstMoleculePtr& mol, 
    ConstRectMatrixPtr fc
)
{
    makeFCFile(mol, fc, "file16");
}

IntderWriter::IntderWriter(
             const ConstMoleculePtr& mol,
             bool statpt,
             TransformationType type,
             int nder,
             const ConstDerivativeIteratorPtr& iter,
             const Set<ConstSimpleInternalCoordinatePtr>& simples,
             const Set<ConstInternalCoordinatePtr>&  coords,
             bool includexyz)
    : InputFileWriter(mol), simples_(simples), coords_(coords)
{
    type_ = type;
    statpt_ = statpt;
    if (!coords.size() || !simples.size())
    {
        except("No coordinates given to intder");
    }

    iter_ = iter;
    nder_ = nder;
    includexyz_ = includexyz;

    ConstSymmetryInternalCoordinatePtr coord = boost::dynamic_pointer_cast<const SymmetryInternalCoordinate,const InternalCoordinate>(coords_[0]);
    has_symms_ = coord;
}

void
IntderWriter::writeHeader()
{
    stringstream out;
    out << "# FILES ##################" << endl << endl;
    out << "TITLE" << endl << endl;
    out << "# INTDER #################" << endl;

    //set the intder options
    int* options = new int[16 + 1]; //16 intder options... start counting at 1 for readaility
    options[1] = mol_->natoms(); //number of atoms
    options[2] = simples_.size(); //number of simples
    options[3] = has_symms_ ? coords_.size() : 0; //number of symm internals
    options[4] = nder_; //derivative level
    options[5] = statpt_ ? 0 : 1; //whether at stat point
    options[6] = 2000; //print options
    options[7] = (int) type_; //option telling intder to transform internals to cartesians
    options[8] = mol_->ndummies(); //number of dummy atoms 
    options[9] = 0; //some testing flag
    options[10] = includexyz_ ? 1 : 0; //1 says read geom from input file
    options[11] = nder_ >= 2 ? 3 : 0; //perform freq analysis in carts and internals, if we have 2nd or higher derivatives
    options[12] = 0;
    options[13] = 0;
    options[14] = 0;
    options[15] = 0;
    options[16] = 0;
    for (int iopt=1; iopt <= 16; iopt++)
        out << stream_printf("%5d", options[iopt]);

    add(out.str());

    delete[] options;
}

void
IntderWriter::writeFooter()
{
}

bool
intder_sort(int i, int j){return i > j;}

void
IntderWriter::writeDerivatives()
{
    if (!iter_) //ther are no internal coordinate derivatives to write
        return;

    stringstream out;
    //now print the derivatives
    int nstart = statpt_ ? 2 : 1; //depnding on whether or not we have a stationary point, we might include 1st derivs
    for (int n=nstart; n <= nder_; n++) //dep
    {
        //loop through the derivatives
        int nderiv = 0;
        DerivativeIterator::const_iterator it;
        for (it = iter_->begin(); it != iter_->end(); ++it, nderiv++)
        {
            ConstDerivativePtr der = *it;
            if (der->level() == n)
            {
                vector<int> indices = der->indices();
                sort(indices.begin(), indices.end(), intder_sort);
                for (int i=0; i < indices.size(); i++)
                    out << stream_printf("%5d", indices[i] + 1); //intder uses 1 based counting
                //fill out the spaces to be four
                for (int i=indices.size(); i < 4; i++)
                    out << stream_printf("%5s", "");
                out << stream_printf("%20.10f", der->value());
                out << endl;
            }
        }
        out << stream_printf("%5d", 0) << endl;
    }
    add(out.str());
}

void
IntderWriter::writeSimpleCoordinates()
{
    int lin1_counter = mol_->natoms();
    stringstream out;
    //print the simple internals information
    for (int i=0; i < simples_.size(); i++)
    {
        out << endl;
        ConstSimpleInternalCoordinatePtr coord = simples_[i];
        out << stream_printf("%5s", coord->type().c_str());
        vector<int> connect; coord->connectivity(connect);
        for (int atom=0; atom < connect.size(); atom++) 
            out << stream_printf("%5d", connect[atom]);
        vector<int> dummies; coord->dummies(dummies);
        for (int atom=0; atom < dummies.size(); atom++) 
            out << stream_printf("%5d", dummies[atom]);
    }
    add(out.str());
}

void
IntderWriter::writeSymmetryCoordinates()
{
    stringstream out;
    for (int i=0; i < coords_.size(); i++)
    {
        ConstSymmetryInternalCoordinatePtr coord = boost::dynamic_pointer_cast<const SymmetryInternalCoordinate,const InternalCoordinate>(coords_[i]);
        if (!coord) //these aren't symmetry internals, these are simple internals
        {
            addEndl();
            return;
        }

        int colwidth = 100; //start it as too large
        int coordnum = i + 1;
        ConstVectorPtr coeffs = coord->getCoefficients();
        for (int simple=0; simple < coeffs.n(); ++simple)
        {
            int simplenum = coord->getSimpleCoordinateNumber(simple, simples_) + 1;
            if (colwidth > 62)
            {
                out << endl << stream_printf("%5d%4d", coordnum, simplenum);
                colwidth = 8;
            }
            else
            {
                out << stream_printf("%6d", simplenum);
                colwidth += 5;
            }
            out << stream_printf("%12.8f", coeffs.get_element(simple));
            colwidth += 12;
        }
    }
    out << endl;
    add(out.str());
    closeSection();
}

void
IntderWriter::writeXYZ()
{
    if (!includexyz_)
        return;

    stringstream out;
    //now print out the atom coordinates
    double unit_conversion = Molecule::getConversion("bohr");
    ConstRectMatrixPtr xyz = mol_->getXYZ(); //must be in bohr
    for (mdim_t n=0; n < xyz.nrow(); n++)
    {
        for (int x=0; x < XYZ_DIM; x++)
            out << stream_printf("%20.10f", xyz.get_element(n,x) * unit_conversion);
        out << endl;
    }

    for (int n=0; n < mol_->ndummies(); ++n)
    {
        ConstVectorPtr xyz = mol_->getDummy(n)->getXYZ();
        for (int x=0; x < XYZ_DIM; x++)
            out << stream_printf("%20.10f", xyz.get_element(x) * unit_conversion);
        out << endl;
    }

    add(out.str());
}

void
IntderWriter::writeMasses()
{
    if (!includexyz_)
        return;

    stringstream out;
    //now print the masses
    int colwidth = 12;
    for (int n=0; n < mol_->natoms(); n++)
    {
        out << stream_printf("%12.6f", mol_->getAtom(n)->mass());
        colwidth += 12;
        if (colwidth > 80) //don't exceed 80 char width
        {
            colwidth = 12;
            out << endl;
        }
    }
    out << endl;
    add(out.str());
}

void
IntderWriter::closeSection()
{
    add("    0\n");
}

void
IntderWriter::commit(const std::string& filename)
{
    writeHeader();
    writeSimpleCoordinates();
    writeSymmetryCoordinates();
    writeXYZ();
    writeMasses();
    writeDerivatives();
    writeFooter();
    InputFileWriter::commit(filename);
}

void
IntderWriter::clearFiles()
{
    remove("file11");
    remove("file12");
    remove("file15");
    remove("file16");
    remove("intder.inp");
    remove("intder.out");
}

VectorPtr
IntderWriter::transformGradientsToInternals(
    ConstRectMatrixPtr xyzgrads,
    const ConstMoleculePtr& mol,
    const Set<ConstSimpleInternalCoordinatePtr>& simples,
    const Set<ConstInternalCoordinatePtr>& coords,
    bool validate
)
{
    clearFiles();

    //transform the gradients using intder
    int derlevel = 1;
    bool statpt = false;
    bool incxyz = false;
    IntderWriter iderwriter(mol, statpt, IntderWriter::TransformCartesians, derlevel, 0, simples, coords, incxyz);
    iderwriter.commit(DEFAULT_INTDER_INPUT_FILE);
    InputFileWriter::makeFile11(mol, xyzgrads);

    stringstream sstr;
    sstr << INTDER_PATH << " < intder.inp > intder.out";
    int status = system(sstr.str().c_str());
    if (status != 0)
    {
        except("Intder returned a status error on gradient transformation to internals");
    }

    string file12 = getFileText("file12", ""); //pass blank comment, no comments
    string regexp = "[-]?\\d+[.]\\d+[ ]+([-]?\\d+[.]\\d+)\\n";
    size_t nummatches;
    double* gradvals = get_regexp_double_array(regexp, file12, nummatches);
    VectorPtr gradients(nummatches); 
    gradients.assign(gradvals);
    delete[] gradvals;

    if (validate)
    {
        RectMatrixPtr check = transformGradientsToCartesian(gradients, mol, simples, coords, false);

        if (!equals(check, xyzgrads, 1e-8))
        {
            cout << "Internal gradient transform failed" << endl;
            xyzgrads.print("Correct");
            check.print("Incorrect");
            abort();
        }
    }

    return gradients;
}

RectMatrixPtr
IntderWriter::transformForceConstantsToInternals(
    ConstRectMatrixPtr xyzfc,
    ConstRectMatrixPtr xyzgrads,
    const ConstMoleculePtr& mol,
    const Set<ConstSimpleInternalCoordinatePtr>& simples,
    const Set<ConstInternalCoordinatePtr>& coords,
    bool validate
)
{
    clearFiles();

    //transform the gradients using intder
    int derlevel = 2;
    bool statpt = false;
    bool incxyz = false;
    IntderWriter iderwriter(mol, statpt, IntderWriter::TransformCartesians, derlevel, 0, simples, coords, incxyz);
    iderwriter.commit(DEFAULT_INTDER_INPUT_FILE);
    InputFileWriter::makeFile11(mol, xyzgrads);
    InputFileWriter::makeFile15(mol, xyzfc);

    stringstream sstr;
    sstr << INTDER_PATH << " < intder.inp > intder.out";
    int status = system(sstr.str().c_str());
    if (status != 0)
    {
        except("Intder returned a status error on force constant transformation to internals");
    }

    string file16 = getFileText("file16", ""); //pass blank comment, no comments

    string regexp = "([-]?\\d+[.]\\d+)";
    size_t nummatches;
    double* fcvals = get_regexp_double_array(regexp, file16, nummatches);
    int ncoords = coords.size();
    RectMatrixPtr intfc(ncoords, ncoords);
    intfc.assign(fcvals);
    delete[] fcvals;

    //make sure the xyzfc matches up to the original
    if (validate)
    {
        VectorPtr intgrads = transformGradientsToInternals(xyzgrads, mol, simples, coords, false);
        RectMatrixPtr check = transformForceConstantsToCartesian(intfc, intgrads, mol, simples, coords, false);

        if (!equals(check, xyzfc, 1e-5))
        {
            cout << "Internals transform failed" << endl;
            xyzfc.print("Correct FC");
            check.print("Incorrect FC");
            //regenerate original file
            RectMatrixPtr null = transformForceConstantsToInternals(xyzfc, xyzgrads, mol, simples, coords, false); //no validation
            abort();
        }
    }

    return intfc;
}

RectMatrixPtr
IntderWriter::transformGradientsToCartesian(
    ConstVectorPtr intgrads,
    const ConstMoleculePtr& mol,
    const Set<ConstSimpleInternalCoordinatePtr>& simples,
    const Set<ConstInternalCoordinatePtr>& coords,
    bool validate
)
{
    clearFiles();

    //transform the gradients using intder
    int derlevel = 1;
    bool statpt = false;
    bool incxyz = true;
    IntderWriter iderwriter(mol, statpt, IntderWriter::TransformFileInternals, derlevel, 0, simples, coords, incxyz);
    iderwriter.commit(DEFAULT_INTDER_INPUT_FILE);
    InputFileWriter::makeFile12(mol, coords, intgrads);

    stringstream sstr;
    sstr << INTDER_PATH << " < intder.inp > intder.out";
    int status = system(sstr.str().c_str());
    if (status != 0)
    {
        except("Intder returned a status error on gradient transformation to Cartesian");
    }

    string file11 = getFileText("file11", ""); //pass blank comment, no comments
    string regexp = "([-]?\\d+[.]\\d+)";
    size_t nummatches;
    double* gradvals = get_regexp_double_array(regexp, file11, nummatches);

    //some versions have different offsets so we compute the offset here
    int offset = nummatches - mol->natoms() * 3;
    double* gradarr = gradvals + offset;

    RectMatrixPtr gradients(mol->natoms(), 3); 
    gradients.assign(gradarr);
    delete[] gradvals;
    

    if (validate)
    {
        VectorPtr check = transformGradientsToInternals(gradients, mol, simples, coords, false);

        if (!equals(check, intgrads, 1e-8))
        {
            cout << "Cartesian gradient transform failed" << endl;
            intgrads.print("Correct");
            check.print("Incorrect");
            abort();
        }
    }

    return gradients;
}


RectMatrixPtr
IntderWriter::transformForceConstantsToCartesian(
    ConstRectMatrixPtr intfc,
    ConstVectorPtr intgrads,
    const ConstMoleculePtr& mol,
    const Set<ConstSimpleInternalCoordinatePtr>& simples,
    const Set<ConstInternalCoordinatePtr>& coords,
    bool validate
)
{
    clearFiles();

    //transform the gradients using intder
    int derlevel = 2;
    bool statpt = false;
    bool incxyz = true;
    IntderWriter iderwriter(mol, statpt, IntderWriter::TransformFileInternals, derlevel, 0, simples, coords, incxyz);
    iderwriter.commit(DEFAULT_INTDER_INPUT_FILE);
    InputFileWriter::makeFile12(mol, coords, intgrads);
    InputFileWriter::makeFile16(mol, intfc);

    stringstream sstr;
    sstr << INTDER_PATH << " < intder.inp > intder.out";
    int status = system(sstr.str().c_str());
    if (status != 0)
    {
        except("Intder returned a status error on force constant transformation to Cartesian");
    }

    string file15 = getFileText("file15", ""); //pass blank comment, no comments
    string regexp = "([-]?\\d+[.]\\d+)";
    size_t nummatches;
    double* fcvals = get_regexp_double_array(regexp, file15, nummatches);
    int nxyz = 3 * mol->natoms();
    RectMatrixPtr xyzfc(nxyz, nxyz);
    xyzfc.assign(fcvals);
    delete[] fcvals;

    //make sure the xyzfc matches up to the original
    if (validate)
    {
        RectMatrixPtr xyzgrads = transformGradientsToCartesian(intgrads, mol, simples, coords, false);
        RectMatrixPtr check = transformForceConstantsToInternals(xyzfc, xyzgrads, mol, simples, coords, false);

        if (!equals(check, intfc, 1e-5))
        {
            cout << "Cartesian transform failed" << endl;
            intfc.print("Correct FC");
            check.print("Incorrect FC");
            //rewrite original file
            //RectMatrixPtr null = transformForceConstantsToCartesian(intfc, intgrads, mol, simples, coords, false);
            abort();
        }
    }

    return xyzfc;
}

AnharmWriter::AnharmWriter(
    const ConstMoleculePtr& mol,
    const ConstForceFieldPtr& qff) 
 : InputFileWriter(mol), qff_(qff)
{
}

void
AnharmWriter::writeHeader()
{
    stringstream out;
out << "#*ANHARM*###### usage ######################" << endl;
out << "NIsotop    Cubic     ReOrder   Fermi1    Sigma" << endl;
out << "      NAtoms   Quartic   Coriolis  Fermi2   Print" << endl;
out << "[Reordering vector (N(I5))]" << endl;
out << "[Coriolis cutoff]" << endl;
out << "[Fermi1 cutoff]" << endl;
out << "[Fermi2 cutoff]" << endl;
out << "Atom1_x             Atom1_y            Atom1_z       Atom_mass" << endl;
out << "Atom2_x             Atom2_y            Atom2_z       Atom_mass" << endl;
out << ".......             .......            .......       ........." << endl;
out << "AtomN_x             AtomN_y            AtomN_z       Atom_mass" << endl;
out << "[Isotope1 mass1]" << endl;
out << "[Isotope1 mass2]" << endl;
out << "................" << endl;
out << "[Isotope1 massN]" << endl;
out << "[Isotope2 mass1]" << endl;
out << "................" << endl;
out << "[IsotopeN massN]" << endl;
out << "# ANHARM ###### The actual calculation ##############################" << endl;
    
    //just write 1 for the number of isotopes
    out << stream_printf("%5d", 1);
    //write the number of atoms
    out << stream_printf("%5d", mol_->natoms());
    add(out.str());
}

void
AnharmWriter::writeResonanceCutoffs()
{
    stringstream out;
    //file 20 is the cubic force field file, 24 the quartic
    out << stream_printf("%5d%5d", 20, 24);
    //don't reorder, i.e. 0
    out << stream_printf("%5d", 0);
    //include cutoffs for the three resonance cutoffs (i.e. make them all 1)
    out << stream_printf("%5d%5d%5d", 1, 1, 1);
    //don't know what these do, but put them in anway
    out << stream_printf("%5d   00", 0) << endl;
    //and print the resonance cutoffs, which will all have to be varied
    out << "20.0" << endl;
    out << "20.0" << endl;
    out << "20.0" << endl;

    add(out.str());
}

void
AnharmWriter::writeXYZ()
{
    stringstream out;
    //print geometry
    double unit_conversion = Molecule::getConversion("bohr");
    ConstRectMatrixPtr xyz = mol_->getXYZ(); //must be in bohr
    for (int n=0; n < mol_->natoms(); n++)
    {
        for (int x=0; x < XYZ_DIM; x++)
            out << stream_printf("%20.10f", xyz.get_element(n,x) * unit_conversion);
        out << stream_printf("%14.8f", mol_->getAtom(n)->mass()) << endl;
    }
    out << endl; //add a new line for good fortran measure
    add(out.str());
}

void
AnharmWriter::commit(const std::string& filename)
{
    writeHeader();
    writeResonanceCutoffs();
    writeXYZ();
    writeFooter();
    InputFileWriter::commit(filename);
}

RectMatrixPtr 
gigide::getMatrix(const std::string& text, mdim_t nrow, mdim_t ncol)
{
    //now that we have the geometry section, get the xyz coordinates
    string regexp = "([-]?\\d+[.]\\d+)\\s+([-]?\\d+[.]\\d+)\\s+([-]?\\d+[.]\\d+)";
    size_t length=0;
    double* values = get_regexp_double_array(regexp, text, length);

    //check consistency of the match and given dimensions
    if (length != nrow * ncol)
    {
        except("get_matrix:  Number of rows and columns inconsistent with number of values found");
    }

    RectMatrixPtr matrix(nrow, ncol);
    matrix.assign(values);

    delete[] values;

    return matrix;
}

RectMatrixPtr gigide::getXYZMatrix(vector<string>& atomlist,
                                   const std::string& geomsection,
                                   const std::string& inputunits,
                                   const std::string& outputunits)
{
    //get the atom list
    string regexp = "([a-zA-Z]{1,3})\\s+[-]?\\d+[.]\\d*"; 
	//the reg exp looks for all combinations of 3 letters followed by a space followed by a decimal number
	//e.g. He 3.25
    findmatch(atomlist, regexp, geomsection, FindAll | IncludeLineBreaks | UpperCase);

    //now get the xyz coordinates
    RectMatrixPtr xyzmat;
    try{
        xyzmat = getMatrix(geomsection, atomlist.size(), XYZ_DIM);
    }
    catch (GigideException e)
    {
        stringstream sstr;
        sstr << "XYZ Matrix not formmated properly: " << endl << geomsection << endl;
        except(sstr.str());
    }

    if (inputunits == angstrom and outputunits == bohr)
        xyzmat.scale(1.0/BOHR_TO_ANGSTROM);
    else if (inputunits == bohr and outputunits == angstrom)
        xyzmat.scale(BOHR_TO_ANGSTROM);

    return xyzmat;
}

VectorPtr gigide::getVector(const std::string& text)
{
    string regexp = "[-]?\\d+[.]\\d+";
    size_t length=0;
    double* values = get_regexp_double_array(regexp, text, length);

    //now get the xyz coordinates
    VectorPtr vec = new Vector(length);
    vec.assign(values);

    return vec;
}
