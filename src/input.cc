#include <src/input.h>
#include <src/exception.h>
#include <src/utilities.h>
#include <src/symmetry.h>
#include <config.h>

#define X 0
#define Y 1
#define Z 2

using namespace gigide;
using namespace smartptr;
using namespace std;

InputFile::InputFile(const std::string& text)
{
    filetext_ = text;
}

KeywordIteratorPtr
InputFile::getKeywords(const std::string& section_name) const
{
    string section = getSection(section_name);
    KeywordIteratorPtr keyset = new KeywordIterator(section);
    return keyset;
}

string
InputFile::getSection(const std::string& section_name) const
{
    if (!hasSection(section_name))
        return "";

    string regexp = section_name + "\n(.*?)[$]end";
    string section = get_regexp_string(regexp, filetext_, IncludeLineBreaks);
    return section;
}

bool
InputFile::hasSection(const std::string& section_name) const
{
    string regexp = section_name + "\n(.*?)[$]end";
    bool test = has_regexp_match(regexp, filetext_, IncludeLineBreaks);
    return test;
}

RectMatrixPtr 
InputFile::getGeometry(vector<string>& atomlist, const std::string& section_name) const
{
    string geomsection = getSection(section_name);

    //check if angstrom have been specified for the input geometry
    string ang_regexp = "(angstrom)";
    bool ang_units = has_regexp_match(ang_regexp, geomsection, LowerCase);
    string inputunits = bohr;
    if (ang_units) 
        inputunits = angstrom;

    //check if bohr or angstrom should be used for the displacements
    string outputunits = KeywordSet::getKeyword("bond units")->getValueString();

    return getXYZMatrix(atomlist, geomsection, inputunits, outputunits);
}

void 
InputFile::getSectionNames(std::vector<std::string>& names) const
{
    string regexp = "[$](.*?)\\n";
    findmatch(names, regexp, filetext_, FindAll);
}

GigideInputFile::GigideInputFile(const std::string& filetext)
    : InputFile(filetext)
{
}

KeywordSetPtr
GigideInputFile::getOptions() const
{
    string inputfile = GIGDATADIR "/defaults";
    string defaults = getFileText(inputfile, "#");
    string section = getSection("options");
    KeywordSetPtr keyset = new KeywordSet(section, defaults);
    return keyset;
}

void
GigideInputFile::validateInputFile()
{
    map<string,int> allowedsections;
    allowedsections["options"] = 1;
    allowedsections["symmetry"] = 1;
    allowedsections["internals"] = 1;
    allowedsections["geometry"] = 1;
    allowedsections["end"] = 1;
    allowedsections["symmops"] = 1;
    allowedsections["classes"] = 1;
    allowedsections["points"] = 1;
    allowedsections["axes"] = 1;
    allowedsections["subgroups"] = 1;
    allowedsections["degeneracy"] = 1;

    vector<string> matches;
    getSectionNames(matches);
    for (int i=0; i < matches.size(); i++)
    {
        string section = matches[i];
        if (!allowedsections[section])
        {
            allowedsections.erase(allowedsections.find(section));
            except(stream_printf("Unrecognized section %s", section.c_str()));
        }
    }
}

void
GigideInputFile::addSymmetryInternalCoordinates(
    const ConstMoleculePtr& mol,
    vector<string>& simples,
    map<string, SimpleInternalCoordinatePtr>& label_map,
    vector<InternalCoordinatePtr>& coords
) 
{
    vector<SimpleInternalCoordinatePtr> blksimples;
    for (int c=0; c < simples.size(); c++)
    {
        SimpleInternalCoordinatePtr coord = label_map[simples[c]];
        if (!coord)
        {
            except(stream_printf("%s is not a valid name in $symmetry", simples[c].c_str()));
        }
        blksimples.push_back(coord);
    }

    //form delocalized internals on just the block
    int oldsize = coords.size();
    vector<SymmetryInternalCoordinatePtr> delocals;
    VectorPtr evals = InternalCoordinate::addDelocalizedInternalCoordinates(mol, blksimples, delocals);
    vector<CoordinateSubspacePtr> subspaces;
    mol->getPointGroup()->formBasis(delocals, subspaces);
    for (int s=0; s < subspaces.size(); ++s)
    {
        for (int n=0; n < subspaces[s]->order(); ++n)
        {
            InternalCoordinatePtr next(subspaces[s]->getCoordinate(n));
            coords.push_back(next);
        }
    }
    int blksize = coords.size() - oldsize;
    cout << stream_printf("Formed %d symmetry internals from %d simples", blksize, blksimples.size()) << endl;
}

void
GigideInputFile::appendCoordinates(
    const std::string& name,
    const KeywordValuePtr& keyval,
    const ConstMoleculePtr& mol,
    vector<InternalCoordinatePtr>& coords,
    map<string, SimpleInternalCoordinatePtr>& label_map
)
{
    string symm_regexp = "set\\d";
    bool symmint = has_regexp_match(symm_regexp, name, LowerCase);
    if (symmint) 
    {
        vector<string> simples; keyval->getValueVectorString(simples);
        addSymmetryInternalCoordinates(mol, simples, label_map, coords);
        return;
    }
    else; //standard read

    //we need a matrix kit for build coefficient vectors
    string label_regexp = "([a-zA-Z][a-zA-Z\\d]*)";
    //to construct a symmetry internal coordinate, we will need a vector of coefficients
    VectorPtr coeffs;
    //and a set of simple internal coordinates
    vector< SimpleInternalCoordinatePtr > simples_in_symm;

    string coord_text = keyval->getValueString();
    vector<string> simple_matches; findmatch(simple_matches, label_regexp, coord_text, FindAll);
    //now that we know how many simples are involved, build an empty vector of coefficients
    coeffs = new Vector(simple_matches.size()); 
    for (int n=0; n < simple_matches.size(); ++n)
    {
        string simple_name = simple_matches[n];
        string coeff_regexp = "([\\-+]?\\s*\\d+[.]?\\d*)\\s+" ; coeff_regexp += simple_name;
        string sign_regexp = "([\\-+])\\s*" ; sign_regexp += simple_name;
        bool has_coeff = has_regexp_match(coeff_regexp, coord_text);
        bool has_sign = has_regexp_match(sign_regexp, coord_text);
        double coeff;
        if (has_coeff)
        {
            coeff = get_regexp_double(coeff_regexp, coord_text);
        }
        else if (has_sign) //the coefficient is either plus or minus 1
        {
            string sign = get_regexp_string(sign_regexp, coord_text);
            if      (sign == "-") coeff = -1.0;
            else if (sign == "+") coeff = 1.0;
        }
        else //must be the first coordinate it must +1
        {
            coeff = 1.0;
        }
        //now that we have figured out the coefficient, set the element in the vector
        coeffs.set_element(n, coeff);
        //figure out from the string map which simple internal coordinate this corresponds to
        string simple_label = simple_matches[n];
        SimpleInternalCoordinatePtr simple = label_map[simple_label];
        if (!simple)
        {
            except(stream_printf("Unrecognized coordinate label %s", simple_label.c_str()));
        }
        //and add it the set of simples for this symmetry internal coordinates
        simples_in_symm.push_back(simple);
    }
    //and finally... build the symmetry internal coordinate
    //but first normalize the coefficients for intder
    //coeffs.normalize();
    SymmetryInternalCoordinatePtr symm = new SymmetryInternalCoordinate(coeffs, simples_in_symm, mol);
    coords.push_back(symm);
}

void
GigideInputFile::readSymmetryCoordinates(
    KeywordIteratorPtr keyset,
    map<string, SimpleInternalCoordinatePtr>& label_map,
    vector<InternalCoordinatePtr>& symm_coords,
    const ConstMoleculePtr& mol
)
{
    //start the iterator
    for (keyset->start(); !keyset->finished(); keyset->next())
    {
        string name = keyset->getName();
        KeywordValuePtr keyval = keyset->getKeyword();
        appendCoordinates(name, keyval, mol, symm_coords, label_map);
    }

}

RectMatrixPtr
GigideInputFile::getGeometry(
    vector<string>& atomlist
)
{
    return InputFile::getGeometry(atomlist, "geometry");
}

void
GigideInputFile::getSimpleCoordinates(
    MoleculePtr mol,
    vector< SimpleInternalCoordinatePtr >& simples
)
{
    //for each line in the section, build a new coordinate
    KeywordIteratorPtr keyset = getKeywords("internals");
    for (keyset->start(); !keyset->finished(); keyset->next())
    {
        KeywordValuePtr keyval = keyset->getKeyword();
        string coord_type = keyval->popString();
        string name = keyset->getName();
    
        SimpleInternalCoordinatePtr new_coord;
        if      (coord_type == "bond" || coord_type == "stretch" || coord_type == "stre")
        {
            vector<int> connect = keyval->popVectorInteger(2);
            new_coord = new BondLength(connect, mol);
        }
        else if (coord_type == "angle" || coord_type == "bend")
        {
            vector<int> connect = keyval->popVectorInteger(3);
            new_coord = new BondAngle(connect, mol);
        }
        else if (coord_type == "tors" || coord_type == "torsion" || coord_type == "dih" || coord_type == "dihedral")
        {
            vector<int> connect = keyval->popVectorInteger(4);
            new_coord = new Torsion(connect, mol);
        }
        else if (coord_type == "oop" || coord_type == "out" || coord_type == "oopbend")
        {
            vector<int> connect = keyval->popVectorInteger(4);
            new_coord = new OutOfPlaneBend(connect, mol);
        }
        else if (coord_type == "linx")
        {
            vector<int> connect = keyval->popVectorInteger(4);
            new_coord = new LinX(connect, mol);
        }
        else if (coord_type == "liny")
        {
            vector<int> connect = keyval->popVectorInteger(4);
            new_coord = new LinY(connect, mol);
        }
        else if (coord_type == "lin1")
        {
            vector<int> connect = keyval->popVectorInteger(3);
            vector<double> linvec = keyval->popVectorDouble(3);
            new_coord = new Lin1(connect, mol, linvec);
        }
        else
        {
            except(stream_printf("Unrecognized coordinate type %s", coord_type.c_str()));
        }
        label_map_[name] = new_coord;
        simples.push_back(new_coord);
    }
}

void
GigideInputFile::getSymmetryCoordinates(
    const ConstMoleculePtr& mol,
    vector<SimpleInternalCoordinatePtr>& simples,
    vector<InternalCoordinatePtr>& final_coords
)
{
    if (hasSection("symmetry"))
    {
        readSymmetryCoordinates(getKeywords("symmetry"), label_map_, final_coords, mol);
    }
    else
    {
        VectorPtr one(1); one.set_element(0, 1.0);
        for (int i=0; i < simples.size(); ++i)  //build symmetry coordinates that are entirely the single coordinate
        {
            vector<SimpleInternalCoordinatePtr> one_simple; one_simple.push_back(simples[i]);
            SymmetryInternalCoordinatePtr coord(new SymmetryInternalCoordinate(one, one_simple, mol));
            final_coords.push_back(coord);
        }
    }

    //read in symmetry relations
    KeywordIteratorPtr keyset = getKeywords("degeneracy");
    for (keyset->start(); !keyset->finished(); keyset->next())
    {
        KeywordValuePtr keyval = keyset->getKeyword();
        string name = keyset->getName();
        int i = keyval->popInteger() - 1; //1 based to zero based
        int j = keyval->popInteger() - 1;
        
        final_coords[i]->addDegenerateCoordinate(final_coords[j]);
        final_coords[j]->addDegenerateCoordinate(final_coords[i]);
    }
}


void
GigideInputFile::readXYZPoints(
    const ConstMoleculePtr& mol,
    map<string, XYZPointPtr >& points
)
{
    //you get all individual atoms for free
    for (int i=0; i < mol->natoms(); ++i)
    {
        stringstream sstr;
        sstr << (i+1);
        ConstAtomPtr atom = mol->getAtom(i);
        XYZPointPtr pt = new XYZPoint(atom);
        points[sstr.str()] = pt;
    }
    //get origin for free
    ConstVectorPtr zero = new Vector(3);
    XYZPointPtr origin = new XYZPoint(zero);
    points["0"] = origin;

    if (!hasSection("points"))
        return;

    KeywordIteratorPtr iter = getKeywords("points");
    for (iter->start(); !iter->finished(); iter->next())
    {
        KeywordValuePtr keyval = iter->getKeyword();
        string name = iter->getName();
        string type = keyval->popString();

        if      (type == "point") //read a vector of xyz coordinates
        {
            vector<double> xyz = keyval->popVectorDouble(3);
            VectorPtr xyzvec = new Vector(3);
            for (int i=0; i < 3; ++i)
                xyzvec.set_element(i, xyz[i]);
            XYZPointPtr pt = new XYZPoint( ConstVectorPtr(xyzvec) ); //cast so compiler doesn't freak
            points[name] = pt;
        }
        else if (type == "midpoint")
        {
            vector<int> atomlist = keyval->getValueVectorInteger();
            vector<ConstAtomPtr > atoms;
            for (int i=0; i < atomlist.size(); ++i)
                atoms.push_back(mol->getAtom(atomlist[i] - 1));
            XYZPointPtr pt = new Midpoint(atoms);
            points[name] = pt;
        }
        else
        {
            except(stream_printf("Unrecognized point type %s", type.c_str()));
        }

    }
}

void
GigideInputFile::geometryError(
    const std::string& value,
    const std::string& spec
)
{
    except(stream_printf("Invalid geometry specifier %s for %s", spec.c_str(), value.c_str()));
}

void
GigideInputFile::readAxes(
    const ConstMoleculePtr& mol,
    map<string, XYZPointPtr >& points,
    map<string, AxisPtr >& axes
)
{
    //you get all bonds for free
    for (int i=0; i < mol->natoms(); ++i)
    {
        ConstAtomPtr ati = mol->getAtom(i);
        XYZPointPtr pti = new XYZPoint(ati);
        for (int j=0; j < mol->natoms(); ++j)
        {
            if (i==j)
                continue;

            ConstAtomPtr atj = mol->getAtom(j);
            XYZPointPtr ptj = new XYZPoint(atj);
            AxisPtr axis = new Axis(pti, ptj);
            stringstream sstr;
            sstr << "r" << (i+1) << (j+1);
            axes[sstr.str()] = axis;
        }
    }
    //get xyz axes for free
    VectorPtr vec = new Vector(3);
    vec.zero(); vec.set_element(X, 1.0); axes["x"] = new Axis(vec);
    vec.zero(); vec.set_element(Y, 1.0); axes["y"] = new Axis(vec);
    vec.zero(); vec.set_element(Z, 1.0); axes["z"] = new Axis(vec);

    //get origin to all atoms for free
    XYZPointPtr origin = points["0"];
    for (int i=0; i < mol->natoms(); ++i)
    {
        stringstream sstr;
        sstr << (i+1);
        XYZPointPtr pti = points[sstr.str()];

        stringstream sstr_axis;
        sstr_axis << "r" << 0 << (i+1);
        AxisPtr axis = new Axis(origin, pti);
        axes[sstr_axis.str()] = axis;
    }

    if (!hasSection("axes"))
        return;

    KeywordIteratorPtr iter = getKeywords("axes");
    for (iter->start(); !iter->finished(); iter->next())
    {
        KeywordValuePtr keyval = iter->getKeyword();
        string name = iter->getName();
        string type = keyval->popString();

        if      (type == "axis") //read 2 points
        {
            XYZPointPtr pt1 = points[keyval->popString()];
            XYZPointPtr pt2 = points[keyval->popString()];
            AxisPtr axis = new Axis(pt1, pt2);
            axes[name] = axis;
        }
        else if (type == "bisect")
        {
            string name1 = keyval->popString();
            AxisPtr ax1 = axes[name1];
            if (!ax1) geometryError(name, name1);
                
            string name2 = keyval->popString();
            AxisPtr ax2 = axes[name2];
            if (!ax2) geometryError(name, name2);

            VectorPtr vec = ax2->getVector() + ax1->getVector();
            AxisPtr axis = new Axis(vec);
            axes[name] = axis;
        }
        else
        {
            except(stream_printf("Unrecognized axis type %s", type.c_str()));
        }

    }
}

void
GigideInputFile::readSymmetryOperations(
    const ConstMoleculePtr& mol,
    map<std::string, XYZPointPtr >& points,
    map<std::string, AxisPtr >& axes,
    map<string, SymmetryOperationPtr >& symmops
)
{
    //strip const modifier
    MoleculePtr unconstmol = boost::const_pointer_cast<Molecule, const Molecule>(mol);

    KeywordIteratorPtr iter = getKeywords("symmops");
    for (iter->start(); !iter->finished(); iter->next())
    {
        KeywordValuePtr keyval = iter->getKeyword();
        string name(iter->getName());
        string type(keyval->popString());

        SymmetryOperationPtr op;
        if      (type == "c")
        {
            int order = keyval->popInteger();
            string axname = keyval->popString();
            AxisPtr axis = axes[axname];
            if (!axis) geometryError(name, axname);
            
            op = new Rotation(axis->getVector(), order);
        }
        else if (type == "s")
        {
            int order = keyval->popInteger();
            string axname = keyval->popString();
            AxisPtr axis(axes[axname]);
            if (!axis) geometryError(name, axname);
            
            op = new ImproperRotation(axis->getVector(), order);
        }
        else if (type == "sigma")
        {
            int count = keyval->count();
            if (count == 1) //axis specified
            {
                string axname = keyval->popString();
                AxisPtr axis(axes[axname]);
                if (!axis) geometryError(name, axname);

                op = new Reflection(axis->getVector());
            }
            else if (count == 3) //3 points specified
            {
                string ptname;
                ptname = keyval->popString();
                XYZPointPtr pt1(points[ptname]);
                if (!pt1) geometryError(name, ptname);

                ptname = keyval->popString();
                XYZPointPtr pt2(points[ptname]);
                if (!pt2) geometryError(name, ptname);

                ptname = keyval->popString();
                XYZPointPtr pt3(points[ptname]);
                if (!pt3) geometryError(name, ptname);

                op = new Reflection(pt1, pt2, pt3);
            }
            else
            {
                except("Symmetry plane requires either 3 points or an axis");
            }
        }
        else
        {
            except(stream_printf("Unrecognized symmetry operation %s", type.data()));
        }

        if (!op)
        {
            cerr << "invalid specification of " << type << endl;
        }

        symmops[name] = op;
        unconstmol->addSymmetryOperation(op);
    }
}

void
GigideInputFile::readSubgroups(
    const ConstMoleculePtr& mol,
    map<string, SymmetryOperationPtr>& symmops
)
{
    MoleculePtr unconstmol = boost::const_pointer_cast<Molecule, const Molecule>(mol);

    KeywordIteratorPtr iter = getKeywords("subgroups");
    for (iter->start(); !iter->finished(); iter->next())
    {

        KeywordValuePtr keyval = iter->getKeyword();
        string name = iter->getName();

        PointGroupPtr subgroup(new PointGroup(mol, name));

        vector<string> subops; keyval->getValueVectorString(subops);
        vector<string>::const_iterator itvec;
        map<string, SymmetryOperationPtr>::const_iterator itmap;
        for (itvec = subops.begin(); itvec != subops.end(); ++itvec)
        {
            string opname = *itvec;
            itmap = symmops.find(opname);
            if (itmap == symmops.end())
                except(stream_printf("%s is not a valid symmetry operation name for subgroup %s", opname.c_str(), name.c_str()));
            
            subgroup->addOperation(itmap->second);
        }
        unconstmol->getPointGroup()->addSubgroup(subgroup);
    }
}

void
GigideInputFile::readSymmetryOperations(
    const ConstMoleculePtr& mol
)
{
    map<string, XYZPointPtr > points;
    readXYZPoints(mol,points);

    map<string, AxisPtr > axes;
    readAxes(mol,points,axes);

    map<string, SymmetryOperationPtr> symmops;
    map<string, SymmetryOperationPtr>::iterator it;
    readSymmetryOperations(mol, points, axes, symmops);
    readSubgroups(mol, symmops);
}
