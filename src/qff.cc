#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <fcntl.h>
#include <errno.h>
#include <config.h>
#include <src/fileio.h>
#include <src/fit.h>
#include <src/defines.h>
#include <src/pyregexp.h>
#include <src/exception.h>
#include <src/symmetry.h>
#include <src/utilities.h>
#include <src/qff.h>
#include <src/derivative.h>
#include <src/displacement.h>
#include <src/permutation.h>
#include <src/units.h>

#define USE_SVD 0

#define ZERODISP 1e-3


#define dbg cout

#define X 0
#define Y 1
#define Z 2

using namespace gigide;
using namespace smartptr;
using namespace std;

SerialDeclare(ForceField);

ForceField::ForceField(
   const ConstMoleculePtr& center_mol,
   const Set<ConstInternalCoordinatePtr>& coords,
   const Set<ConstSimpleInternalCoordinatePtr>& simples
)
    : coords_(coords), simples_(simples), mol_(center_mol), energy_(0), 
    disp_sizes_(GigideKeyword::getDisplacementSizes(coords.size()))
{
    SetRuntime(ForceField);

    //the number of displacement sizes must match number of coordinates
    if (coords.size() != disp_sizes_.size())
    {
        except(stream_printf("Have %d coordinates, but %d disp sizes",
                                  coords.size(), disp_sizes_.size()));
    }

    //make sure equivalent displacements have the same displacement sizes

    int robustlevel = KeywordSet::getKeyword("robustness")->getValueInteger() - 1;

    //the final level of derivative wanted
    nderiv_ = KeywordSet::getKeyword("nderiv")->getValueInteger();
    //the level of derivative the user will provide
    nvalue_ = KeywordSet::getKeyword("nvalue")->getValueInteger();
    int disp_level = nderiv_ - nvalue_;
    disp_iter_ = new DisplacementIterator(coords, simples, center_mol, disp_sizes_);
    value_iter_ = new DerivativeIterator(nvalue_, coords, center_mol);
    compute_iter_ = new DerivativeIterator(disp_level, coords, center_mol);
    deriv_iter_ = new DerivativeIterator(nderiv_, coords, center_mol);
    fit_iter_ = new DerivativeIterator(disp_level + robustlevel, coords, center_mol);

    if (nvalue_ == 0)
    {
        DerivativeIterator::iterator it;
        for (it = compute_iter_->begin(); it != compute_iter_->end(); ++it)
            (*it)->mapEquivalentDerivatives(compute_iter_);
        for (it = deriv_iter_->begin(); it != deriv_iter_->end(); ++it)
            (*it)->mapEquivalentDerivatives(deriv_iter_);
    }
    
    vector<int> zeroes(coords.size(), 0);
    vector<DisplacementPtr> uniqcoords;
    map<DisplacementPtr, double> dispsizechecks;
    for (int i=0; i < coords.size(); i++)
    {
        ConstInternalCoordinatePtr coord = coords[i];
        for (int j=i+1; j < coords.size(); ++j)
        {
            ConstInternalCoordinatePtr test = coords[j];

            if (!coord->isEquivalent(test))
                continue;

            //found equivalent.  make sure they have the same displacement size
            double diff = fabs(disp_sizes_[j] - disp_sizes_[i]);
            if (diff > 1e-12)
            {
                except(stream_printf("Equivalent coordinates %d and %d have different displacement sizes. Please correct this.", i+1,j+1));
            }
        }
    }

    /** At this point, we may or may not be able to use symmetry. If the user is providing energies,
        then the quantity being differentiated is totally symmetric, so this is easy.  If not,
        then the quantity being differentiated may or may not be totally symmetric, which makes
        things much more difficult.  For now, if nvalue is not zero, don't make any use of symmetry */
    if (nvalue_)
        include_zeros_ = true;
    else
        include_zeros_ = false;

    //now that we have the derivative and displacement iterators, build the fitting matrix
    buildFit();

    disp_iter_->mapEquivalentDisplacements();
}

ForceField::ForceField(const ArchivePtr& arch)
    : Serializable(arch)
{
    SetRuntime(ForceField);

    serial_load(include_zeros);
    serial_load(disp_sizes);
    serial_load(disp_iter);
    serial_load(fit_iter);
    serial_load(value_iter);
    serial_load(compute_iter);
    serial_load(deriv_iter);
    serial_load(mol);
    serial_load(fc);
    serial_load(grads);
    serial_load(simples);
    serial_load(coords);
    serial_load(energy);
    serial_load(nvalue);
    serial_load(nderiv);
}

void
ForceField::serialize(const ArchivePtr& arch) const
{
    Serializable::serialize(arch);

    serial_save(include_zeros);
    serial_save(disp_sizes);
    serial_save(disp_iter);
    serial_save(fit_iter);
    serial_save(value_iter);
    serial_save(compute_iter);
    serial_save(deriv_iter);
    serial_save(mol);
    serial_save(fc);
    serial_save(grads);
    serial_save(simples);
    serial_save(coords);
    serial_save(energy);
    serial_save(nvalue);
    serial_save(nderiv);
}

void
ForceField::validateFit()
{
    DerivativeIterator::iterator it;
    for (it = compute_iter_->begin(); it != compute_iter_->end(); ++it)
    {
        DerivativePtr deriv = *it;
        if (deriv->zero()) 
        {
            cout << "Rigorously zero derivative: "; deriv->print();
            continue;
        }

        int robustlevel = robustness(deriv->level());
        if (deriv->nonzero() && deriv->isUnique())
            deriv->validateFit(robustlevel);
    }
}

void
ForceField::validateRow(
    vector<TaylorTermPtr>& terms,
    ConstVectorPtr coefs
)
{
    vector<TaylorTermPtr >::iterator it;
    for (it = terms.begin(); it != terms.end(); ++it)
        (*it)->reset();

    int i = 0;

    DisplacementIterator::iterator itdisp;
    for (itdisp = disp_iter_->begin(); itdisp != disp_iter_->end(); ++itdisp, ++i)
    {
        double coef = coefs.get_element(i);
        for (it = terms.begin(); it != terms.end(); ++it)
            (*it)->accumulate(*itdisp, coef);
    }

    for (it = terms.begin(); it != terms.end(); ++it)
    {
        double coef = (*it)->coef();
        if ( fabs(coef) > 0.02 )
            (*it)->print();
    }
}

RectMatrixPtr
ForceField::buildB() const
{
    return InternalCoordinate::formBMatrix(coords_);
}

int
ForceField::ncoords() const
{
    return coords_.size();
}

int
ForceField::robustness(int derivlevel) const
{
    int robustness = KeywordSet::getKeyword("robustness")->getValueInteger();
    int nderiv_comp = derivlevel_compute();
    //assume a 2n-1 rule for robustness
    int extra_robustness = (nderiv_comp - derivlevel) / 2;
    int final_robustness = robustness + extra_robustness;
    return final_robustness;
}

int
ForceField::derivlevel_compute() const
{
	int nvalue = KeywordSet::getKeyword("nvalue")->getValueInteger();
    int nder = KeywordSet::getKeyword("nderiv")->getValueInteger();
    int nderiv_comp = nder - nvalue;

    return nderiv_comp;
}

void
ForceField::buildNumericalFit()
{
	int debug = KeywordSet::getKeyword("fit debug")->getValueInteger();
    int nderiv_comp = compute_iter_->level();

    //figure out the dimensions of the fitting matrix
    int ndisps = disp_iter_->ndisps();
    int nderivs_fit = fit_iter_->nderivs();
    //this will be the matrix that generates the "Taylor" series
    RectMatrixPtr TMat(ndisps, nderivs_fit); 

    double* scratch = new double[TMat.nrow() * TMat.ncol()];
    double* scratchptr = scratch;
    DisplacementIterator::iterator itdisp;
    DerivativeIterator::iterator itder;
    for (itdisp = disp_iter_->begin(); itdisp != disp_iter_->end(); ++itdisp)
    {
        DisplacementPtr disp = *itdisp;
        vector<double> increments = disp->getIncrements();
        for (itder = fit_iter_->begin(); itder != fit_iter_->end(); ++itder, ++scratchptr)
        {
            DerivativePtr deriv = *itder;
            double coeff = deriv->taylorCoeff(increments);
            (*scratchptr) = coeff;
        }
    }
    TMat.assign(scratch);
    delete[] scratch;

    SymmMatrixPtr lsqMat(TMat.ncol());
    RectMatrixPtr XT = TMat.t();
    lsqMat.accumulate_symmetric_product(XT);
    SymmMatrixPtr inv = lsqMat.i();
    RectMatrixPtr fit_mat = inv * XT;

    int n = 0;

    for (itder = compute_iter_->begin(); itder != compute_iter_->end(); ++itder, ++n)
    {
        DerivativePtr deriv = *itder;
        VectorPtr coefs = fit_mat.get_row(n);
        FitPtr fit = new Fit(disp_iter_, coefs);
        deriv->assignFit(fit);
    }
}

void
ForceField::buildFormulaFit()
{
    int nvalue = KeywordSet::getKeyword("nvalue")->getValueInteger();
    bool usezero = nvalue == 0; 
    DerivativeIterator::iterator itder;

    for (itder = compute_iter_->begin(); itder != compute_iter_->end(); ++itder)
    {
        DerivativePtr deriv = *itder;
        if (deriv->zero() && usezero)
        {
            deriv->assignFit(0); //just assign an empty fit, this deriv is zero
        }
        else if (!deriv->isUnique())
        {
            deriv->assignFit(0); //just assign an empty fit, this deriv is zero
        }
        else
        {
            int robustlevel = robustness(deriv->level());
            Fit::assignFit(deriv, disp_iter_, robustlevel); //just use 2 for now
        }
    }
}

void
ForceField::buildFit()
{
    string fittype = KeywordSet::getKeyword("fittype")->getValueString();
    if (fittype == "formula")
        buildFormulaFit();
    else if (fittype == "numerical")
        buildNumericalFit();
    else
    {
        except(stream_printf("Invalid fit type %s. Choose formula or numerical.", fittype.c_str()));
    }
}

void
ForceField::computeDerivatives(DerivativePtr value)
{
    int debug = KeywordSet::getKeyword("fit debug")->getValueInteger();
    int index = 0;

    DerivativeIterator::iterator itder;
    for (itder = compute_iter_->begin(); itder != compute_iter_->end(); ++itder)
    {
        DerivativePtr deriv = *itder;
        DerivativePtr total_deriv = deriv_iter_->findCombination(deriv, value);
        if (total_deriv.get() == NULL)
        {
            stringstream sstr;
            sstr << "Could not find combination derivative" << endl;
            deriv->print(sstr);
            value->print(sstr);
            except(sstr.str());
        }

        if (total_deriv->nonzero() && total_deriv->isUnique())
        {
            double denom = deriv->displacementCoeff(disp_sizes_);
            double derval = deriv->compute(value) / denom;
            if (debug)
            {
                cout << "Assign value " << derval << " to derivative ";
                total_deriv->print();
            }
            total_deriv->setValue(derval);
        }
        else if (total_deriv->zero())//set rigorously to zero
        {
            total_deriv->setValue(0.0); 
        }
    }
}

SymmMatrixPtr
ForceField::formGMatrix(
    const Set<ConstInternalCoordinatePtr>& coords
)
{
    //now we must build the G matrix
    RectMatrixPtr B = InternalCoordinate::formBMatrix(coords);
    ConstMoleculePtr mol = coords[0]->mol();
    int natoms = mol->natoms();
    int ncoords = coords.size();
    SymmMatrixPtr G(ncoords);
    for (int i=0; i < ncoords; i++)
    {
        VectorPtr Bi = B.get_row(i);
        for (int j=i; j < ncoords; j++)
        {
            VectorPtr Bj = B.get_row(j);
            double Gij = 0.0;
            for (int k=0; k < XYZ_DIM * natoms; k++)
            {
                double mass = mol->getAtom(k/3)->mass();
                Gij += Bi[k] * Bj[k] / mass;
            }
            G.set_element(i, j, Gij);
        }
    }
    return G;
}

void
ForceField::buildArrays()
{
    //first set up the harmonic force constant arrays
    fc_ = SymmMatrixPtr(ncoords());
    grads_ = VectorPtr(ncoords()); 

    //iterate through the derivative
    int derivnum = 0;
    DerivativeIterator::iterator itder;
    for (itder = deriv_iter_->begin(); itder != deriv_iter_->end(); ++itder)
    {
        DerivativePtr deriv = *itder;
        int level = deriv->level();
        if (level==1) //gradient
        {
            vector<int> indices = deriv->indices();
            int i = indices[0];
            grads_.set_element(i, deriv->value());
        }
        else if (level==2) //force constant
        {

            vector<int> indices = deriv->indices();
            int i = indices[0]; int j = indices[1];
            fc_.set_element(i, j, deriv->value());
        }
    }

    G_ = formGMatrix(coords_);
    grads_.print("Gradients");
}

void
ForceField::computeNormalModes(
    ConstSymmMatrixPtr F,
    ConstSymmMatrixPtr G,
    VectorPtr& freqs,
    RectMatrixPtr& evecs
)
{
    VectorPtr gvals;
    RectMatrixPtr gvecs;
    G.eigen(gvals, gvecs);
    gvals.sort();
    //compute G^(-1/2)
    double tol = 1e-8;
    SymmMatrixPtr G_m12 = G.invsqrt_matrix(); 
    SymmMatrixPtr G_12 = G.sqrt_matrix();

    SymmMatrixPtr GFG = G.clone();
    GFG.accumulate_transform(G_12, F);

    VectorPtr evals;
    //get the eigenvalues
    GFG.eigen(evals, evecs);
    //transform the eigenvectors back to the original basis
    evecs = G_12 * evecs;

    evals.sort();

    //if we are in attajoule units, convert back to hartree
    if (KeywordSet::getKeyword("energy units")->getValueString() == "aj")
        evals.scale(1.0/HARTREE_TO_AJ);

    //and convert to mks units
    evals.scale(AU_TO_MKS);

    //take the square root since the eigenvalue is actually the angular frequency squared
    freqs = evals.clone();
    for (int i=0; i < evals.n(); ++i)
    {
        double val = evals[i];
        if (val < 0)
            freqs.set_element(i, -1.0 * sqrt(-val));
        else
            freqs.set_element(i, sqrt(val));
    }

    //convert from hertz to wavenumber
    freqs.scale(HERTZ_TO_WN);

    //finally, if we are in angstrom units, convert back to bohr
    if (KeywordSet::getKeyword("bond units")->getValueString() != "bohr")
        freqs.scale(BOHR_TO_ANGSTROM);
}

void
ForceField::computeFrequencies()
{
    RectMatrixPtr L;
    VectorPtr freqs;
    computeNormalModes(fc_, G_, freqs, L);

    freqs.print("Frequencies");
}

vector<double>
ForceField::getDisps(vector<int> disp_numbers)
{
    vector<double> disps;
    for (int i=0; i < disp_numbers.size(); i++)
        disps.push_back( disp_sizes_[i] * disp_numbers[i] );
    return disps;
}

void
ForceField::writeMoleculeToDispcart(
    ofstream& dispcart,
    const ConstMoleculePtr& mol,
    int dispnum
) const
{
    int indent = 1;
    string xyz_string = mol->getXYZString(bohr, indent);
    dispcart << stream_printf("%d", dispnum);

    vector<string> details; KeywordSet::getKeyword("mol details")->getValueVectorString(details);
    vector<string>::iterator it;
    string molinfo;
    for (it = details.begin(); it != details.end(); ++it)
    {
        molinfo = *it;
        if    (molinfo == "pg") dispcart << " " << mol->getPointGroup()->name();
    }
    dispcart << endl;

    dispcart << xyz_string << endl;
}

void
ForceField::writeDisplacementsToFile(string filename) const
{
    ofstream dispcart;
    dispcart.open(filename.c_str());

    ConstDisplacementPtr disp;
    int dispnum = 1;
    for (DisplacementIterator::const_iterator it(disp_iter_->begin()); it != disp_iter_->end(); ++it)
    {
        disp = *it;
        if (disp->isUnique())
        {
            ConstMoleculePtr mol = disp->getDisplacementMolecule();
            writeMoleculeToDispcart(dispcart, mol, dispnum);
            ++dispnum;
        }
    }


    //if we are computing things by displacement of higher derivatives, we should always include the zero point
    //that way, we don't compute exact derivatives by finite differences, but rather get them exactly from the center point
    int nvalue = KeywordSet::getKeyword("nvalue")->getValueInteger();
    if (nvalue && !disp_iter_->hasZeroDisplacement())
    {
        cout << "Adding center molecule for exact values" << endl;
        writeMoleculeToDispcart(dispcart, mol_, dispnum);
    }

    dispcart.close();
}


void
ForceField::generateDisplacements()
{
    //get to generating!

    cout << stream_printf("Generating %d unique displacements from %d total displacements", 
                          disp_iter_->nunique(),
                          disp_iter_->ndisps()
                          ) << endl;

    int nuniq = disp_iter_->nunique();
    if (nuniq << 1 == GGF1 || nuniq << 1 == GGF2)
        cout << endl << "Jackpot.  Gigide." << endl;

    int dispnum = 1;
    DisplacementIterator::iterator it;
    for (it = disp_iter_->begin(); it != disp_iter_->end(); ++it)
    {
        DisplacementPtr disp = *it;
        if (disp->isUnique())
        {
            vector<double> increments = disp->getIncrements();
            cout << "Generating "; disp->print();
            disp->generateDisplacement();
            cout << " ...generated displacement " << dispnum << endl;
            ++dispnum;

            disp->validate();
        }
    }
}

bool
ForceField::getXMLEnergy(const XMLParser& xml, double& energy)
{
    /** Get the energy */
    XMLParser enode = xml->getChildByAttribute("energy", "type", "molecular");
    if (enode.get() == NULL)
        return false;

    if (!enode->hasAttribute("value")) return false;
    string etext = enode->getAttribute("value");
    stringstream sstr(etext);
    sstr >> energy;

    return true; //energy found
}

bool
ForceField::getXMLXYZ(const XMLParser& xml, RectMatrixPtr& geom)
{
    /** Get the xyz coordinates */
    XMLParser xyznode = xml->getChildElement("xyz");
    if (xyznode.get() == NULL)
        return false;

    string xyztext = xyznode->getText();

    vector<string> atoms;
    //we want specific units
    string outunits = KeywordSet::getKeyword("bond units")->getValueString();
    //the input units are not necessarily given
    string inunits = "bohr";
    if (xyznode->hasAttribute("units"))
        inunits = xyznode->getAttribute("units");

    geom = getXYZMatrix(atoms, xyztext, inunits, outunits);
    return true;
}

bool
ForceField::getXMLGradients(const XMLParser& xml, RectMatrixPtr& gradients)
{
    XMLParser gradnode = xml->fetchNode("gradient");
    if (gradnode.get() == NULL)
        return false;

    string gradtype = "xyz";
    if (gradnode->hasAttribute("type"))     
        gradtype = gradnode->getAttribute("type");
    if (gradtype == "internals")
    {
        except("Cannot read in internals gradients. Please give only xyz gradients.");
    }

    if (!HAS_INTDER)
    {
        except("No intder path given. Cannot transform coordinates.");
    }

    //if not internals specified, do xyz 
    string gradunits = "hartree/bohr";
    if (gradnode->hasAttribute("units"))
        gradunits = gradnode->getAttribute("units");
    if (gradunits != "hartree/bohr")
    {
        except("Please input gradients in XML file in hartree/bohr units");
    }

    string gradtext = gradnode->getText();
    gradients = getMatrix(gradtext, mol_->natoms(), 3);

    return true;
}

bool
ForceField::getXMLForceConstants(const XMLParser& xml, RectMatrixPtr& fc)
{
    XMLParser fcnode = xml->fetchNode("fc");
    if (fcnode.get() == NULL)
        return false;

    string fctype = "xyz";
    if (fcnode->hasAttribute("type"))
        fcnode->getAttribute(fctype, "type");

    if (fctype != "xyz")
    {
        except("Cannot read in internals force constants. Please give only xyz force constants.");
    }

    string fcunits = "hartree/bohr";
    if (fcnode->hasAttribute("units"))
        fcnode->getAttribute(fcunits, "units");

    if (fcunits != "hartree/bohr")
    {
        except("Please input force constants in XML file in hartree/bohr units");
    }

    string fctext = fcnode->getText();
    int ncoords = mol_->natoms() * 3;
    fc = getMatrix(fctext, ncoords, ncoords);

    return true;
}

void
ForceField::readXMLDisplacement(const XMLParser& node, bool& extrazero)
{
    if (node.get() == NULL)
        except("How the fuck did this happen?");

    int nvalue = KeywordSet::getKeyword("nvalue")->getValueInteger();
    string energyunits = KeywordSet::getKeyword("energy units")->getValueString();
    
    int dispnumber = -1;
    if (node->hasAttribute("number"))
    {
        string ntext; node->getAttribute(ntext, "number");
        stringstream sstr(ntext);
        sstr >> dispnumber;
    }
    XMLParser molnode = node->fetchNode("molecule");

    RectMatrixPtr geom;
    RectMatrixPtr gradients;
    RectMatrixPtr fc;

    double energy;

    bool foundxyz = getXMLXYZ(molnode, geom);

    if (!foundxyz)
    {
        cerr << "No xyz coordinates for displacement";
        if (dispnumber != -1)
            cerr << " " << dispnumber << ".";
        cerr << " Must have xyz coordinates to determine displacement." << endl;
    }

    bool founde = getXMLEnergy(node, energy);
    bool foundgrads = getXMLGradients(node, gradients);
    bool foundfc = getXMLForceConstants(node, fc);

    int geomdebug = KeywordSet::getKeyword("geometry debug")->getValueInteger();

    //from the geometry determine the displacement
    vector<double> increments = disp_iter_->computeIncrements(geom);

    bool iszero = true;
    for (int i=0; i < increments.size(); i++)
    {
        if ( fabs(increments[i]) > ZERODISP )
        {
            iszero = false;
            break;
        }
    }

    DisplacementPtr disp = disp_iter_->findDisplacement(increments);

    if (disp.get() == NULL && nvalue > 0) //see if it is the zero displacement... this gets tacked on after
    {
        if (iszero)
        {
            disp_iter_->addZeroDisplacement();
            disp = disp_iter_->findDisplacement(increments);
            disp->generateDisplacement();
            extrazero = true;
        }
    }

    if (disp.get() == NULL)
    {
        cerr << "Invalid displacement read in" << endl;
        print_vector(increments);
        geom.print("Invalid geometry");
        return;
    }

    cout << "Assigning values for displacement "; 
    if (dispnumber != -1)
        cout << dispnumber << " ";
    print_vector(increments, true, "%12.8f");

    if (disp->isAssigned())
    {
        stringstream sstr;
        sstr << "Multiple values read in for displacement" << endl;
        print_vector(increments, sstr);
        geom.print("Invalid geometry", sstr);
        except(sstr.str());
    }


    if (nvalue >= 1 && geomdebug >= 1)
    {
        ConstRectMatrixPtr exact_xyz = disp->getDisplacementMolecule()->getXYZ(); 
        vector<double> exact_increments = disp_iter_->computeIncrements(exact_xyz);
        print_vector(exact_increments, false, "%12.8f"); 
    }


    if (nvalue == 1 && !foundgrads)
    {
        stringstream sstr;
        sstr << "nvalue = 1 but no gradients found for displacement";
        if (dispnumber == -1) //no displacment, print xyz
        {
            geom.print("XYZ", sstr);
        }
        else
        {
            sstr << " " << dispnumber << endl;
        }
        except(sstr.str());
    }

    if (KeywordSet::getKeyword("energy units")->getValueString() == "aj")
        energy = energy * HARTREE_TO_AJ;
    
    if (founde)
    {
        disp->assignEnergy(energy);
        if (nvalue >= 1 && iszero) //assign these gradients rigorously as the gradients for this force field
            assignEnergy(energy);
    }
    else if (nvalue == 0 || energy == 0.0)
    {
        stringstream sstr;
        sstr << "No energy found for displacement";
        if (dispnumber == -1) //no displacment, print xyz
        {
            geom.print("XYZ", sstr);
        }
        else
        {
            sstr << " " << dispnumber << endl;
        }
        except(sstr.str());
    }

    VectorPtr intgrads;
    RectMatrixPtr xyzgrads;
    RectMatrixPtr intfc;
    RectMatrixPtr xyzfc;
    if (foundgrads && nvalue >= 1)
    {
        MoleculePtr dispmol = mol_->copy(); dispmol->setXYZ(geom);

        intgrads = IntderWriter::transformGradientsToInternals(
                                                            gradients, 
                                                            dispmol, 
                                                            simples_, 
                                                            coords_
                                                           );
        
        xyzgrads = IntderWriter::transformGradientsToCartesian(
                                                            intgrads, 
                                                            disp->getDisplacementMolecule(), 
                                                            simples_, 
                                                            coords_);
        disp->assignGradients(xyzgrads);

        if (nvalue >= 1 && iszero) //assign these gradients rigorously as the gradients for this force field
            assignGradients(xyzgrads);
    }
    if (foundfc && nvalue >= 2)
    {
        MoleculePtr dispmol = mol_->copy(); 
        dispmol->setXYZ(geom);

        intfc = IntderWriter::transformForceConstantsToInternals(
                                                            fc, 
                                                            gradients,
                                                            dispmol, 
                                                            simples_, 
                                                            coords_);

        xyzfc = IntderWriter::transformForceConstantsToCartesian(
                                                            intfc, 
                                                            intgrads,
                                                            disp->getDisplacementMolecule(), 
                                                            simples_, 
                                                            coords_);
        disp->assignForceConstants(xyzfc, xyzgrads);

        if (nvalue >= 2 && iszero) //assign these rigorously as the fcs for this force field
            assignForceConstants(xyzfc, xyzgrads);

    }
}

void
ForceField::readXMLData(string filename)
{
    cout << stream_printf("Reading in data for %d unique displacements, %d total displacements",
                     disp_iter_->nunique(), disp_iter_->ndisps()) << endl << endl;

    
    XMLParser xml(new pyxml::PyXMLDomParser(filename));
    vector<XMLParser> nodes; 
    xml->getElementsByTagName(nodes, "displacement");

    int nvalue = KeywordSet::getKeyword("nvalue")->getValueInteger();
    if (nvalue >= 1) //the displacements have to know their geometries
    {
        int idisp = 1;
        int ndisp = disp_iter_->ndisps();
        for (DisplacementIterator::iterator it(disp_iter_->begin()); it != disp_iter_->end(); ++it, ++idisp)
        {
            cout << stream_printf("Generating displacement %d out of %d", idisp, ndisp) << endl;
            (*it)->generateDisplacement();
        }
    }
    
    int ntot = 0;
    bool extrazero = false;
    vector<XMLParser>::iterator it;
    for (it=nodes.begin(); it != nodes.end(); ++it, ++ntot)
    {
        XMLParser disp = *it;
        if (disp.get() == NULL)
            except(stream_printf("XML node %d is null", ntot + 1));
        readXMLDisplacement(disp, extrazero);
    }

    int nunique = disp_iter_->nunique();
    if (extrazero)
        nunique++;

    if (ntot != nunique)
    {
        cout << stream_printf("WARNING! Read in %d displacements, but there are %d unique displacements",
                         ntot, disp_iter_->nunique()) << endl << endl;
    }

    //loop through and figure out which displacements didn't get assigned a value
    stringstream sstr;
    vector<double> disps = GigideKeyword::getDisplacementSizes(coords_.size());

    vector<DisplacementPtr> missing;
    for (DisplacementIterator::iterator it (disp_iter_->begin()); it != disp_iter_->end(); ++it)
    {
        DisplacementPtr disp = *it;
        if (disp->energyAssigned())
            continue; //nothing to do

        missing.push_back(disp);
        cerr << disp->label() << " not assigned energy" << endl;
    }

    if (missing.size())
    {
        ofstream dispcart;
        dispcart.open("dispcart_missing");
        for (int i=0; i < missing.size(); ++i)
        {
            DisplacementPtr disp(missing[i]);
            disp->generateDisplacement();
            writeMoleculeToDispcart(dispcart, disp->getDisplacementMolecule(), i + 1);
        }
        dispcart.close();

        //generate a dispcart file of the missing displacements
        except(stream_printf("%d displacements not assigned energy", missing.size()));
    }
}

ConstVectorPtr
ForceField::getDerivativeValues(const ConstDerivativePtr& deriv) const
{
    int ndisps = disp_iter_->ndisps();

    //the number of displacements should be computed from fit as extra points can get tacked on to
    //the displacement iterator and they should not be included

    VectorPtr derivs = new Vector(ndisps);
    int num = 0;
    for (DisplacementIterator::iterator it(disp_iter_->begin()); it != disp_iter_->end(); ++it, ++num)
    {
        DisplacementPtr disp = *it;
        if (disp.get() == NULL)
        {
            except("Got null displacement in derivative iteration.  This should not be possible");
        }
        double val = disp->getDerivativeValue(deriv);
        derivs.set_element(num, val);
    }
    return derivs;
}

void
ForceField::compute()
{
    double econv = 1.0;
    if (KeywordSet::getKeyword("energy units")->getValueString() == "aj")
        econv = 1.0/HARTREE_TO_AJ;

    //loop throug the energy points and try to determine if any are
    //wiggidy wack
    double emin = 0;
    DisplacementPtr mindisp;
    for (DisplacementIterator::iterator it(disp_iter_->begin()); it != disp_iter_->end(); ++it)
    {
        DisplacementPtr disp = *it;
        double edisp = disp->getEnergy() * econv;
        if (edisp < emin)
        {
            mindisp = disp;
            emin = edisp;
        }
    }

    if (mindisp.get() != NULL)
    {
            cout << "Minimum energy displacement is "; mindisp->print(); cout << endl;
            //cout << stream_printf("Energy =  %16.12f", emin) << endl << endl;
    }

    DisplacementPtr center = disp_iter_->findZeroDisplacement();
    if (center.get() != NULL)
    {
        double ecenter = center->getEnergy() * econv;
        if ( fabs(ecenter - emin) > 1e-12 )
            cout << stream_printf("WARNING! Center displacement has energy %14.10f but lowest energy is %14.10f",
                             ecenter, emin) << endl << endl;
    }
    else
    {
        cout << "No central displacement.  Cannot verify validity of minimum energy point" << endl << endl;
    }

    //arbitrary cutoff for too large a finite displacement difference
    double maxdiff = fabs(emin * 0.0002);
    double mindiff = fabs(emin * 1e-8);

    bool pointstoolarge = false;
    bool pointstoosmall = false;
    for (DisplacementIterator::iterator it(disp_iter_->begin()); it != disp_iter_->end(); ++it)
    {
        DisplacementPtr disp = *it;
        if (disp->getRefcount() == 0)
            continue; //this displacement is not really needed

        double edisp = disp->getEnergy() * econv;
        double diff = edisp - emin;
        if (!disp->energyAssigned()) //never got assigned an energy
        {
            cout << "WARNING! Displacement never got assigned energy "; disp->print();
            ConstMoleculePtr mol  = disp->getDisplacementMolecule();
            if (mol.get() == NULL)
                cout << mol->getXYZString() << endl;
            else
                cout << endl;
            
        }
        else if (diff > maxdiff)
        {
            cout << "WARNING! Displacement might be wrong "; disp->print();
            cout << stream_printf("Energy %16.12f differs by %6.4f from minimum.  Check convergence or consider smaller displacements.",
                             edisp, diff) << endl << endl;
            pointstoolarge = true;
        }
        else if (diff < mindiff && fabs(diff) > 1e-14)
        {
            cout << "WARNING! Displacement might give bad results "; disp->print();
            cout << stream_printf("Energy %16.12f only differs by %6.4e from minimum.  Check convergence or consider larger displacements.",
                             edisp, diff) << endl << endl;
            pointstoosmall = true;
        }
    }
    
    if (!pointstoolarge && !pointstoosmall)
    {
        cout << "Force field appears valid.  No points appear to have erroneous energies." << endl;
    }
    if (!pointstoolarge)
        cout << stream_printf("All energy differences are smaller than %12.8f", maxdiff) <<endl;
    if (!pointstoosmall)
        cout << stream_printf("All energy differences are larger than %12.8f", mindiff) <<endl;

    //loop through all the provided values and differentiate them by finite difference
    //the nvalue + 1 says that we should only compute by finite difference those values
    //that we can't get exactly from the midpoint
    int nvalue = KeywordSet::getKeyword("nvalue")->getValueInteger();
    int fitdebug = KeywordSet::getKeyword("fit debug")->getValueInteger();
    for (DerivativeIterator::iterator it(value_iter_->begin()); it != value_iter_->end(); ++it)
    {
        DerivativePtr value = *it;
        if (value->level() != nvalue) //not something to differentiate
            continue;

        if (fitdebug)
        {
            cout << "Computing values for derivative" << endl;
            value->print();
        }
        computeDerivatives(value);
    }

    if (disp_iter_->hasZeroDisplacement())
    {
        //assign rigorously accurate values from the center displacement, if applicable
        for (DerivativeIterator::iterator it(value_iter_->begin()); it != value_iter_->end(); ++it)
        {
            DisplacementPtr disp = disp_iter_->findZeroDisplacement();
            DerivativePtr deriv = *it;

            if (deriv->level() != nvalue) //not something to differentiate
                continue;

            double value = disp->getDerivativeValue(deriv);
            if (fitdebug)
            {
                cout << "Assigning value " << value << " for derivative" << endl;
                deriv->print();
            }
            deriv->setValue(value);
        }
    }

    //now compute the taylor series approximations
    vector<TaylorSeriesEnergyPtr> energies;
    TaylorSeriesEnergy::buildEnergyApproximations(energies, disp_iter_, deriv_iter_);

    vector<TaylorSeriesEnergyPtr>::const_iterator it(energies.begin());
    double error = 0;
    for ( ; it != energies.end(); ++it)
    {
        TaylorSeriesEnergyPtr energy(*it);
        double pterror = energy->error();
        error += pterror * pterror;
        //energy->print();
    }

    error = sqrt(error / energies.size());
    error = Units::convert(error, Units::Hartree) * 1e6; //go to microhartree
    cout << stream_printf("RMS Error in Fit for all points: %12.8f uH", error) << endl;

    //now that all derivatives have been computed, build the important arrays
    buildArrays();

    //if we have computed second derivatives, do a frequency analysis
    if (deriv_iter_->level() >= 2)
        computeFrequencies();

    bool statpt = KeywordSet::getKeyword("statpt")->getValueBoolean();
    //and finally write an intder file
    vector<InternalCoordinatePtr> symms;
    
    IntderWriter iderarch(mol_, statpt, IntderWriter::TransformInternals, deriv_iter_->level(), deriv_iter_, simples_, coords_);
    iderarch.commit(DEFAULT_INTDER_INPUT_FILE);

    if (nderiv() >= 3) //if we are doing 3rd or higher derivatives, write an anharm file
    {
        AnharmWriter anharch(mol_, this);
        anharch.commit(DEFAULT_ANHARM_INPUT_FILE);
    }


}

void
ForceField::assignGradients(ConstRectMatrixPtr xyzgrads)
{
    VectorPtr intgrads = IntderWriter::transformGradientsToInternals(xyzgrads, mol_, simples_, coords_);
    int ngrads = intgrads.n();
    for (int i=0; i < ngrads; ++i)
    {
        vector<int> indices(ngrads, 0);
        indices[i] = 1;
        DerivativePtr deriv = deriv_iter_->getDerivative(indices);
        deriv->setValue( intgrads.get_element(i) );
    }
}

void
ForceField::assignEnergy(double energy)
{
    energy_ = energy;
}

int
ForceField::nderiv() const
{
    return deriv_iter_->level();
}

void
ForceField::assignForceConstants(ConstRectMatrixPtr xyzfc, ConstRectMatrixPtr xyzgrads)
{
    RectMatrixPtr fc = IntderWriter::transformForceConstantsToInternals(xyzfc, xyzgrads, mol_, simples_, coords_);
    int n = fc.nrow();
    for (int i=0; i < n; i++)
    {
        for (int j=0; j < n; j++)
        {
            vector<int> indices(n,0);
            indices[i] += 1;
            indices[j] += 1;
            DerivativePtr deriv = deriv_iter_->getDerivative(indices);
            deriv->setValue( fc.get_element(i,j) );
        }
    }
}

void
ForceField::print(ostream& os) const
{
    os << stream_printf("Force field for %dth derivatives from %dth derivative displacements", nderiv_, nvalue_) << endl;
    mol_->print(os);
    disp_iter_->print(os); os <<endl;
    if (grads_.nonnull())
        grads_.print("Gradients", os);
    if (fc_.nonnull())
        fc_.print("Force constants", os);
}

ConstVectorPtr
ForceField::gradients() const
{
    return grads_;
}

double
ForceField::energy() const
{
    return energy_;
}

ConstSymmMatrixPtr
ForceField::fc() const
{
    return fc_;
}

vector<string> Derivative::letters_;
bool Derivative::initdone_ = false;
map<string, Derivative::DerivativeType> Derivative::typemap_;
