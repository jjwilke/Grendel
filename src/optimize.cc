#include <src/optimize.h>
#include <src/utilities.h>
#include <src/derivative.h>
#include <src/exception.h>

using namespace gigide;
using namespace smartptr;
using namespace std;

SerialDeclare(Optimization);
SerialDeclare(OptimizationStep);

Optimization::Optimization(
   const ConstMoleculePtr& mol,
   const Set<ConstInternalCoordinatePtr>& coords,
   const Set<ConstSimpleInternalCoordinatePtr>& simples
) : mol_(mol), coords_(coords), simples_(simples)
{
    SetRuntime(Optimization);

#if 0
    //check for hessian
    ifstream ifile("hessian.xml");
    if (ifile) //we have a hessian file
    {
        ArchivePtr arch(new smartptr::Archive("hessian.xml"));
        serial_load(hessian);
        precond_ = hessian_.i();
    }
    else
    {
        //use identity
        hessian_ = SymmMatrixPtr(coords_.size());
        for (int i=0; i < hessian_.n(); ++i) 
            hessian_.set_element(i,i,1.0);
        precond_ = hessian_.copy();
    }
#endif
}

Optimization::Optimization(const ArchivePtr& arch)
    : Serializable(arch)
{
    SetRuntime(Optimization);
    serial_load(hessian);
    serial_load(precond);
    serial_load(simples);
    serial_load(coords);
}

void
Optimization::serialize(const ArchivePtr& arch) const
{
    Serializable::serialize(arch);
    serial_save(hessian);
    serial_save(precond);
    serial_save(simples);
    serial_save(coords);
}

void
Optimization::writeDisplacements()
{
    OptimizationStepPtr next_step(new OptimizationStep(mol_, coords_, simples_));
    steps_.push_back(next_step);
    next_step->writeDisplacements();
}

void
Optimization::takeStep()
{
    if (steps_.size() == 0)
    {
        except("Cannot take step.  No gradients have been generated.");
    }

    OptimizationStepPtr step = *(steps_.end() - 1);
    VectorPtr disp = precond_ * step->gradients();
    OptimizationStepPtr next = step->takeStep(disp);
    steps_.push_back(step);
}

OptimizationStep::OptimizationStep(
    const ConstMoleculePtr& mol,
    const Set<ConstInternalCoordinatePtr>& coords,
    const Set<ConstSimpleInternalCoordinatePtr>& simples
) : computed_(false), qff_(new ForceField(mol, coords, simples)),
    mol_(mol),
    coords_(coords),
    simples_(simples)
{
    SetRuntime(OptimizationStep);
}

OptimizationStep::OptimizationStep(const ArchivePtr& arch)
{
    SetRuntime(OptimizationStep);
    serial_load(qff);
}

void
OptimizationStep::serialize(const ArchivePtr& arch) const
{
    serial_save(qff);
}

void
OptimizationStep::writeDisplacements()
{
    qff_->generateDisplacements();
    qff_->writeDisplacementsToFile("dispcart");
}

ConstVectorPtr
OptimizationStep::gradients()
{
    if (!computed_)
    {
        qff_->compute(); 
        computed_ = true;
    }

    return qff_->gradients();
}

double
OptimizationStep::energy()
{
    if (!computed_)
    {
        qff_->compute(); 
        computed_ = true;
    }

    return qff_->energy();
}

OptimizationStepPtr
OptimizationStep::takeStep(ConstVectorPtr dispvec) const
{
    vector<double> disps;
    for (int i=0; i < dispvec.n(); ++i) 
        disps.push_back(dispvec[i]);

    vector<InternalCoordinatePtr> newcoords;
    vector<SimpleInternalCoordinatePtr> newsimples;
    MoleculePtr displmol = Displacement::displaceGeometry(disps, mol_, coords_, simples_, newcoords, newsimples);
}

