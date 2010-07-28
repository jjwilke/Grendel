#ifndef gigide_deftypes_h
#define gigide_deftypes_h

#include <src/smartptr/ref.h>

#define gigstr const std::string&

#define InternalCoordinatePtr boost::intrusive_ptr<InternalCoordinate>
#define SimpleInternalCoordinatePtr  boost::intrusive_ptr<SimpleInternalCoordinate>
#define SymmetryInternalCoordinatePtr  boost::intrusive_ptr<SymmetryInternalCoordinate>
#define LinXPtr  boost::intrusive_ptr<LinX>
#define LinYPtr  boost::intrusive_ptr<LinY>
#define Lin1Ptr  boost::intrusive_ptr<Lin1>
#define BondAnglePtr  boost::intrusive_ptr<BondAngle>
#define BondLengthPtr  boost::intrusive_ptr<BondLength>   
#define TorsionPtr  boost::intrusive_ptr<Torsion>
#define OutOfPlaneBendPtr  boost::intrusive_ptr<OutOfPlaneBend>  
#define PeriodicCoordinatePtr  boost::intrusive_ptr<PeriodicCoordinate>   

#define DisplacementMappingPtr boost::intrusive_ptr<DisplacementMapping>   
#define DisplacementPtr boost::intrusive_ptr<Displacement>
#define DisplacementIteratorPtr boost::intrusive_ptr<DisplacementIterator>   

#define DerivativePtr boost::intrusive_ptr<Derivative>
#define DerivativeIteratorPtr boost::intrusive_ptr<DerivativeIterator>  

#define FitPtr boost::intrusive_ptr<Fit>
#define FitFactoryPtr boost::intrusive_ptr<FitFactory>
#define FitPointPtr boost::intrusive_ptr<FitPoint>

#define XYZPointPtr boost::intrusive_ptr<XYZPoint>
#define AxisPtr boost::intrusive_ptr<Axis>
#define MidpointPtr boost::intrusive_ptr<Midpoint>

#define EmpiricalHessianTermPtr boost::intrusive_ptr<EmpiricalHessianTerm>
#define ConstantTermPtr boost::intrusive_ptr<ConstantTerm>
#define BondBondPtr boost::intrusive_ptr<BondBond_EmpiricalHessian>
#define BendBendPtr boost::intrusive_ptr<BendBend_EmpiricalHessian> 
#define BondBendPtr boost::intrusive_ptr<BondBend_EmpiricalHessian>  
#define TorsTorsPtr boost::intrusive_ptr<TorsTors_EmpiricalHessian>   

#define InputFilePtr boost::intrusive_ptr<InputFile>
#define GigideInputFilePtr boost::intrusive_ptr<GigideInputFile>

#define KeywordValuePtr boost::intrusive_ptr<KeywordValue>
#define KeywordSetPtr boost::intrusive_ptr<KeywordSet>
#define KeywordIteratorPtr boost::intrusive_ptr<KeywordIterator>

#define MoleculePtr boost::intrusive_ptr<Molecule>  
#define ConstMoleculePtr boost::intrusive_ptr<const Molecule>  
#define AtomPtr boost::intrusive_ptr<Atom>  

#define PermutationGeneratorIntPtr boost::intrusive_ptr<PermutationGenerator<int> >  
#define PermutationGeneratorDoublePtr boost::intrusive_ptr<PermutationGenerator<double> >  



#define ForceFieldPtr boost::intrusive_ptr<ForceField>  

#define SymmetryOperationPtr boost::intrusive_ptr<SymmetryOperation>
#define ImproperRotationPtr boost::intrusive_ptr<ImproperRotation>
#define RotationPtr boost::intrusive_ptr<Rotation>
#define ReflectionPtr boost::intrusive_ptr<Reflection>
#define InversionPtr boost::intrusive_ptr<Inversion>
#define IdentityElementPtr boost::intrusive_ptr<IdentityElement>
#define PointGroupClassPtr boost::intrusive_ptr<PointGroupClass>
#define PointGroupPtr boost::intrusive_ptr<PointGroup>
#define CoordinateSubspacePtr boost::intrusive_ptr<CoordinateSubspace>  

#define TaylorTermPtr boost::intrusive_ptr<TaylorTerm>  

#define OptimizationPtr boost::intrusive_ptr<Optimization>  
#define OptimizationStepPtr boost::intrusive_ptr<OptimizationStep>  

//#define MatrixPtr boost::intrusive_ptr<Matrix>    
//#define VectorPtr boost::intrusive_ptr<Vector>  

#define ConstInternalCoordinatePtr boost::intrusive_ptr<const InternalCoordinate>
#define ConstSimpleInternalCoordinatePtr  boost::intrusive_ptr<const SimpleInternalCoordinate>
#define ConstSymmetryInternalCoordinatePtr  boost::intrusive_ptr<const SymmetryInternalCoordinate>
#define ConstLinXPtr  boost::intrusive_ptr<const LinX>
#define ConstLinYPtr  boost::intrusive_ptr<const LinY>
#define ConstLin1Ptr  boost::intrusive_ptr<const Lin1>
#define ConstBondAnglePtr  boost::intrusive_ptr<const BondAngle>
#define ConstBondLengthPtr  boost::intrusive_ptr<const BondLength>   
#define ConstTorsionPtr  boost::intrusive_ptr<const Torsion>
#define ConstOutOfPlaneBendPtr  boost::intrusive_ptr<const OutOfPlaneBend>  
#define ConstPeriodicCoordinatePtr  boost::intrusive_ptr<const PeriodicCoordinate>   
#define ConstDisplacementMappingPtr boost::intrusive_ptr<const DisplacementMapping>   
#define ConstDisplacementPtr boost::intrusive_ptr<const Displacement>
#define ConstDisplacementIteratorPtr boost::intrusive_ptr<const DisplacementIterator>   
#define ConstDerivativePtr boost::intrusive_ptr<const Derivative>
#define ConstDerivativeIteratorPtr boost::intrusive_ptr<const DerivativeIterator>  
#define ConstFitPtr boost::intrusive_ptr<const Fit>
#define ConstFitFactoryPtr boost::intrusive_ptr<const FitFactory>
#define ConstFitPointPtr boost::intrusive_ptr<const FitPoint>
#define ConstXYZPointPtr boost::intrusive_ptr<const XYZPoint>
#define ConstAxisPtr boost::intrusive_ptr<const Axis>
#define ConstMidpointPtr boost::intrusive_ptr<const Midpoint>
#define ConstEmpiricalHessianTermPtr boost::intrusive_ptr<const EmpiricalHessianTerm>
#define ConstConstantTermPtr boost::intrusive_ptr<const ConstantTerm>
#define ConstBondBondPtr boost::intrusive_ptr<const BondBond_EmpiricalHessian>
#define ConstBendBendPtr boost::intrusive_ptr<const BendBend_EmpiricalHessian> 
#define ConstBondBendPtr boost::intrusive_ptr<const BondBend_EmpiricalHessian>  
#define ConstTorsTorsPtr boost::intrusive_ptr<const TorsTors_EmpiricalHessian>   
#define ConstInputFilePtr boost::intrusive_ptr<const InputFile>
#define ConstGigideInputFilePtr boost::intrusive_ptr<const GigideInputFile>
#define ConstKeywordValuePtr boost::intrusive_ptr<const KeywordValue>
#define ConstKeywordSetPtr boost::intrusive_ptr<const KeywordSet>
#define ConstKeywordIteratorPtr boost::intrusive_ptr<const KeywordIterator>
#define ConstMoleculePtr boost::intrusive_ptr<const Molecule>  
#define ConstConstMoleculePtr boost::intrusive_ptr<const const Molecule>  
#define ConstAtomPtr boost::intrusive_ptr<const Atom>  
#define ConstPermutationGeneratorIntPtr boost::intrusive_ptr<const PermutationGenerator<int> >  
#define ConstPermutationGeneratorDoublePtr boost::intrusive_ptr<const PermutationGenerator<double> >  
#define ConstForceFieldPtr boost::intrusive_ptr<const ForceField>  
#define ConstSymmetryOperationPtr boost::intrusive_ptr<const SymmetryOperation>
#define ConstImproperRotationPtr boost::intrusive_ptr<const ImproperRotation>
#define ConstRotationPtr boost::intrusive_ptr<const Rotation>
#define ConstReflectionPtr boost::intrusive_ptr<const Reflection>
#define ConstInversionPtr boost::intrusive_ptr<const Inversion>
#define ConstIdentityElementPtr boost::intrusive_ptr<const IdentityElement>
#define ConstPointGroupClassPtr boost::intrusive_ptr<const PointGroupClass>
#define ConstPointGroupPtr boost::intrusive_ptr<const PointGroup>
#define ConstCoordinateSubspacePtr boost::intrusive_ptr<const CoordinateSubspace>  
#define ConstTaylorTermPtr boost::intrusive_ptr<const TaylorTerm>  
//#define ConstMatrixPtr boost::intrusive_ptr<const Matrix>    
//#define ConstRectMatrixPtr boost::intrusive_ptr<const RectMatrix>  
//#define ConstSymmMatrixPtr boost::intrusive_ptr<const SymmMatrix>  
//#define ConstVectorPtr boost::intrusive_ptr<const Vector>  

#endif
