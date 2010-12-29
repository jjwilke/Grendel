#ifndef gigide_hessian_hpp
#define gigide_hessian_hpp

namespace gigide {

class EmpiricalHessianTerm;
class ConstantTerm;
class BondBond_EmpiricalHessian;
class BendBend_EmpiricalHessian;
class BondBend_EmpiricalHessian;
class TorsTors_EmpiricalHessian;

typedef boost::intrusive_ptr<EmpiricalHessianTerm> EmpiricalHessianTermPtr;
typedef boost::intrusive_ptr<ConstantTerm> ConstantTermPtr;
typedef boost::intrusive_ptr<BondBond_EmpiricalHessian> BondBondPtr;
typedef boost::intrusive_ptr<BendBend_EmpiricalHessian>  BendBendPtr;
typedef boost::intrusive_ptr<BondBend_EmpiricalHessian>   BondBendPtr;
typedef boost::intrusive_ptr<TorsTors_EmpiricalHessian>    TorsTorsPtr;

typedef boost::intrusive_ptr<const EmpiricalHessianTerm> ConstEmpiricalHessianTermPtr;
typedef boost::intrusive_ptr<const ConstantTerm> ConstConstantTermPtr;
typedef boost::intrusive_ptr<const BondBond_EmpiricalHessian> ConstBondBondPtr;
typedef boost::intrusive_ptr<const BendBend_EmpiricalHessian>  ConstBendBendPtr;
typedef boost::intrusive_ptr<const BondBend_EmpiricalHessian>   ConstBondBendPtr;
typedef boost::intrusive_ptr<const TorsTors_EmpiricalHessian>    ConstTorsTorsPtr;

}

#endif

