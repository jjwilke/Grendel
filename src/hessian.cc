
#include <src/hessian.h>
#include <iostream>

using namespace gigide;
using namespace std;

#undef TEST
#define TEST(x1,x2,t1,t2,classname) if (x1->type() == t1 && x2->type() == t2) return new classname(x1,x2)
#define TEST_swap(x1,x2,t1,t2,classname) if (x1->type() == t1 && x2->type() == t2) return new classname(x2,x1)

EmpiricalHessianTerm::EmpiricalHessianTerm(
    const ConstSimpleInternalCoordinatePtr& coord1,
    const ConstSimpleInternalCoordinatePtr& coord2
) : c1_(coord1), c2_(coord2)
{
}

EmpiricalHessianTerm::EmpiricalHessianTerm(
)
{
}

EmpiricalHessianTermPtr
EmpiricalHessianTerm::getTerm(
    const ConstSimpleInternalCoordinatePtr& c1,
    const ConstSimpleInternalCoordinatePtr& c2
)
{
    TEST(c1,c2,"STRE","STRE",BondBond_EmpiricalHessian);
    TEST(c1,c2,"STRE","BEND",BondBend_EmpiricalHessian);
    TEST_swap(c1,c2,"BEND","STRE",BondBend_EmpiricalHessian);
    TEST(c1,c2,"BEND","BEND",BendBend_EmpiricalHessian);
    TEST(c1,c2,"TORS","TORS",TorsTors_EmpiricalHessian);

    if (c1 == c2) //just return default diagonal term
        return new ConstantTerm(0.5);

    return new ConstantTerm(0.0);
}

bool
EmpiricalHessianTerm::matches(
    int id1, 
    int id2,
    ConstAtomPtr& at
) const
{
    ConstAtomPtr at1 = c1_->getAtom(id1);
    ConstAtomPtr at2 = c2_->getAtom(id2);

    if (at1 == at2)
    {
        at = at1;
        return true;
    }
    else
        return false;
}

ConstantTerm::ConstantTerm(
    double val
)
{
    val_ = val;
}

double
ConstantTerm::getValue() const
{
    return val_;
}

#undef CONSTRUCT
#define CONSTRUCT(x1) x1::x1(const ConstSimpleInternalCoordinatePtr& c1, const ConstSimpleInternalCoordinatePtr& c2) : EmpiricalHessianTerm(c1,c2) {}

CONSTRUCT(BondBond_EmpiricalHessian)
CONSTRUCT(BondBend_EmpiricalHessian)
CONSTRUCT(BendBend_EmpiricalHessian)
CONSTRUCT(TorsTors_EmpiricalHessian)

double
BondBond_EmpiricalHessian::getValue() const
{
    ConstAtomPtr at1, at2;
    if (matches(0,0,at1) && matches(1,1,at2))
        return xxtype(at1,at2);
    if (matches(0,1,at1) && matches(1,0,at2))
        return xxtype(at1,at2);
    if (matches(0,0,at1))
        return xytype(at1,c1_->getAtom(1),c2_->getAtom(1));
    if (matches(0,1,at1))
        return xytype(at1,c1_->getAtom(1),c2_->getAtom(0));
    if (matches(1,0,at1))
        return xytype(at1,c1_->getAtom(0),c2_->getAtom(1));
    if (matches(1,1,at1))
        return xytype(at1,c1_->getAtom(0),c2_->getAtom(0));

    //no other matches
    return 0.0;
}

double
BondBond_EmpiricalHessian::xxtype(
    const ConstAtomPtr& at1,
    const ConstAtomPtr& at2
) const
{
    return 10; //just return 10 for now
}

double
BondBond_EmpiricalHessian::xytype(
    const ConstAtomPtr& at1,
    const ConstAtomPtr& at2,
    const ConstAtomPtr& at3
) const
{
    return 0; //just return 1 for now
}

double
BondBend_EmpiricalHessian::getValue() const
{
    ConstAtomPtr at1, at2, at3;
    if (matches(0,0,at1) && matches(1,1,at2))
    {
        at3 = c2_->getAtom(2);
        return xxytype(at1,at2,at3);
    }
    if (matches(0,2,at1) && matches(1,1,at2))
    {
        at3 = c2_->getAtom(0);
        return xxytype(at1,at2,at3);
    }

    //no other matches
    return 0.0;
}

double
BondBend_EmpiricalHessian::xxytype(
    const ConstAtomPtr& at1,
    const ConstAtomPtr& at2,
    const ConstAtomPtr& at3
) const
{
    return 0.1;
}

double
BendBend_EmpiricalHessian::getValue() const
{
    ConstAtomPtr at1, at2, at3;
    if (matches(0,0,at1) && matches(1,1,at2) && matches(2,2, at3))
        return xxxtype(at1,at2,at3);
    if (matches(2,0,at1) && matches(1,1,at2) && matches(0,2, at3))
        return xxxtype(at1,at2,at3);

    //no other matches
    return 0.0;
}

double
BendBend_EmpiricalHessian::xxxtype(
    const ConstAtomPtr& at1,
    const ConstAtomPtr& at2,
    const ConstAtomPtr& at3
) const
{
    return 0.5;
}

double
TorsTors_EmpiricalHessian::getValue() const
{
    ConstAtomPtr at1, at2, at3, at4;
    if (matches(0,0,at1) && matches(1,1,at2) && matches(2,2, at3) && matches(3,3,at4))
        return xxxxtype(at1,at2,at3,at4);
    if (matches(3,0,at1) && matches(2,1,at2) && matches(1,2, at3) && matches(0,3,at4))
        return xxxxtype(at1,at2,at3,at4);

    //no other matches
    return 0.0;
}

double
TorsTors_EmpiricalHessian::xxxxtype(
    const ConstAtomPtr& at1,
    const ConstAtomPtr& at2,
    const ConstAtomPtr& at3,
    const ConstAtomPtr& at4
) const
{
    return 0.25;
}

