#ifndef _gigide_optimize_h_
#define _gigide_optimize_h_

#include "gigide.hpp"

namespace gigide {

class OptimizationStep : public smartptr::Serializable {

    private:
        ForceFieldPtr qff_;
        ConstMoleculePtr mol_;
        smartptr::Set<ConstInternalCoordinatePtr> coords_;
        smartptr::Set<ConstSimpleInternalCoordinatePtr> simples_;

        bool computed_;

    public:
        OptimizationStep(
            const ConstMoleculePtr& mol,
            const smartptr::Set<ConstInternalCoordinatePtr>& coords,
            const smartptr::Set<ConstSimpleInternalCoordinatePtr>& simples
        );

        OptimizationStep(const XMLArchivePtr& arch);

        void serialize(const XMLArchivePtr& arch) const;

        void writeDisplacements();

        double energy();

        ConstVectorPtr gradients();

        ConstForceFieldPtr getForceField() const;

        OptimizationStepPtr takeStep(ConstVectorPtr disps) const;
};

class Optimization : public smartptr::Serializable {

    private:
        std::vector<OptimizationStepPtr> steps_;

        /** The force constant matrix */
        SymmMatrixPtr hessian_;

        /** The preconditioner, usually the inverse hessian */
        SymmMatrixPtr precond_;

        double tol_;

        ConstMoleculePtr mol_;
        smartptr::Set<ConstInternalCoordinatePtr> coords_;
        smartptr::Set<ConstSimpleInternalCoordinatePtr> simples_;

    public:
        Optimization(
            const ConstMoleculePtr& mol,
            const smartptr::Set<ConstInternalCoordinatePtr>& coords,
            const smartptr::Set<ConstSimpleInternalCoordinatePtr>& simples
        );

        Optimization(const XMLArchivePtr& arch);

        void serialize(const XMLArchivePtr& arch) const;

        void takeStep();

        void writeDisplacements();

};

}

#endif

