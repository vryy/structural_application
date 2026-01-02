//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 7/3/2022 $
//
//

#if !defined(KRATOS_PRESCRIBED_OBJECT_H_INCLUDED )
#define  KRATOS_PRESCRIBED_OBJECT_H_INCLUDED

#include "includes/define.h"

namespace Kratos
{

/**
 * An interface that supports to prescribe the primal value, i.e. nodal displacement via scheme
 * TODO put this interface to the element abstract interface
 */
template<typename TNodeType>
class BasePrescribedObject
{

public:
    ///@name Type Definitions
    ///@{

    typedef typename TNodeType::DofType::DataType DataType;
    typedef typename MatrixVectorTypeSelector<DataType>::MatrixType MatrixType;
    typedef typename MatrixVectorTypeSelector<DataType>::VectorType VectorType;

    KRATOS_CLASS_POINTER_DEFINITION( BasePrescribedObject );

    ///@}
    ///@name Operation
    ///@{

    /// Modify the RHS to account for prescribed values at nodes
    virtual void ApplyPrescribedDofs(const MatrixType& LHS_Contribution,
        VectorType& RHS_Constribution,
        const ProcessInfo& CurrentProcessInfo) const = 0;

    /// Compute the forces induces by prescribed dofs
    virtual void ComputePrescribedForces(const MatrixType& LHS_Contribution,
        VectorType& Force,
        const ProcessInfo& CurrentProcessInfo) const = 0;

    ///@}
}; // Class BasePrescribedObject

typedef BasePrescribedObject<RealNode> PrescribedObject;

}  // namespace Kratos.

#endif // KRATOS_PRESCRIBED_OBJECT_H_INCLUDED defined
