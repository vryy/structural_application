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
class PrescribedObject
{

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( PrescribedObject );

    ///@}
    ///@name Operation
    ///@{

    /// Modify the RHS to account for prescribed values at nodes
    virtual void ApplyPrescribedDofs(const Matrix& LHS_Contribution, Vector& RHS_Constribution, const ProcessInfo& CurrentProcessInfo) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    /// Compute the forces induces by prescribed dofs
    virtual void ComputePrescribedForces(const Matrix& LHS_Contribution, Vector& Force, const ProcessInfo& CurrentProcessInfo) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

}; // Class PrescribedObject

}  // namespace Kratos.

#endif // KRATOS_PRESCRIBED_OBJECT_H_INCLUDED defined
