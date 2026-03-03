//   Project Name:        KratosStructuralApplication
//   Modified by:         $Author: hbui $
//   Date:                $Date: 2/3/2026 $
//
//

#if !defined(KRATOS_TIME_INTEGRABLE_OBJECT_H_INCLUDED )
#define  KRATOS_TIME_INTEGRABLE_OBJECT_H_INCLUDED

#include "includes/define.h"

namespace Kratos
{

enum class NonlinearMassDampingType
{
    NO_MASS_DAMPING                 = 0,
    LINEAR_MASS_DAMPING             = 1,
    NONLINEAR_MASS_DAMPING          = 2,
    REDUCED_NONLINEAR_MASS_DAMPING  = 3,
    TOTAL_NONLINEAR_MASS_DAMPING    = 4,
    UNKNOWN                         = -1
};

inline std::ostream& operator <<(std::ostream& rOStream, const NonlinearMassDampingType& rThis)
{
    rOStream << static_cast<int>(rThis);
    return rOStream;
}

/**
 * An interface to select the method to compute the local contribution for time integration scheme
 */
class TimeIntegrableObject
{

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( TimeIntegrableObject );

    ///@}
    ///@name Inquiry
    ///@{

    /// Return the elemental nonlinear mass damping approach, i.e., the way the element calculates
    /// the residual and stiffness with regards to the time integration
    virtual NonlinearMassDampingType GetNonlinearMassDampingApproach() const = 0;

    ///@}
}; // Class TimeIntegrableObject

}  // namespace Kratos.

#endif // KRATOS_TIME_INTEGRABLE_OBJECT_H_INCLUDED defined
