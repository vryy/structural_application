/*
*/
/* *********************************************************
*
*   Last Modified by:    $Author: hbui$
*   Date:                $Date: 6 Feb 2024 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(MEAN_DISPLACEMENT_CONSTRAINT )
#define  MEAN_DISPLACEMENT_CONSTRAINT



// System includes


// External includes


// Project includes
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

/**
 * This condition can be used to enforce the zero average displacement along an edge using global Lagrange multiplier.
 * It is useful to circumvent the point-wise boundary condition, which can lead to stress concentration.
 * Reference:
 * - https://math.stackexchange.com/questions/2112280/zero-average-in-solution-to-laplace-problem-with-lagrange-multipliers
 */
template<int IDimIndex>
class KRATOS_API(STRUCTURAL_APPLICATION) MeanDisplacementConstraint : public Condition
{

public:
    typedef Condition BaseType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef BaseType::MatrixType LHS_ContributionType;

    typedef Condition::GeometryType::Pointer PointerGeometryType;

    // Counted pointer of PointPointContactLink
    KRATOS_CLASS_POINTER_DEFINITION(MeanDisplacementConstraint);

    /**
     * Default constructor.
     */
    MeanDisplacementConstraint( IndexType NewId, GeometryType::Pointer pGeometry);

    MeanDisplacementConstraint( IndexType NewId, GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties);

    /**
     * Destructor.
     */
    virtual ~MeanDisplacementConstraint();

    /**
     * Operations.
     */
    Condition::Pointer Create( IndexType NewId,
                               NodesArrayType const& ThisNodes,
                               PropertiesType::Pointer pProperties) const final;

    Condition::Pointer Create( IndexType NewId,
                               GeometryType::Pointer pGeom,
                               PropertiesType::Pointer pProperties) const final;

    /**
     * Calculates the local system contributions for this contact element
     */
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo) final;

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo) const final;

    void GetDofList( DofsVectorType& ConditionalDofList,
                     const ProcessInfo& CurrentProcessInfo) const final;

    std::string Info() const final
    {
        std::stringstream ss;
        ss << "MeanDisplacementConstraint<" << IDimIndex << ">";
        return ss.str();
    }

private:

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag);

    friend class Serializer;

    // A private default constructor necessary for serialization
    MeanDisplacementConstraint() {};

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

}; // Class MeanDisplacementConstraint

}  // namespace Kratos.

#endif // MEAN_DISPLACEMENT_CONSTRAINT  defined
