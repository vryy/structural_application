//
//   Project Name:        KratosStructuralApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Jun 20 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_TIMOSHENKO_LINEAR_BEAM_ELEMENT_H_INCLUDED )
#define  KRATOS_TIMOSHENKO_LINEAR_BEAM_ELEMENT_H_INCLUDED


// External includes

// Project includes
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

/**
 * Another implementation of Timoshenko beam element
 */
class TimoshenkoLinearBeamElement : public Element
{
public:
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    typedef ConstitutiveLaw ConstitutiveLawType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    // Counted pointer of TimoshenkoLinearBeamElement
    KRATOS_CLASS_POINTER_DEFINITION(TimoshenkoLinearBeamElement);

    /**
     * Default constructor.
     */
    TimoshenkoLinearBeamElement();
    TimoshenkoLinearBeamElement( IndexType NewId, GeometryType::Pointer pGeometry);
    TimoshenkoLinearBeamElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
     * Destructor.
     */
    ~TimoshenkoLinearBeamElement() override;

    /**
     * Operations.
     */

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix( MatrixType& rMassMatrix,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix( MatrixType& rDampMatrix,
                                 const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList( DofsVectorType& ElementalDofList,
                     const ProcessInfo& CurrentProcessInfo) const override;

    void CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& Output, const ProcessInfo& rCurrentProcessInfo ) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "TimoshenkoLinearBeamElement";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TimoshenkoLinearBeamElement #" << Id();
    }

private:

    Matrix mInitialDisp;
    Matrix mInitialRot;

    friend class Serializer;

    void save ( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Element )
    }

    void load ( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Element )
    }

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       const bool& CalculateStiffnessMatrixFlag,
                       const bool& CalculateResidualVectorFlag);

    void GetCurrentDisplacement( Vector& Disp ) const;

    void CalculateInitialLocalCS(Matrix& transformation_matrix) const;

    void CreateElementStiffnessMatrix_Material(Matrix& local_stiffness_matrix) const;

}; // Class TimoshenkoLinearBeamElement

}  // namespace Kratos.

#endif // KRATOS_TIMOSHENKO_LINEAR_BEAM_ELEMENT_H_INCLUDED defined
