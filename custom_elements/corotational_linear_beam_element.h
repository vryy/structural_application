//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 16 Jun 20 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_COROTATIONAL_LINEAR_BEAM_ELEMENT_H_INCLUDED )
#define  KRATOS_COROTATIONAL_LINEAR_BEAM_ELEMENT_H_INCLUDED


// External includes

// Project includes
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

/**
 */
class CorotationalLinearBeamElement : public Element
{
public:
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    typedef ConstitutiveLaw ConstitutiveLawType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    // Counted pointer of CorotationalLinearBeamElement
    KRATOS_CLASS_POINTER_DEFINITION(CorotationalLinearBeamElement);

    /**
     * Default constructor.
     */
    CorotationalLinearBeamElement();
    CorotationalLinearBeamElement( IndexType NewId, GeometryType::Pointer pGeometry);
    CorotationalLinearBeamElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
     * Destructor.
     */
    ~CorotationalLinearBeamElement() override;

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

    void CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Turn back information as a string.
     * (DEACTIVATED)
     */
    //std::string Info();

    /**
     * Print information about this object.
     * (DEACTIVATED)
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;

    /**
     * Print object's data.
     * (DEACTIVATED)
     */
    //virtual void PrintData(std::ostream& rOStream) const;

protected:


private:

    // double mE, mG, mA, mNU, mI, mk;

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

    double CalculatePsi(const double& E, const double& G,
        const double& L, const double& I, const double& A) const;

}; // Class CorotationalLinearBeamElement

}  // namespace Kratos.

#endif // KRATOS_COROTATIONAL_LINEAR_BEAM_ELEMENT_H_INCLUDED defined
