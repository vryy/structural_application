
//   Project Name:        Kratos
//   Last Modified by:    $Author: nelson $
//   Date:                $Date: 2008-12-10 11:10:16 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_BEAM_ELEMENT_INCLUDED )
#define  KRATOS_BEAM_ELEMENT_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

class BeamElement : public Element
{

    typedef GeometryData::IntegrationMethod IntegrationMethod;

private:
    ///@name Static Member Variables

    Matrix mInitialDisp;
    Matrix mInitialRot; // mInitialRot and mInitialDisp are used to account for the stress-free reactivation
    Vector mPreForces; // mPreForces is used to prescribe the PRESTRESS for beam element
    Vector mCurrentForces;

    double mArea;                            // Area de la seccion tranversal de la viga.
    double mInertia_x;                       // Momento de Inercia alredor del eje Ix local.
    double mInertia_y;                       // Momento de Inercia alrededor del eje Iy local.
    double mInertia_Polar;                   // Momento Polar de Inercia
    double mlength;                          // Longitud del Elemento.


    void CalculateSectionProperties();

    void CalculateLocalMatrix(Matrix& LocalMatrix);

    void CalculateTransformationMatrix(Matrix& Rotation);

    void CalculateBodyForce(const Matrix& Rotation, Vector& LocalBody, Vector& GlobalBody);

    void CalculateLocalNodalStress(Vector& Stress);

    double CalculateInternalAxil(   const double& Ao, const double& Load, const double& X);
    double CalculateInternalShear(  const double& Vo, const double& Load, const double& X);
    double CalculateInternalMoment( const double& Mo, const double& Vo,   const double& Load, const double& X);

    void CalculateDistributedBodyForce(const int Direction, Vector& Load);

public:

    KRATOS_CLASS_POINTER_DEFINITION(BeamElement);

    /// Default constructor.
    BeamElement(IndexType NewId, GeometryType::Pointer pGeometry);
    BeamElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Constructor for 2-node beam
    BeamElement(IndexType NewId, GeometryType::PointType::Pointer pNode1,
        GeometryType::PointType::Pointer pNode2, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~BeamElement();

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const final;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const final;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) final;

    void ResetConstitutiveLaw() final;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) final;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const final;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const final;

    void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo) final;

    void FinalizeSolutionStep(const ProcessInfo& CurrentProcessInfo) final;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateDampingMatrix(MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateRHS(Vector& rRightHandSideVector);

    void CalculateLHS(Matrix& rLeftHandSideMatrix);

    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo,
                      bool CalculateStiffnessMatrixFlag,
                      bool CalculateResidualVectorFlag);

    void GetValuesVector(Vector& values, int Step) const final;
    void GetFirstDerivativesVector(Vector& values, int Step) const final;
    void GetSecondDerivativesVector(Vector& values, int Step) const final;

    void SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      const std::vector<Vector>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
                                      std::vector< array_1d<double,3> >& Output,
                                      const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>& Output,
                                      const ProcessInfo& rCurrentProcessInfo) final;

    IntegrationMethod GetIntegrationMethod() const final;

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const final;

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        std::stringstream buffer;
        buffer << "BeamElement #" << Id();
        return buffer.str();
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << "BeamElement #" << Id();
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const final
    {
        Element::PrintData(rOStream);
    }

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    BeamElement() {};


    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
    }

}; // Class BeamElement

} // Namespace Kratos.


#endif

