/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
giang.bui@rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hurga $
//   Date:                $Date: 2007-10-18 16:23:41 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_TOTAL_LAGRANGIAN_ELEMENT_H_INCLUDED )
#define  KRATOS_TOTAL_LAGRANGIAN_ELEMENT_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "custom_elements/prescribed_object.h"


namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Total Lagrangian element for 2D and 3D geometries.

/**
 * Implements a total Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 2D and 3D
 * In this element, the strain measure is Green-Lagrange strain and
 * the stress measure is 2nd Piola Kirchhoff stress.
 * Because the compatibility of Infinitesimal strain and Green-Lagrange strain,
 * the constitutive law supporting Infinitesimal strain (i.e. linear elastic) can
 * also be used
 * It is noted that, this element is only compatible with symmetric constitutive law (i.e C_ijkl=C_jikl=C_ijlk=C_klij)
 * Reference: D. Kuhl, Computational Dynamics lecture note
 */
class KRATOS_API(STRUCTURAL_APPLICATION) TotalLagrangian : public Element, public PrescribedObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of TotalLagrangian
    KRATOS_CLASS_POINTER_DEFINITION(TotalLagrangian);

    typedef Element BaseType;

    typedef typename BaseType::ElementType ElementType;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef typename BaseType::ConstitutiveLawType ConstitutiveLawType;

    typedef typename ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;

    typedef typename BaseType::DofsVectorType DofsVectorType;

    typedef typename BaseType::IndexType IndexType;

    typedef typename BaseType::DataType DataType;

    typedef typename BaseType::ValueType ValueType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::ZeroVectorType ZeroVectorType;

    typedef typename BaseType::MatrixType MatrixType;

    typedef typename BaseType::ZeroMatrixType ZeroMatrixType;

    typedef typename BaseType::GeometryType GeometryType;

    typedef typename BaseType::PropertiesType PropertiesType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /// local Jacobian container
    struct MyJacobians
    {
        GeometryType::JacobiansType J0Values, JnValues, JValues;
        GeometryType::JacobiansType *J0p, *Jnp, *Jp;

        MyJacobians() : J0p(nullptr), Jnp(nullptr), Jp(nullptr) {}

        const GeometryType::JacobiansType& J0() const {return *J0p;}
        const GeometryType::JacobiansType& Jn() const {return *Jnp;}
        const GeometryType::JacobiansType& J()  const {return *Jp;}

        void Check() const
        {
            if (J0p == nullptr || Jnp == nullptr || Jp == nullptr)
                KRATOS_ERROR << "The Jacobians are not fully defined";
        }
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TotalLagrangian(IndexType NewId, GeometryType::Pointer pGeometry);
    TotalLagrangian(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    ~TotalLagrangian() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    IntegrationMethod GetIntegrationMethod() const override;

    typename ElementType::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    typename ElementType::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void ResetConstitutiveLaw() override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void ApplyPrescribedDofs(const MatrixType& LHS_Contribution, VectorType& RHS_Constribution, const ProcessInfo& CurrentProcessInfo) const override;

    void ComputePrescribedForces(const MatrixType& LHS_Contribution, VectorType& Force, const ProcessInfo& CurrentProcessInfo) const override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo) override;

    void InitializeNonLinearIteration( const ProcessInfo& CurrentProcessInfo ) override;

    void FinalizeNonLinearIteration( const ProcessInfo& CurrentProcessInfo ) override;

    void FinalizeSolutionStep(const ProcessInfo& CurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<DataType>& rVariable, std::vector<DataType>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<DataType, 3> >& rVariable, std::vector<array_1d<DataType, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& values, int Step = 0) const override;
    void GetFirstDerivativesVector(Vector& values, int Step = 0) const override;
    void GetSecondDerivativesVector(Vector& values, int Step = 0) const override;

    void Calculate(const Variable<DataType>& rVariable, DataType& Output, const ProcessInfo& rCurrentProcessInfo);

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) const final;

    //std::string Info() const;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "TotalLagrangian";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << " #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {}

    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    /**
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;

    ///@}
    ///@name Protected Operators
    ///@{
    TotalLagrangian() : Element()
    {
    }

    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    virtual void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo,
                              bool CalculateStiffnessMatrixFlag,
                              bool CalculateResidualVectorFlag);

    /// Compute the Jacobian in reference and current configuration
    virtual void CalculateJacobians( MyJacobians& Jacobians, const Matrix& DeltaPosition, const Matrix& CurrentDisp ) const;
    void CalculateJacobians( MatrixType& J0, MatrixType& J, const CoordinatesArrayType& rCoordinates,
            const Matrix& DeltaPosition, const Matrix& CurrentDisp ) const;

    /**
     * Calculation of the Geometric Stiffness Matrix. Kg = dB * S
     */
    virtual void CalculateAndAddKg(
        MatrixType& K,
        const Vector& N,
        const Matrix& DN_DX,
        const Vector& StressVector,
        DataType Weight
    ) const;

    void CalculateBodyForces(
        Vector& BodyForce,
        const ProcessInfo& CurrentProcessInfo
    ) const;

    void AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, const DataType& Weight ) const;

    void CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        const Vector& BodyForce,
        VectorType& mResidualVector,
        DataType Weight
    ) const;

    virtual unsigned int GetStrainSize( unsigned int dim ) const
    {
        return dim * (dim + 1) / 2;
    }

    virtual unsigned int GetFSize( unsigned int dim ) const
    {
        return dim;
    }

    // Calculate Green-Lagrange strain, providing C
    virtual void CalculateStrain( const Matrix& C,
                                  Vector& StrainVector ) const;

    virtual void CalculateF( Matrix& F,
                             const Vector& N,
                             const Matrix& DN_DX,
                             const Matrix& CurrentDisp ) const;

    virtual void CalculateB( Matrix& B,
                             const Matrix& F,
                             const Vector& N,
                             const Matrix& DN_DX ) const
    {
        this->CalculateB(B, F, DN_DX);
    }

    virtual DataType GetIntegrationWeight( DataType Weight, const Vector& N ) const
    {
        return Weight;
    }

    void Comprobate_State_Vector(Vector& Result) const;

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    /*  static Matrix msB;
    static Matrix msF;
    static Matrix msD;
    static Matrix msC;
    static Vector msStrainVector;
    static Vector msStressVector;
    static Matrix msDN_DX;
     */
    ///@}
    ///@name Member Variables
    ///@{

    DataType mTotalDomainInitialSize;
    bool mIsInitialized;

    ///@}
    ///@name Private Operators
    ///@{

//    void CalculateAndAddKm(
//        MatrixType& K,
//        Matrix& B,
//        Matrix& D,
//        DataType weight);

    void InitializeVariables();

    virtual void InitializeMaterial( const ProcessInfo& CurrentProcessInfo );

    void CalculateB(Matrix& B,
                    const Matrix& F,
                    const Matrix& DN_DX) const;

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;
//        {
//            rSerializer.save("Name", "TotalLagrangian");
//            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
//        }

    void load(Serializer& rSerializer) override;
//        {
//            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
//        }



    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //TotalLagrangian& operator=(const TotalLagrangian& rOther);
    /// Copy constructor.
    //TotalLagrangian(const TotalLagrangian& rOther);
    ///@}

}; // Class TotalLagrangian

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
TotalLagrangian& rThis);
 */
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
const TotalLagrangian& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}*/
///@}

} // namespace Kratos.
#endif // KRATOS_TOTAL_LAGRANGIAN_ELEMENT_H_INCLUDED  defined
