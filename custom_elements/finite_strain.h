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
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Aug 2021 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_FINITE_STRAIN_ELEMENT_H_INCLUDED )
#define  KRATOS_FINITE_STRAIN_ELEMENT_H_INCLUDED



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

/**
 * Implements a finite strain definition for structural analysis.
 * This works for arbitrary geometries in 2D and 3D
 * In this element, the stress measure is Cauchy and the strain measure is Hencky.
 * Reference: De Souza Neto, Computational Plasticity, Box 14.3
 */
class KRATOS_API(STRUCTURAL_APPLICATION) FiniteStrain : public Element, public PrescribedObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of FiniteStrain
    KRATOS_CLASS_POINTER_DEFINITION(FiniteStrain);

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

    typedef typename BaseType::IdentityMatrixType IdentityMatrixType;

    typedef typename BaseType::GeometryType GeometryType;

    typedef typename BaseType::PropertiesType PropertiesType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FiniteStrain(IndexType NewId, GeometryType::Pointer pGeometry);
    FiniteStrain(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    ~FiniteStrain() override;

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

    void SetValuesOnIntegrationPoints( const Variable<DataType>& rVariable, const std::vector<DataType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<int>& rVariable, const std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints(const Variable<VectorType>& rVariable, std::vector<VectorType>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void SetValuesOnIntegrationPoints(const Variable<MatrixType>& rVariable, std::vector<MatrixType>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void SetValuesOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, const std::vector<ConstitutiveLaw::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<VectorType>& rVariable, std::vector<VectorType>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<MatrixType>& rVariable, std::vector<MatrixType>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(VectorType& values, int Step = 0) const override;
    void GetFirstDerivativesVector(VectorType& values, int Step = 0) const override;
    void GetSecondDerivativesVector(VectorType& values, int Step = 0) const override;

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
        return "FiniteStrain";
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

    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    ///@}
    ///@name Protected Operators
    ///@{
    FiniteStrain() : Element()
    {
    }

    /// Get the size of the strain vector
    virtual unsigned int GetStrainSize( unsigned int dim ) const
    {
        return dim * (dim + 1) / 2;
    }

    /// Get the size of unsymmetric second order tensor
    virtual unsigned int GetFSize( unsigned int dim ) const
    {
        return dim;
    }

    /// Get the size of unsymmetric second order tensor
    virtual unsigned int GetGSize( unsigned int dim ) const
    {
        return dim * dim;
    }

    /// Calculate the B operator
    virtual void CalculateB( MatrixType& B_Operator, const Vector& N, const MatrixType& DN_DX, const MatrixType& CurrentDisp ) const;

    /// Calculate the G operator
    virtual void CalculateG( MatrixType& G_Operator, const Vector& N, const MatrixType& DN_DX ) const;

    /// Calculate the G operator
    virtual void CalculateG( MatrixType& G_Operator, const Vector& N, const MatrixType& DN_DX, const MatrixType& CurrentDisp ) const;

    /// Calculate the deformation gradient
    virtual void CalculateF( MatrixType& F, const MatrixType& G_Operator, const MatrixType& CurrentDisp ) const;

    /// Get the integration weight
    virtual double GetIntegrationWeight( double Weight, const Vector& N, const MatrixType& CurrentDisp ) const
    {
        return Weight;
    }

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    virtual void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo,
                              bool CalculateStiffnessMatrixFlag,
                              bool CalculateResidualVectorFlag);

    void AddBodyForcesToRHS( VectorType& R, const Vector& N_DISP, const double Weight ) const;

    void CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        const VectorType& BodyForce,
        VectorType& mResidualVector,
        const double weight
    ) const;

    /// Calculate Almansi strain, providing B
    virtual void CalculateStrain( const MatrixType& B, VectorType& StrainVector ) const;

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
    /*  static MatrixType msB;
    static MatrixType msF;
    static MatrixType msD;
    static MatrixType msC;
    static Vector msStrainVector;
    static Vector msStressVector;
    static MatrixType msDN_DX;
     */
    ///@}
    ///@name Member Variables
    ///@{

    double mTotalDomainInitialSize;
    bool mIsInitialized;

    ///@}
    ///@name Private Operators
    ///@{

    void CalculateBodyForces(
        VectorType& BodyForce,
        const ProcessInfo& CurrentProcessInfo
    );

    void InitializeVariables();

    virtual void InitializeMaterial(const ProcessInfo& CurrentProcessInfo);

    double CalculateIntegrationWeight(const double GaussPointWeight, const double DetJ0);

    // void CalculateB(MatrixType& B,
    //                 const MatrixType& F,
    //                 const MatrixType& DN_DX,
    //                 unsigned int StrainSize);

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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //FiniteStrain& operator=(const FiniteStrain& rOther);
    /// Copy constructor.
    //FiniteStrain(const FiniteStrain& rOther);
    ///@}

}; // Class FiniteStrain

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_FINITE_STRAIN_ELEMENT_H_INCLUDED  defined
