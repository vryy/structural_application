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
//   Modified by:         $Author: janosch $
//   Date:                $Date: 2009-01-14 09:30:38 $
//   Modified by:         $Author: hbui $
//   Date:                $Date: 2013-02-22 16:16:48 $
//   Modified by:         $Author: hbui $
//   Date:                $Date: 2023-07-24 16:16:48 $
//
//


#if !defined(KRATOS_KINEMATIC_LINEAR_INCLUDED )
#define  KRATOS_KINEMATIC_LINEAR_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/serializer.h"
#include "includes/element.h"
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

/// Short class definition.
/** Detail class definition.
 * Define a small strain element with strain measure as Infinitesimal strain and stress measure as Cauchy stress.
 * The inertial and viscous forces are assumed linear.
 * This element is designed to always compute Jacobian in the reference configuration. It is not influenced
 * by the MoveMeshFlag.
 */
class KRATOS_API(STRUCTURAL_APPLICATION) KinematicLinear : public Element, public PrescribedObject
{

public:
    ///@name Type Definitions
    ///@{
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef ConstitutiveLaw ConstitutiveLawType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    KRATOS_CLASS_POINTER_DEFINITION( KinematicLinear );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry );
    KinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~KinematicLinear();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    IntegrationMethod GetIntegrationMethod() const override;

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const override;

    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void ResetConstitutiveLaw() override;

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void ApplyPrescribedDofs(const MatrixType& LHS_Contribution, VectorType& RHS_Constribution, const ProcessInfo& CurrentProcessInfo) const override;

    void ComputePrescribedForces(const MatrixType& LHS_Contribution, VectorType& Force, const ProcessInfo& CurrentProcessInfo) const override;

    void CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const override;

    void GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo ) const override;

    void MassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void DampMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    ///@brief Routines to enable the element to use with nonlinear mass damping time integration scheme
    ///@{

    void AddInertiaForces(VectorType& rRightHandSideVector, double coeff, const ProcessInfo& rCurrentProcessInfo) override;

    void AddDampingForces(VectorType& rRightHandSideVector, double coeff, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalAccelerationContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

    void FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo ) override;

    void InitializeSolutionStep( const ProcessInfo& CurrentProcessInfo ) override;

    void InitializeNonLinearIteration( const ProcessInfo& CurrentProcessInfo ) override;

    void FinalizeNonLinearIteration( const ProcessInfo& CurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints( const Variable<MatrixType>& rVariable, std::vector<MatrixType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints( const Variable<VectorType>& rVariable, std::vector<VectorType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    #ifndef SD_APP_FORWARD_COMPATIBILITY
    void CalculateOnIntegrationPoints( const Variable<std::string>& rVariable, std::vector<std::string>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;
    #endif

    void CalculateOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, std::vector<ConstitutiveLaw::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<double>& rVariable, const std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<int>& rVariable, const std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<bool>& rVariable, const std::vector<bool>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<MatrixType>& rVariable, const std::vector<MatrixType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<VectorType>& rVariable, const std::vector<VectorType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, const std::vector<ConstitutiveLaw::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void GetValuesVector( VectorType& values, int Step = 0 ) const override;

    void GetFirstDerivativesVector( VectorType& values, int Step = 0 ) const override;

    void GetSecondDerivativesVector( VectorType& values, int Step = 0 ) const override;

    int Check( const ProcessInfo& rCurrentProcessInfo ) const override;

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
        return "KinematicLinear";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << " #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        Element::PrintData(rOStream);
        rOStream << "mConstitutiveLawVector::size = " << mConstitutiveLawVector.size() << std::endl;
        rOStream << "IntegrationMethod: " << static_cast<std::underlying_type<IntegrationMethod>::type>(mThisIntegrationMethod);
    }

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

    MatrixType mInitialDisp;
    IntegrationMethod mThisIntegrationMethod;
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void InitializeMaterial(const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateBoperator( MatrixType& B_Operator, const VectorType& N, const MatrixType& DN_DX ) const;

    virtual void CalculateBBaroperator( MatrixType& B_Operator, const MatrixType& DN_DX, const MatrixType& Bdil_bar ) const;

    virtual unsigned int GetStrainSize( const unsigned int& dim ) const
    {
        return dim * (dim + 1) / 2;
    }

    virtual double GetIntegrationWeight( double Weight, const VectorType& N ) const
    {
        return Weight;
    }

    virtual void CalculateJacobian( GeometryType::JacobiansType& J ) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    KinematicLinear() {}

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer,  Element );
        rSerializer.save( "mInitialDisp", mInitialDisp );
        rSerializer.save( "mThisIntegrationMethod", static_cast<int>(mThisIntegrationMethod) );
        rSerializer.save( "mConstitutiveLawVector", mConstitutiveLawVector );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer,  Element );
        rSerializer.load( "mInitialDisp", mInitialDisp );
        int tmp;
        rSerializer.load( "mThisIntegrationMethod", tmp );
        mThisIntegrationMethod = static_cast<IntegrationMethod>(tmp);
        rSerializer.load( "mConstitutiveLawVector", mConstitutiveLawVector );
    }

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

    //************************************************************************************
    //************************************************************************************
    //DEBUGGING
    //************************************************************************************
    //************************************************************************************

    void CalculateNumericalStiffness(Matrix& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo, const double epsilon);

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{
    /** K += weight*Btrans*D*B */
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //CALCULATE FORCEVECTORS DISPLACEMENT

    void AddBodyForcesToRHS( VectorType& R, const VectorType& N_DISP, double Weight, double detJ ) const;

    void CalculateAndAdd_ExtForceContribution( const VectorType& N, const ProcessInfo& CurrentProcessInfo,
                                               const VectorType& BodyForce, VectorType& rRightHandSideVector,
                                               double weight, double detJ) const;

    void AddInternalForcesToRHS( VectorType& R, const MatrixType& B_Operator, VectorType& StressVector, double Weight, double detJ ) const;

    void CalculateStiffnesMatrix( MatrixType& K, const MatrixType& tan_C, const MatrixType& B_Operator, double Weight, double detJ ) const;

    void CalculateStressAndTangentialStiffness( VectorType& StressVector, MatrixType& tanC_U,
                                                VectorType& StrainVector, const MatrixType& B_Operator,
                                                int PointNumber, const ProcessInfo& CurrentProcessInfo ) const;

    void CalculateStrain( const MatrixType& B, const MatrixType& Displacements, VectorType& StrainVector ) const;

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //KinematicLinear& operator=(const KinematicLinear& rOther);

    /// Copy constructor.
    //KinematicLinear(const KinematicLinear& rOther);


    ///@}

}; // Class KinematicLinear

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                   KinematicLinear& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                   const KinematicLinear& rThis)
           {
                   rThis.PrintInfo(rOStream);
                   rOStream << std::endl;
                   rThis.PrintData(rOStream);

                   return rOStream;
}*/
///@}

}  // namespace Kratos.

#endif // KRATOS_KINEMATIC_LINEAR_INCLUDED defined


