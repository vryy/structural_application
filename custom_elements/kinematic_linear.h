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
template<typename TNodeType>
class KRATOS_API(STRUCTURAL_APPLICATION) BaseKinematicLinear : public BaseElement<TNodeType>, public BasePrescribedObject<TNodeType>
{

public:
    ///@name Type Definitions
    ///@{

    typedef BaseElement<TNodeType> BaseType;

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

    KRATOS_CLASS_POINTER_DEFINITION( BaseKinematicLinear );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BaseKinematicLinear( IndexType NewId, typename GeometryType::Pointer pGeometry );
    BaseKinematicLinear( IndexType NewId, typename GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties );

    /// Destructor.
    ~BaseKinematicLinear() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    IntegrationMethod GetIntegrationMethod() const override;

    typename BaseType::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties ) const override;

    typename BaseType::Pointer Create( IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties ) const override;

    typename BaseType::Pointer Create( IndexType NewId, std::vector<typename GeometryType::Pointer> pGeom, typename PropertiesType::Pointer pProperties ) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void ResetConstitutiveLaw() override;

    void RewindConstitutiveLaw() override;

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

    void AddInertiaForces(VectorType& rRightHandSideVector, DataType coeff, const ProcessInfo& rCurrentProcessInfo) override;

    void AddDampingForces(VectorType& rRightHandSideVector, DataType coeff, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalAccelerationContribution(MatrixType& rMassMatrix, MatrixType& rMassInducedStiffnessMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalVelocityContribution(MatrixType& rDampMatrix, MatrixType& rDampInducedStiffnessMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

    void FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo ) override;

    void InitializeSolutionStep( const ProcessInfo& CurrentProcessInfo ) override;

    void InitializeNonLinearIteration( const ProcessInfo& CurrentProcessInfo ) override;

    void FinalizeNonLinearIteration( const ProcessInfo& CurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints( const Variable<MatrixType>& rVariable, std::vector<MatrixType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints( const Variable<VectorType>& rVariable, std::vector<VectorType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<DataType, 3> >& rVariable, std::vector<array_1d<DataType, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints( const Variable<DataType>& rVariable, std::vector<DataType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    #ifndef SD_APP_FORWARD_COMPATIBILITY
    void CalculateOnIntegrationPoints( const Variable<std::string>& rVariable, std::vector<std::string>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;
    #endif

    void CalculateOnIntegrationPoints( const Variable<typename ConstitutiveLawType::Pointer>& rVariable, std::vector<typename ConstitutiveLawType::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<DataType>& rVariable, const std::vector<DataType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<int>& rVariable, const std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<bool>& rVariable, const std::vector<bool>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<MatrixType>& rVariable, const std::vector<MatrixType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<VectorType>& rVariable, const std::vector<VectorType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<array_1d<DataType, 3> >& rVariable, const std::vector<array_1d<DataType, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints( const Variable<typename ConstitutiveLawType::Pointer>& rVariable, const std::vector<typename ConstitutiveLawType::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

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
        std::stringstream ss;
        ss << "KinematicLinear<" << DataTypeToString<DataType>::Get() << ">";
        return ss.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << " #" << this->Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
        const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();
        rOStream << "mConstitutiveLawVector::size = " << mConstitutiveLawVector.size() << std::endl;
        rOStream << "IntegrationMethod: " << static_cast<std::underlying_type<IntegrationMethod>::type>(ThisIntegrationMethod);
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
    std::vector<typename ConstitutiveLawType::Pointer> mConstitutiveLawVector;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void InitializeMaterial(const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateBoperator( MatrixType& B_Operator, const Vector& N, const MatrixType& DN_DX ) const;

    virtual void CalculateBBaroperator( MatrixType& B_Operator, const MatrixType& DN_DX, const MatrixType& Bdil_bar ) const;

    virtual unsigned int GetStrainSize( const unsigned int dim ) const
    {
        return dim * (dim + 1) / 2;
    }

    virtual ValueType GetIntegrationWeight( ValueType Weight, const Vector& N ) const
    {
        return Weight;
    }

    virtual void CalculateJacobian( typename GeometryType::JacobiansType& J, const IntegrationMethod ThisIntegrationMethod ) const;

    /** K += weight*Btrans*D*B */
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag );

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    BaseKinematicLinear() {}

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save( "mInitialDisp", mInitialDisp );
        rSerializer.save( "mConstitutiveLawVector", mConstitutiveLawVector );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load( "mInitialDisp", mInitialDisp );
        int tmp;
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

    void CalculateNumericalStiffness(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo, const ValueType epsilon);

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    void AddBodyForcesToRHS( VectorType& R, const Vector& N_DISP, DataType Weight, DataType detJ ) const;

    void CalculateAndAdd_ExtForceContribution( const Vector& N, const ProcessInfo& CurrentProcessInfo,
                                               const VectorType& BodyForce, VectorType& rRightHandSideVector,
                                               DataType Weight, DataType detJ) const;

    void AddInternalForcesToRHS( VectorType& R, const MatrixType& B_Operator, VectorType& StressVector, DataType Weight, DataType detJ ) const;

    void CalculateStiffnesMatrix( MatrixType& K, const MatrixType& tan_C, const MatrixType& B_Operator, DataType Weight, DataType detJ ) const;

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
    //BaseKinematicLinear& operator=(const BaseKinematicLinear& rOther);

    /// Copy constructor.
    //BaseKinematicLinear(const BaseKinematicLinear& rOther);

    ///@}

}; // Class BaseKinematicLinear

///@}
///@name Type Definitions
///@{

typedef BaseKinematicLinear<RealNode> KinematicLinear;
typedef BaseKinematicLinear<ComplexNode> ComplexKinematicLinear;
typedef BaseKinematicLinear<GComplexNode> GComplexKinematicLinear;

///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_KINEMATIC_LINEAR_INCLUDED defined
