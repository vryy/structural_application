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
 * In this element, the strain measure is Green-Lagrange strain and
 * the stress measure is Kirchhoff stress.
 * Because the compatibility of Infinitesimal strain and Green-Lagrange strain,
 * the constitutive law supporting Infinitesimal strain (i.e. linear elastic) can
 * also be used
 * Reference: A. Bower, Applier Mechanics of Solids
 * TODO: this element is incomplete
 */
class FiniteStrain : public Element
{
public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of FiniteStrain
    KRATOS_CLASS_POINTER_DEFINITION(FiniteStrain);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FiniteStrain(IndexType NewId, GeometryType::Pointer pGeometry);
    FiniteStrain(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~FiniteStrain();

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
    IntegrationMethod GetIntegrationMethod() const;

    virtual Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    virtual Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const;

    void Initialize(const ProcessInfo& rCurrentProcessInfo);

    void ResetConstitutiveLaw();

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const;

    void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo);

    void InitializeNonLinearIteration( const ProcessInfo& CurrentProcessInfo );

    void FinalizeNonLinearIteration( const ProcessInfo& CurrentProcessInfo );

    void FinalizeSolutionStep(const ProcessInfo& CurrentProcessInfo);

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo);

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo);

    void SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValuesVector(Vector& values, int Step = 0) const;
    void GetFirstDerivativesVector(Vector& values, int Step = 0) const;
    void GetSecondDerivativesVector(Vector& values, int Step = 0) const;

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
    std::string Info() const final
    {
        return "FiniteStrain";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info() << " #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
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

    ///@}
    ///@name Protected Operators
    ///@{
    FiniteStrain() : Element()
    {
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
    /**
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;
    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    double mTotalDomainInitialSize;
    bool mIsInitialized;

    ///@}
    ///@name Private Operators
    ///@{

//    void CalculateAndAddKm(
//        MatrixType& K,
//        Matrix& B,
//        Matrix& D,
//        double weight);

    /**
     * Calculation of the Geometric Stiffness Matrix. Kg = dB * S
     */
    void CalculateAndAddKg(
        MatrixType& K,
        const Matrix& DN_DX,
        const Vector& StressVector,
        const double& weight
    );

    void CalculateBodyForces(
        Vector& BodyForce,
        const ProcessInfo& CurrentProcessInfo
    );

    void InitializeVariables();

    virtual void InitializeMaterial();

    double CalculateIntegrationWeight(const double& GaussPointWeight, const double& DetJ0);

    void AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, const double& Weight );

    void CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        const Vector& BodyForce,
        VectorType& mResidualVector,
        const double& weight
    );

    void CalculateStrain(const Matrix& C,
                         Vector& StrainVector);

    // void CalculateB(Matrix& B,
    //                 const Matrix& F,
    //                 const Matrix& DN_DX,
    //                 unsigned int StrainSize);

    void CalculateB( Matrix& B_Operator, const Matrix& DN_DX );

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



    virtual void save(Serializer& rSerializer) const;
//        {
//            rSerializer.save("Name", "FiniteStrain");
//            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
//        }

    virtual void load(Serializer& rSerializer);
//        {
//            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
//        }



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
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
FiniteStrain& rThis);
 */
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
const FiniteStrain& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}*/
///@}

} // namespace Kratos.
#endif // KRATOS_FINITE_STRAIN_ELEMENT_H_INCLUDED  defined
