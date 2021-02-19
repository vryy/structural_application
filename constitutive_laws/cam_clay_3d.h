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
/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 9 Nov 2015 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

#if !defined(KRATOS_CAM_CLAY_3D_H_INCLUDED )
#define  KRATOS_CAM_CLAY_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
//#include "custom_utilities/yield_plot_utility.h"

//#define ENABLE_YIELD_PLOT

namespace Kratos
{

/**
 * Defines a cam-clay constitutive law in 3D space.
 * REF: + Msc thesis, Enes Siljak, Bochum, 2010
 *      + Borja et al, Cam-Clay plasticity, Part 1: Implicit integration of elasto-plastic constitutive relations
 */

class CamClay3D : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;

    #ifdef BOOST_NO_CXX11_CONSTEXPR
    const static double TOL;
    #else
//    constexpr static double TOL = 1.0e-5;
    constexpr static double TOL = 1.0e-8;
    #endif

    const static double unit4thSym3D[][6];

    static const double unit2nd3D[6];

    /**
     * Counted pointer of CamClay3D
     */
    KRATOS_CLASS_POINTER_DEFINITION(CamClay3D);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    CamClay3D();

    virtual  ConstitutiveLaw::Pointer Clone() const
    {
        ConstitutiveLaw::Pointer p_clone( new CamClay3D() );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~CamClay3D();

    /**
     * Operators
     */
    /**
     * Operations
     */
    bool Has( const Variable<int>& rThisVariable );
    bool Has( const Variable<double>& rThisVariable );
    bool Has( const Variable<Vector>& rThisVariable );
    bool Has( const Variable<Matrix>& rThisVariable );

    int& GetValue( const Variable<int>& rThisVariable, int& rValue );
    double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );

    void SetValue( const Variable<int>& rThisVariable, const int& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo );

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& props,
                             const GeometryType& geom,
                             const Vector& ShapeFunctionsValues );

    /**
     * As this constitutive law describes only linear elastic material properties
     * this function is rather useless and in fact does nothing
     */
    void InitializeSolutionStep( const Properties& props,
                                 const GeometryType& geom, //this is just to give the array of nodes
                                 const Vector& ShapeFunctionsValues,
                                 const ProcessInfo& CurrentProcessInfo );

    void InitializeNonLinearIteration( const Properties& props,
                                       const GeometryType& geom, //this is just to give the array of nodes
                                       const Vector& ShapeFunctionsValues,
                                       const ProcessInfo& CurrentProcessInfo );

    void ResetMaterial( const Properties& props,
                        const GeometryType& geom,
                        const Vector& ShapeFunctionsValues );

    void FinalizeNonLinearIteration( const Properties& props,
                                     const GeometryType& geom, //this is just to give the array of nodes
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo );

    void FinalizeSolutionStep( const Properties& props,
                               const GeometryType& geom, //this is just to give the array of nodes
                               const Vector& ShapeFunctionsValues,
                               const ProcessInfo& CurrentProcessInfo );

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    virtual int Check( const Properties& props,
                       const GeometryType& geom,
                       const ProcessInfo& CurrentProcessInfo );

    void CalculateMaterialResponse( const Vector& StrainVector,
                                    const Matrix& DeformationGradient,
                                    Vector& StressVector,
                                    Matrix& AlgorithmicTangent,
                                    const ProcessInfo& CurrentProcessInfo,
                                    const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues,
                                    bool CalculateStresses = true,
                                    int CalculateTangent = true,
                                    bool SaveInternalVariables = true
                                  );

    /**
     * returns the size of the strain vector of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType GetStrainSize()
    {
        return 6;
    }

    /**
     * converts a strain vector styled variable into its form, which the
     * deviatoric parts are no longer multiplied by 2
     */
    //             void Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //virtual String Info() const;
    /**
     * Print information about this object.
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;
    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;

protected:
    /**
     * there are no protected class members
     */
private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw );
        rSerializer.save( "CurrentStress", mCurrentStress );
        rSerializer.save( "LastStress", mLastStress );
        rSerializer.save( "NU", mNU );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw );
        rSerializer.load( "CurrentStress", mCurrentStress );
        rSerializer.load( "LastStress", mLastStress );
        rSerializer.load( "NU", mNU );
    }

    /**
     * Member Variables
     */

    Vector mPrestress;
    Vector mCurrentStress;
    Vector mLastStress;
    Vector mCurrentStrain;
    Vector mLastStrain;
    Vector mLastStrainIncr;
    double mKm, mGm, mNU;
    double mM, mLambda, mKappa, mVoidRatio, mPc, mTheta, mPk, mQk, mPck, mDGamma;

    std::size_t mParentElementId;
    std::size_t mIntegrationPointIndex;

    bool mIsYielded;

    #ifdef ENABLE_YIELD_PLOT
    YieldPlotUtility::Pointer mPlotUtil;
    #endif

    /**
     * Un accessible methods
     */
    int returnMapping(double pTr, double qTr);
    int calculateF(std::vector<double>& rResults, double pTr, double qTr);
    void getFactors(Vector& rResults, double devS);
    int solveG(double& Pc, double pTr);
    double getP();

    void CalculateElasticTangent(Matrix& C, double K, double G);
    Vector& getDeviatoricComp(Vector& devStress, const Vector& stress);
    double getJ2(const Vector& devStress);

    void StressIntegration(const Vector& StrainVector, Vector& StressVector);
    void ComputeTangent(Matrix& AlgorithmicTangent);

    /**
     * Assignment operator.
     */
    //CamClay3D& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
    /**
     * Copy constructor.
     */
    //CamClay3D(const IsotropicPlaneStressWrinklingNew& rOther);
}; // Class CamClay3D


} // namespace Kratos.

#undef ENABLE_YIELD_PLOT

#endif // KRATOS_CAM_CLAY_3D_H_INCLUDED defined
