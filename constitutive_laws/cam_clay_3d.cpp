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

// System includes
#include <iostream>
#include <cmath>
#include <math.h>

// External includes
#include <boost/math/special_functions/fpclassify.hpp>

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "utilities/kratos_log.h"
#include "constitutive_laws/cam_clay_3d.h"
#include "structural_application_variables.h"

#define DEBUG_CAM_CLAY
#define DEBUG_ELEMENT_ID 1
#define DEBUG_POINT_ID 0
// #define LOG_CAM_CLAY
//#define ENABLE_YIELD_PLOT

#define RETURN_MAPPING_NOT_CONVERGED -100
#define LOCAL_PC_SOLVE_INVALID_INPUT -200
#define LOCAL_PC_SOLVE_NOT_CONVERGED -300
#define LOCAL_PC_SOLVE_NEGATIVE_PTRIAL 1
#define LOCAL_PC_SOLVE_SMALL_INPUT 2

namespace Kratos
{

#ifdef BOOST_NO_CXX11_CONSTEXPR
//const double CamClay3D::TOL = 1.0e-5;
const double CamClay3D::TOL = 1.0e-8;
#endif

const double CamClay3D::unit2nd3D[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};

const double CamClay3D::unit4thSym3D[][6] = {
                {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.5, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.5, 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.5}
            };

/**
 * TO BE TESTED!!!
 */
CamClay3D::CamClay3D()
    : ConstitutiveLaw()
{
    mParentElementId = -1;
    mIntegrationPointIndex = -1;
}

/**
 * TO BE TESTED!!!
 */
CamClay3D::~CamClay3D()
{
}

bool CamClay3D::Has( const Variable<int>& rThisVariable )
{
    return false;
}

bool CamClay3D::Has( const Variable<double>& rThisVariable )
{
    if ( rThisVariable == YOUNG_MODULUS || rThisVariable == POISSON_RATIO )
        return true;

    return false;
}

bool CamClay3D::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == PRESTRESS )
        return true;

    if ( rThisVariable == INSITU_STRESS )
        return true;

    if ( rThisVariable == STRESSES )
        return true;

//    if ( rThisVariable == CURRENT_STRAIN_VECTOR )
//        return true;

    return false;
}

bool CamClay3D::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}


int& CamClay3D::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return rValue;
}

double& CamClay3D::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if ( rThisVariable == POISSON_RATIO ) {
        rValue = mNU;
        return rValue;
    }

    if (rThisVariable == PRESSURE_P) // this follows the soil mechanics's hydrostatic pressure convention
    {
        rValue = (mCurrentStress[0] + mCurrentStress[1] + mCurrentStress[2]) / 3;
        return rValue;
    }

    if (rThisVariable == PRESSURE_Q)
    {
        double p = (mCurrentStress[0] + mCurrentStress[1] + mCurrentStress[2]) / 3;
        double sxx = mCurrentStress[0] - p;
        double syy = mCurrentStress[1] - p;
        double szz = mCurrentStress[2] - p;
        double sxy = mCurrentStress[3];
        double syz = mCurrentStress[4];
        double sxz = mCurrentStress[5];

        rValue = sqrt( 3.0 * ( 0.5*(sxx*sxx + syy*syy + szz*szz) + sxy*sxy + syz*syz + sxz*sxz ) );
        return rValue;
    }

    return rValue;
}

Vector& CamClay3D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == PRESTRESS )
    {
        if(rValue.size() != mPrestress.size())
            rValue.resize(mPrestress.size(), false);
        noalias(rValue) = mPrestress;
        return rValue;
    }

    if ( rThisVariable == STRESSES )
    {
        if(rValue.size() != mCurrentStress.size())
            rValue.resize(mCurrentStress.size(), false);
        noalias(rValue) = -mCurrentStress;
        return rValue;
    }

    if ( rThisVariable == PLASTIC_STRAIN_VECTOR )
    {
        // TODO
        rValue = ZeroVector( 6 );
        return( rValue );
    }

    if ( rThisVariable == INTERNAL_VARIABLES )
    {
        // TODO
        rValue = ZeroVector( 1 );
        return( rValue );
    }

    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}

Matrix& CamClay3D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
        KRATOS_THROW_ERROR( std::logic_error, "Matrix Variable case not considered", "" );
}

void CamClay3D::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
    {
        mParentElementId = rValue;
    }
    if(rThisVariable == INTEGRATION_POINT_INDEX)
    {
        mIntegrationPointIndex = rValue;
    }
}

void CamClay3D::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;

    if ( rThisVariable == PRECONSOLIDATION_PRESSURE_MIN )
    {
        mPck = this->getP();

        if(rValue < 0.0)
            KRATOS_THROW_ERROR(std::logic_error, "Negative PRECONSOLIDATION_PRESSURE is detected", "")

        if(mPck < rValue)
            mPck = rValue;

        mKm = (1.0 + mVoidRatio) * mPck / mKappa;
        mGm = 3.0 * mKm * (1.0 - 2.0 * mNU) / 2.0 / (1.0 + mNU);

        #ifdef DEBUG_CAM_CLAY_PRECONSOLIDATION
        std::cout << "At PRECONSOLIDATION_PRESSURE_MIN, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << ":" << std::endl;
        std::cout << " mVoidRatio: " << mVoidRatio << std::endl;
        std::cout << " this->getP(): " << this->getP() << std::endl;
        std::cout << " mPck: " << mPck << std::endl;
        std::cout << " mKappa: " << mKappa << std::endl;
        std::cout << " mNU: " << mNU << std::endl;
        std::cout << " mKm: " << mKm << std::endl;
        std::cout << " mGm: " << mGm << std::endl;
        #endif
    }

    if ( rThisVariable == PRECONSOLIDATION_PRESSURE_DEF || rThisVariable == PRECONSOLIDATION_PRESSURE )
    {
        if(rValue < 0.0)
            KRATOS_THROW_ERROR(std::logic_error, "Negative PRECONSOLIDATION_PRESSURE is detected", "")

        mPck = rValue;

        mKm = (1.0 + mVoidRatio) * mPck / mKappa;
        mGm = 3.0 * mKm * (1.0 - 2.0 * mNU) / 2.0 / (1.0 + mNU);

        #ifdef DEBUG_CAM_CLAY_PRECONSOLIDATION
        std::cout << "At PRECONSOLIDATION_PRESSURE_DEF, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << ":" << std::endl;
        std::cout << " mVoidRatio: " << mVoidRatio << std::endl;
        std::cout << " mPck: " << mPck << std::endl;
        std::cout << " mKappa: " << mKappa << std::endl;
        std::cout << " mNU: " << mNU << std::endl;
        std::cout << " mKm: " << mKm << std::endl;
        std::cout << " mGm: " << mGm << std::endl;
        #endif
    }

    if ( rThisVariable == CSL_SLOPE )
    {
        mM = rValue;
    }

    if ( rThisVariable == VIRGIN_COMPRESSION_INDEX )
    {
        mLambda = rValue;
    }

    if ( rThisVariable == SWELL_INDEX )
    {
        mKappa = rValue;
    }

    if ( rThisVariable == VOID_RATIO )
    {
        mVoidRatio = rValue;
    }
}

void CamClay3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        noalias(mPrestress) = rValue;
        noalias(mLastStress) = rValue;
        noalias(mCurrentStress) = mLastStress;

        mPk = (rValue[0] + rValue[1] + rValue[2]) / 3;
        Vector devStress(6);
        devStress = this->getDeviatoricComp(devStress, rValue);
        double J2 = this->getJ2( devStress );
        mQk = sqrt(3.0 * J2);
//        mQk = sqrt( 0.5*(pow(rValue[0] - rValue[1], 2) + pow(rValue[1] - rValue[2], 2) + pow(rValue[0] - rValue[2], 2))
//                    + 3.0*(pow(rValue[3], 2) + pow(rValue[4], 2) + pow(rValue[5], 2)) );

        // adjust the stiffness parameters
        mKm = (1.0 + mVoidRatio) * this->getP() / mKappa;
        mGm = 3.0 * mKm * (1.0 - 2.0 * mNU) / 2.0 / (1.0 + mNU);
    }

//    if ( rThisVariable == STRESSES )
//    {
//        noalias(mCurrentStress) = rValue;
//    }

//    if ( rThisVariable == CURRENT_STRAIN_VECTOR )
//    {
//        noalias(mCurrentStrain) = -rValue;
//    }
}

void CamClay3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void CamClay3D::InitializeMaterial( const Properties& props,
        const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
    mPrestress = ZeroVector( 6 );
    mCurrentStress = ZeroVector( 6 );
    mLastStress = ZeroVector( 6 );
    mCurrentStrain = ZeroVector( 6 );
    mLastStrain = ZeroVector( 6 );
    mLastStrainIncr = ZeroVector( 6 );
    mNU = props[POISSON_RATIO];
    mM = props[CSL_SLOPE];
    mLambda = props[VIRGIN_COMPRESSION_INDEX];
    mKappa = props[SWELL_INDEX];
    mVoidRatio = props[VOID_RATIO];
    mTheta = (1 + mVoidRatio) / (mLambda - mKappa);
    mDGamma = 0.0;
    mPk = 0.0;
    mQk = 0.0;
    mPck = 0.0;
    mPc = mPck;
    mIsYielded = false;
    #ifdef ENABLE_YIELD_PLOT
    mPlotUtil = YieldPlotUtility::Pointer(new YieldPlotUtility());
    if(mParentElementId == 1 && mIntegrationPointIndex == 18)
    {
        std::stringstream ss;
        ss << "cam_clay_3d_e" << mParentElementId << "_p" << mIntegrationPointIndex;
        mPlotUtil->SetFileOutput(true, ss.str());
    }
    #endif
}

void CamClay3D::ResetMaterial( const Properties& props,
        const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
    noalias(mCurrentStress) = mPrestress;
    noalias(mLastStress) = mPrestress;
    noalias(mCurrentStrain) = ZeroVector( 6 );
    noalias(mLastStrain) = ZeroVector( 6 );
    noalias(mLastStrainIncr) = ZeroVector( 6 );
    mTheta = (1 + mVoidRatio) / (mLambda - mKappa);
    mPc = mPck;
    mDGamma = 0.0;

    // TODO: check if the prestress point still lies in the yield surface. Otherwise, we have to expand the yield surface accordingly.
    double yield_value = mQk * mQk / (mM * mM) + mPk * (mPk - mPc);
    mIsYielded = (yield_value / fabs(mPc)) > 1.0e-6;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
    {
        std::cout << "At ResetMaterial, material point elem = " << DEBUG_ELEMENT_ID << ", point = " << DEBUG_POINT_ID << ":" << std::endl;
        KRATOS_WATCH(mPrestress)
        KRATOS_WATCH(mPk)
        KRATOS_WATCH(mQk)
        KRATOS_WATCH(mPck)
        KRATOS_WATCH(mM)
        KRATOS_WATCH(yield_value)
        KRATOS_WATCH(yield_value / fabs(mPc))
        KRATOS_WATCH(mIsYielded)
    }
    #endif

    if((mIsYielded == true) && (mPc > 0.0))
    {
        std::cout << "----------------ERROR----------------" << std::endl;
        KRATOS_WATCH(mParentElementId)
        KRATOS_WATCH(mIntegrationPointIndex)
        KRATOS_WATCH(mPrestress)
        KRATOS_WATCH(mPk)
        KRATOS_WATCH(mQk)
        KRATOS_WATCH(mPc)
        KRATOS_WATCH(mM)
        KRATOS_WATCH(yield_value)
        KRATOS_WATCH(yield_value / fabs(mPc))
        KRATOS_THROW_ERROR(std::logic_error, "The prestress point violate the yield criteria. It is not handled for now", "")
    }

    #ifdef ENABLE_YIELD_PLOT
    mPlotUtil->Reset();
    if(mParentElementId == 1 && mIntegrationPointIndex == 18)
    {
        std::stringstream ss;
        ss << "cam_clay_3d_e" << mParentElementId << "_p" << mIntegrationPointIndex;
        mPlotUtil->SetFileOutput(true, ss.str());
    }
    #endif
}

void CamClay3D::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mPc = mPck;
    mDGamma = 0.0;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
    {
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "At InitializeSolutionStep, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << std::endl;
        std::cout << " mCurrentStrain: " << mCurrentStrain << std::endl;
        std::cout << " mLastStrain: " << mLastStrain << std::endl;
        std::cout << " mLastStrainIncr: " << mLastStrainIncr << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif
}

void CamClay3D::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
//    #ifdef DEBUG_CAM_CLAY
//    if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
//    {
//        std::cout << "---------------------------------------------" << std::endl;
//        std::cout << "At InitializeNonLinearIteration, material point elem = "
//                  << mParentElementId << ", int = " << mIntegrationPointIndex << std::endl;
////        std::cout << " mCurrentStrain: " << mCurrentStrain << std::endl;
////        std::cout << " mLastStrain: " << mLastStrain << std::endl;
//        std::cout << "---------------------------------------------" << std::endl;
//    }
//    #endif
}

void CamClay3D::CalculateMaterialResponse( const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo,
        const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        bool CalculateStresses,
        int CalculateTangent,
        bool SaveInternalVariables )
{
    if(CalculateStresses)
        StressIntegration(StrainVector, StressVector);
    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == 902 && mIntegrationPointIndex == 0)
    {
        double norm_stress = norm_2(StressVector);
        if((norm_stress != norm_stress) || std::isnan(norm_stress) /*|| isnan(norm_stress)*/ || boost::math::isnan(norm_stress))
        {
            std::stringstream ss;
            ss << "NaN stress at element " << mParentElementId << ", point " << mIntegrationPointIndex << ": " << StressVector;
            KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
        }
        std::stringstream ss;
        ss << "elem " << mParentElementId << ", point " << mIntegrationPointIndex << ": stress = " << StressVector << ", norm_stress: " << norm_stress;
        std::cout << ss.str() << std::endl;
    }
    #endif
    if(CalculateTangent)
        ComputeTangent(AlgorithmicTangent);
}

void CamClay3D::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void CamClay3D::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mKm = (1.0 + mVoidRatio) * this->getP() / mKappa;
    mGm = 3.0 * mKm * (1.0 - 2.0 * mNU) / 2.0 / (1.0 + mNU);

    noalias(mLastStrain) = mCurrentStrain;
    noalias(mLastStress) = mCurrentStress;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
    {
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "At FinalizeSolutionStep, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << std::endl;
//        std::cout << " mCurrentStrain: " << mCurrentStrain << std::endl;
//        std::cout << " mLastStrain: " << mLastStrain << std::endl;
//        std::cout << " mKm: " << mKm << std::endl;
//        std::cout << " mGm: " << mGm << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif
}

//**********************************************************************
int CamClay3D::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    // NO CHECK AT ALL
    // TODO SHALL WE CHECK?

    return 0;

    KRATOS_CATCH( "" );
}

void CamClay3D::StressIntegration(const Vector& StrainVector, Vector& StressVector)
{
    ///////////////////////////////////////////////////////////////////////
    // stres integration

    // save the strain
    noalias(mCurrentStrain) = StrainVector;

    // compute the strain increment
    Vector strainIncr = mLastStrain - mCurrentStrain;
    noalias(mLastStrainIncr) = strainIncr;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
    {
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "At StressIntegration, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << std::endl;
        std::cout << " mLastStrain: " << mLastStrain << std::endl;
        std::cout << " mCurrentStrain: " << mCurrentStrain << std::endl;
        std::cout << " strainIncr: " << strainIncr << std::endl;
        std::cout << " mLastStress: " << mLastStress << std::endl;
        std::cout << " mKm: " << mKm << std::endl;
        std::cout << " mGm: " << mGm << std::endl;
    }
    #endif

    // compute the elasticity matrix
    Matrix Ce(6, 6);
    this->CalculateElasticTangent(Ce, mKm, mGm);

    Vector stressTr = prod(Ce, strainIncr) + mLastStress;
    double pTr = (stressTr[0] + stressTr[1] + stressTr[2]) / 3;
    Vector devStressTr(6);
    devStressTr = this->getDeviatoricComp(devStressTr, stressTr);
    double J2 = this->getJ2( devStressTr );
    double qTr = sqrt(3.0 * J2);

    double yield_value = qTr * qTr / (mM * mM) + pTr * (pTr - mPc);
    mIsYielded = (yield_value / fabs(mPc)) > 1.0e-6;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
    {
        std::cout << " pTr: " << pTr << std::endl;
        std::cout << " qTr: " << qTr << std::endl;
        std::cout << " mM: " << mM << std::endl;
        std::cout << " mPc: " << mPc << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif
//    if(mIsYielded)
//    {
//        std::cout << "yielded at element " << mParentElementId << ", point " << mIntegrationPointIndex << std::endl;
//        std::cout << " pTr: " << pTr << std::endl;
//        std::cout << " qTr: " << qTr << std::endl;
//        std::cout << " mM: " << mM << std::endl;
//        std::cout << " mPc: " << mPc << std::endl;
//        std::cout << " yield_value: " << yield_value << std::endl;
//        std::cout << "---------------------------------------------" << std::endl;
//    }

//    if(pTr < 0.0)
    if(pTr < -1e-6)
    {
        std::cout << "WARNING: pTr is negative (" << pTr << ") at element " << mParentElementId << ", point " << mIntegrationPointIndex << std::endl;
//        mIsYielded = false;
        #ifdef DEBUG_CAM_CLAY
        if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
        {
            KRATOS_WATCH(strainIncr)
            KRATOS_WATCH(mLastStress)
        }
        #endif
    }

    if(mIsYielded)
    {
        int stat = this->returnMapping(pTr, qTr);

        if(stat >= 0)
        {
            double devNorm = sqrt(2.0 * J2);

            double factor;
            if( devNorm < TOL )
                factor = 0.0;
            else
                factor = mQk * sqrt(2.0 / 3) / devNorm;

            mCurrentStress[0] = factor * devStressTr[0] + mPk;
            mCurrentStress[1] = factor * devStressTr[1] + mPk;
            mCurrentStress[2] = factor * devStressTr[2] + mPk;
            mCurrentStress[3] = factor * devStressTr[3];
            mCurrentStress[4] = factor * devStressTr[4];
            mCurrentStress[5] = factor * devStressTr[5];
        }
        else
        {
            // if the return mapping is failed. The material switches back to linear elastic mode.
            mIsYielded = false;
            mPk = pTr;
            mPck = mPc;
            mDGamma = 0.0;
            mQk = qTr;
            noalias(mCurrentStress) = stressTr;
        }
    }
    else
    {
        mDGamma = 0.0;
        noalias(mCurrentStress) = stressTr;
    }

    // export the stress
    noalias(StressVector) = -1.0 * mCurrentStress;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
    {
        std::cout << " mCurrentStress: " << mCurrentStress << std::endl;
        std::cout << " mIsYielded: " << mIsYielded << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif

    #ifdef ENABLE_YIELD_PLOT
    mPlotUtil->RegisterPoint(pTr, qTr, mPk, mQk, true);
    #endif

    #ifdef LOG_CAM_CLAY
    if(mIsYielded)
    {
        std::cout << "Element " << mParentElementId
                  << " at integration point " << mIntegrationPointIndex
                  << " is yielded" << std::endl;
    }
    #endif
}

void CamClay3D::ComputeTangent(Matrix& AlgorithmicTangent)
{
    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
    {
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "At ComputeTangent, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << std::endl;
//        std::cout << " StrainVector: " << StrainVector << std::endl;
//        std::cout << " mCurrentStrain: " << mCurrentStrain << std::endl;
//        std::cout << " mLastStrain: " << mLastStrain << std::endl;
        std::cout << " mLastStrainIncr: " << mLastStrainIncr << std::endl;
        std::cout << " mIsYielded: " << mIsYielded << std::endl;
    }
    #endif

    // compute the elasticity matrix
    Matrix Ce(6, 6);
    this->CalculateElasticTangent(Ce, mKm, mGm);

    ///////////////////////////////////////////////////////////////////////
    // compute the tangent from the last converged stress
    if(mIsYielded)
    {
        Vector stressTr = prod(Ce, mLastStrainIncr) + mLastStress;
        Vector devStressTr(6);
        devStressTr = this->getDeviatoricComp(devStressTr, stressTr);
        double devNorm = sqrt(2.0 * this->getJ2(devStressTr));

        double devS, invDevNorm;
        if (devNorm < TOL)
        {
            devS = 1.0;
            invDevNorm = 0.0;
        }
        else
        {
            devS = sqrt(2.0 / 3) * mQk / devNorm;
            invDevNorm = 1.0 / devNorm;
        }

        Vector unitDev(6);
        noalias(unitDev) = devStressTr * invDevNorm;

        Vector fact;
        this->getFactors(fact, devS);
        #ifdef DEBUG_CAM_CLAY
        if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
        {
            std::cout << " stressTr: " << stressTr << std::endl;
            std::cout << " devStressTr: " << devStressTr << std::endl;
            std::cout << " devNorm: " << devNorm << std::endl;
            std::cout << " invDevNorm: " << invDevNorm << std::endl;
            std::cout << " unitDev: " << unitDev << std::endl;
            std::cout << " devS: " << devS << std::endl;
            std::cout << " fact: " << fact << std::endl;
        }
        #endif
        for(int i = 0; i < 6; ++i)
            for(int j = 0; j < 6; ++j)
                AlgorithmicTangent(i, j) = fact[0] * unit4thSym3D[i][j] +
                                           fact[1] * unit2nd3D[i] * unit2nd3D[j] +
                                           fact[2] * unit2nd3D[i] * unitDev[j] +
                                           fact[3] * unitDev[i] * unit2nd3D[j] +
                                           fact[4] * unitDev[i] * unitDev[j];
    }
    else
    {
        noalias(AlgorithmicTangent) = Ce;
    }

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
    {
        std::cout << " mKm: " << mKm << std::endl;
        std::cout << " mGm: " << mGm << std::endl;
//        std::cout << " Ce: " << Ce << std::endl;
        std::cout << " mCurrentStress: " << mCurrentStress << std::endl;
        std::cout << " AlgorithmicTangent: " << AlgorithmicTangent << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif
}

double CamClay3D::getP()
{
    double result = (mCurrentStress[0] + mCurrentStress[1] + mCurrentStress[2]) / 3;
//    if (result < 1.0)
//    {
//        result = 1.0;
//        KRATOS_THROW_ERROR(std::logic_error, "getP return value < 1.0", "")
//    }
    return result;
}

int CamClay3D::solveG(double& Pc, double pTr)
{
    // SPECIAL CASE 1
    if(pTr < 0.0)
    {
        Pc = mPc;
        return LOCAL_PC_SOLVE_NEGATIVE_PTRIAL;
    }

    // SPECIAL CASE 2
    if(fabs(mTheta*mDGamma) < TOL)
    {
        Pc = mPc;
        return LOCAL_PC_SOLVE_SMALL_INPUT;
    }

    /**************** METHOD 1 ****************/
    // USING NEWTON RAPHSON ITERATION. THIS METHOD IS QUADRATIC CONVERGENCE
    // BUT NOT ALWAYS CONVERGE TO THE SOLUTION
//    double result = mPc;
//    double dGp;

//    double expFact = exp( mTheta * mDGamma * ( 2 * pTr - result ) / ( 1 + 2 * mDGamma * mKm ) );
//    double Gp = mPc * expFact - result;

//    int it = 0;
//    const int max_it = 30;
//    while(fabs( Gp ) > TOL && it < max_it)
//    {
//        dGp = mPc * expFact * ( -mTheta * mDGamma / ( 1 + 2 * mDGamma * mKm ) ) - 1.0;
//        result -= Gp / dGp;

//        expFact = exp( mTheta * mDGamma * ( 2 * pTr - result ) / ( 1 + 2 * mDGamma * mKm ) );
//        Gp = mPc * expFact - result;

//        ++it;
//    }

//    if(it == max_it)
//    {
//        KRATOS_WATCH(expFact)
//        KRATOS_WATCH(mTheta)
//        KRATOS_WATCH(mDGamma)
//        KRATOS_WATCH(pTr)
//        KRATOS_WATCH(mKm)
//        KRATOS_WATCH(mPc)
//        KRATOS_WATCH(Gp)
//        KRATOS_WATCH(dGp)
//        KRATOS_WATCH(result)
//        KRATOS_WATCH(mParentElementId)
//        KRATOS_WATCH(mIntegrationPointIndex)
//        KRATOS_THROW_ERROR(std::runtime_error, "solveG does not converge in 30 steps", "")
//    }

    /**************** METHOD 2 ****************/
    // USING FIXED POINT ITERATION. THIS METHOD IS NOT ALWAYS CONVERGED.
//    double result = mPc, d_result = 1.0;
//    double old_result;
//    double coef = (1.0 + 2.0*mDGamma*mKm) / (mTheta*mDGamma);
//    const int max_it = 300;
//    int it = 0;
//    while(fabs(d_result/mPc) > TOL && it < max_it)
//    {
//        old_result = result;
//        result = 2.0*pTr - coef*log(result/mPc);
//        d_result = result - old_result;
//    }

//    if(it >= max_it && fabs(d_result/mPc) > TOL)
//    {
//        KRATOS_WATCH(mPc)
//        KRATOS_WATCH(coef)
//        KRATOS_WATCH(pTr)
//        KRATOS_WATCH(result)
//        std::stringstream ss;
//        ss << "solveG does not converge in 300 steps";
//        KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
//    }

//    if(mParentElementId == 2232 && mIntegrationPointIndex == 15)
//    {
//        KRATOS_WATCH(mPc)
//        KRATOS_WATCH(coef)
//        KRATOS_WATCH(pTr)
//        KRATOS_WATCH(mTheta)
//        KRATOS_WATCH(mDGamma)
//        KRATOS_WATCH(mKm)
//        std::cout << "solveG: " << result << std::endl;
//    }

    /**************** METHOD 3 ****************/
    double coef = (1.0 + 2.0*mDGamma*mKm) / (mTheta*mDGamma);
    double result;
    if(coef > 0.0)
    {
        // USING BISECTIVE ALGORITHM. THIS METHOD WILL ALWAYS CONVERGE
        // IF coef, mPc, pTr is positive
        double xleft = 0.5*std::min(mPc, 2.0*pTr);
        double xright = 2.0*std::max(mPc, 2.0*pTr);
        double res = mPc;
        const double tol = 1.0e-10;
        const int maxits = 300;
        int it = 0;
        while((fabs(res/mPc) > tol) && (it < maxits))
        {
            result = 0.5 * (xleft + xright);
            res = result + coef*log(result/mPc) - 2.0*pTr;
            if(res > 0)
                xright = result;
            else
                xleft = result;
            ++it;
        }
        if(it == maxits)
        {
//            KRATOS_WATCH(mPc)
//            KRATOS_WATCH(coef)
//            KRATOS_WATCH(pTr)
//            KRATOS_WATCH(mTheta)
//            KRATOS_WATCH(mDGamma)
//            KRATOS_WATCH(mKm)
//            KRATOS_WATCH(result)
//            std::stringstream ss;
//            ss << "solveG does not converge in 300 steps";
//            KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
            return LOCAL_PC_SOLVE_NOT_CONVERGED;
        }

        Pc = result;

        #ifdef DEBUG_CAM_CLAY
        if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
        {
            KRATOS_WATCH(Pc)
            KRATOS_WATCH(mPc)
            KRATOS_WATCH(coef)
            KRATOS_WATCH(pTr)
            KRATOS_WATCH(mTheta)
            KRATOS_WATCH(mDGamma)
            KRATOS_WATCH(mKm)
            std::cout << "solveG: " << result << std::endl;
        }
        #endif

        return 0;
    }
    else
    {
//        KRATOS_WATCH(mPc)
//        KRATOS_WATCH(coef)
//        KRATOS_WATCH(pTr)
//        KRATOS_WATCH(mTheta)
//        KRATOS_WATCH(mDGamma)
//        KRATOS_WATCH(mKm)
//        std::stringstream ss;
//        ss << "the coefficient is negative";
//        KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
        Pc = mPc;

        #ifdef DEBUG_CAM_CLAY
        if(mParentElementId == DEBUG_ELEMENT_ID && mIntegrationPointIndex == DEBUG_POINT_ID)
        {
            KRATOS_WATCH(Pc)
            KRATOS_WATCH(mPc)
            KRATOS_WATCH(coef)
            KRATOS_WATCH(pTr)
            KRATOS_WATCH(mTheta)
            KRATOS_WATCH(mDGamma)
            KRATOS_WATCH(mKm)
            std::cout << "solveG (invalid input): " << result << std::endl;
        }
        #endif

        return LOCAL_PC_SOLVE_INVALID_INPUT;
    }

    /******************************************/

    return 0;
}

void CamClay3D::getFactors(Vector& rResults, double devS)
{
    double a = 1.0 + 2.0 * mKm * mDGamma + mPck * mTheta * mDGamma;
    double a1 = ( 1.0 + mPck * mTheta * mDGamma ) / a;
    double a2 = -( 2.0 * mPk - mPck ) / a;
    double a3 = 2.0 * mPck * mTheta * mDGamma / a;
    double a4 = mTheta * mPck * (2.0 * mPk - mPck) / (a * mKm);
    double denom = 1.0 + 6.0 * mGm * mDGamma / (mM * mM);
    double a5 = sqrt(3.0 / 2) / denom;
    double a6 = -3.0 * mQk / ( denom * mM * mM );

    double b = -2.0 * mGm * 2.0 * mQk * a6 / (mM * mM) - mKm * ((2.0 * a2 - a4) * mPk - a2 * mPck );
    double b1 = -mKm * ( (a3 - 2.0 * a1) * mPk + a1 * mPck ) / b;
    double b2 = 2.0 * mGm * 2.0 * mQk * a5 / (b * mM * mM);

    if(rResults.size() != 5)
        rResults.resize(5);

    rResults[0] = 2.0 * mGm * devS;
    rResults[1] = mKm * (a1 + a2 * b1) - 2.0 * mGm * devS / 3;
    rResults[2] = mKm * a2 * b2;
    rResults[3] = 2.0 * mGm * sqrt(2.0 / 3) * (a6 * b1);
    rResults[4] = 2.0 * mGm * (sqrt(2.0 / 3) * (a5 + a6 * b2) - devS);
}

int CamClay3D::calculateF(std::vector<double>& rResults, double pTr, double qTr)
{
    mQk = qTr / ( 1 + 6 * mGm * mDGamma / (mM * mM) );
    int stat = this->solveG(mPck, pTr);
    if(stat < 0)
        return stat;
    mPk = (pTr + mDGamma * mKm * mPck) / (1 + 2 * mDGamma * mKm);

    double dFdPk  = 2.0 * mPk - mPck;
    double dFdQk  = 2.0 * mQk / (mM * mM);
    double dFdPck = -mPk;

    double dpdDGa  = -mKm * (2.0 * mPk - mPck) / ( 1.0 + (2.0 * mKm + mTheta * mPck) * mDGamma );
    double dqdDGa  = -mQk / (mDGamma + mM * mM / (6.0 * mGm));
    double dpcdDGa = mTheta * mPck * (2.0 * mPk - mPck) / ( 1.0 + (2.0 * mKm + mTheta * mPck) * mDGamma );

    if(rResults.size() != 2)
        rResults.resize(2);

    rResults[0] = mQk * mQk / (mM * mM) + mPk * (mPk - mPck);
    rResults[1] = dFdPk * dpdDGa + dFdQk * dqdDGa + dFdPck * dpcdDGa;

    return 0;
}

int CamClay3D::returnMapping(double pTr, double qTr)
{
    mDGamma = 0.0;
    std::vector<double> F(2);

    double denom;
    if(pow(fabs(pTr) + fabs(qTr), 2) > TOL)
        denom = pow(fabs(pTr) + fabs(qTr), 2);
    else
        denom = 1.0;

    int stat = this->calculateF(F, pTr, qTr);
    if(stat < 0)
        return stat;
    int it = 0;
    const int max_it = 30;
    while( fabs( F[0] / denom ) > TOL && it < max_it )
    {
        mDGamma -= F[0] / F[1];
        stat = this->calculateF(F, pTr, qTr);
        if(stat < 0)
            return stat;
        ++it;
    }
    if(it == max_it)
    {
//        KRATOS_WATCH(F[0])
//        KRATOS_WATCH(F[1])
//        KRATOS_WATCH(pTr)
//        KRATOS_WATCH(qTr)
//        KRATOS_THROW_ERROR(std::runtime_error, "returnMapping does not converge in 30 steps", "")
        return RETURN_MAPPING_NOT_CONVERGED;
    }

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == 1 && mIntegrationPointIndex == 0)
    {
        std::cout << "  After return mapping:" << std::endl;
        std::cout << "   mDGamma = " << mDGamma << std::endl;
        std::cout << "   mPk = " << mPk << std::endl;
        std::cout << "   mQk = " << mQk << std::endl;
        std::cout << "   mPck = " << mPck << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif

    return 0;
}

void CamClay3D::CalculateElasticTangent(Matrix& C, double K, double G)
{
    double C1 = 2.0 * G;
    double C2 = K - 2.0 * G / 3;

    if(C.size1() != 6 || C.size2() != 6)
        C.resize(6, 6);

    for(int i = 0; i < 6; ++i)
        for(int j = i; j < 6; ++j)
            C(i, j) = C1 * unit4thSym3D[i][j] + C2 * unit2nd3D[i] * unit2nd3D[j];//p.622 De = 2G Is + (K - 2/3 G) I x I

    for(int j = 0; j < 5; ++j)
        for(int i = j+1; i < 6; ++i)
            C(i, j) = C(j, i);
}

Vector& CamClay3D::getDeviatoricComp(Vector& devStress, const Vector& stress)
{
    double p = ( stress[0] + stress[1] + stress[2] ) / 3;
    devStress[0] = stress[0] - p;
    devStress[1] = stress[1] - p;
    devStress[2] = stress[2] - p;
    devStress[3] = stress[3];
    devStress[4] = stress[4];
    devStress[5] = stress[5];

    return devStress;
}

double CamClay3D::getJ2(const Vector& s)
{
    return 0.5 * ( s[0] * s[0] + s[1] * s[1] + s[2] * s[2] ) + s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
}

} // Namespace Kratos

#undef DEBUG_CAM_CLAY
#undef DEBUG_CAM_CLAY_PRECONSOLIDATION
#undef DEBUG_ELEMENT_ID
#undef DEBUG_POINT_ID
#undef LOG_CAM_CLAY
#undef ENABLE_YIELD_PLOT
#undef RETURN_MAPPING_NOT_CONVERGED
#undef LOCAL_PC_SOLVE_INVALID_INPUT
#undef LOCAL_PC_SOLVE_NEGATIVE_PTRIAL
#undef LOCAL_PC_SOLVE_SMALL_INPUT
#undef LOCAL_PC_SOLVE_NOT_CONVERGED

