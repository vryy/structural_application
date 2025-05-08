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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2009-01-14 17:14:12 $
*   Revision:            $Revision: 1.13 $
*
* ***********************************************************/

// System includes
#include <iostream>

// External includes
#include <cmath>

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws/isotropic_3d.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"

namespace Kratos
{

template<class TNodeType>
Isotropic3DImpl<TNodeType>::Isotropic3DImpl()
    : BaseType()
{
}

template<class TNodeType>
Isotropic3DImpl<TNodeType>::~Isotropic3DImpl()
{
}

template<class TNodeType>
bool Isotropic3DImpl<TNodeType>::Has( const Variable<int>& rThisVariable )
{
    return false;
}

template<class TNodeType>
bool Isotropic3DImpl<TNodeType>::Has( const Variable<DataType>& rThisVariable )
{
    if ( rThisVariable == PRESTRESS_FACTOR || rThisVariable == YOUNG_MODULUS || rThisVariable == POISSON_RATIO )
        return true;

    return false;
}

template<class TNodeType>
bool Isotropic3DImpl<TNodeType>::Has( const Variable<VectorType>& rThisVariable )
{
    if ( rThisVariable == INSITU_STRESS )
        return true;

    if ( rThisVariable == PRESTRESS )
        return true;

    if ( rThisVariable == STRESSES )
        return true;

    return false;
}

template<class TNodeType>
bool Isotropic3DImpl<TNodeType>::Has( const Variable<MatrixType>& rThisVariable )
{
    if ( rThisVariable == CAUCHY_STRESS_TENSOR )
        return true;
    if(rThisVariable == ALGORITHMIC_TANGENT || rThisVariable == ELASTIC_TANGENT || rThisVariable == THREED_ALGORITHMIC_TANGENT)
        return true;
    return false;
}

template<class TNodeType>
int& Isotropic3DImpl<TNodeType>::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    if (rThisVariable == IS_SHAPE_FUNCTION_REQUIRED)
        rValue = 0;

    return rValue;
}

template<class TNodeType>
typename Isotropic3DImpl<TNodeType>::DataType& Isotropic3DImpl<TNodeType>::GetValue( const Variable<DataType>& rThisVariable, DataType& rValue )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
    {
        rValue = mPrestressFactor;
        return rValue;
    }

    if(rThisVariable == YOUNG_MODULUS )
    {
       rValue = mE;
       return rValue;
    }

    if ( rThisVariable == POISSON_RATIO )
    {
        rValue = mNU;
        return rValue;
    }

    if(rThisVariable == DAMAGE)
    {
        rValue = 0.00;
        return rValue;
    }

    if (rThisVariable == DELTA_TIME)
    {
        rValue = sqrt(mE/mDE);
        return rValue;
    }

    if (rThisVariable == PRESSURE_P)
    {
        rValue = -(mCurrentStress[0] + mCurrentStress[1] + mCurrentStress[2]) / 3;
        return rValue;
    }

    if (rThisVariable == PRESSURE_Q || rThisVariable == VON_MISES_STRESS)
    {
        DataType p = (mCurrentStress[0] + mCurrentStress[1] + mCurrentStress[2]) / 3;
        DataType sxx = mCurrentStress[0] - p;
        DataType syy = mCurrentStress[1] - p;
        DataType szz = mCurrentStress[2] - p;
        DataType sxy = mCurrentStress[3];
        DataType syz = mCurrentStress[4];
        DataType sxz = mCurrentStress[5];

        rValue = std::sqrt( 3.0 * ( 0.5*(sxx*sxx + syy*syy + szz*szz) + sxy*sxy + syz*syz + sxz*sxz ) );
        return rValue;
    }

    rValue = 0.00;
    return rValue;
}

template<class TNodeType>
typename Isotropic3DImpl<TNodeType>::VectorType& Isotropic3DImpl<TNodeType>::GetValue( const Variable<VectorType>& rThisVariable, VectorType& rValue )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        const unsigned int size = mPrestress.size();
        if(rValue.size() != size)
            rValue.resize(size, false );
        noalias(rValue) = mPrestressFactor * mPrestress;
        return rValue;
    }

    if ( rThisVariable == STRESSES || rThisVariable == THREED_STRESSES )
    {
        const unsigned int size = mCurrentStress.size();
        rValue.resize(size, false );
        noalias(rValue)  = mCurrentStress;
        return rValue;
    }

    if ( rThisVariable == THREED_STRAIN )
    {
        // REF: http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_stress.cfm
        const unsigned int size = mCurrentStress.size();
        rValue.resize(size, false );
        VectorType Stress = mCurrentStress + mPrestressFactor * mPrestress;
        rValue(0) = (Stress(0) - mNU*Stress(1) - mNU*Stress(2)) / mE;
        rValue(1) = (Stress(1) - mNU*Stress(0) - mNU*Stress(2)) / mE;
        rValue(2) = (Stress(2) - mNU*Stress(0) - mNU*Stress(1)) / mE;
        rValue(3) = 2.0*(1.0+mNU)*Stress(3) / mE;
        rValue(4) = 2.0*(1.0+mNU)*Stress(4) / mE;
        rValue(5) = 2.0*(1.0+mNU)*Stress(5) / mE;
        return rValue;
    }

    if ( rThisVariable == PLASTIC_STRAIN_VECTOR )
    {
        rValue = ZeroVectorType( 6 );
        return( rValue );
    }

    if ( rThisVariable == INTERNAL_VARIABLES )
    {
        rValue = ZeroVectorType( 1 );
        return( rValue );
    }

    return rValue;
}

template<class TNodeType>
typename Isotropic3DImpl<TNodeType>::MatrixType& Isotropic3DImpl<TNodeType>::GetValue( const Variable<MatrixType>& rThisVariable, MatrixType& rValue )
{
    if(rThisVariable == ALGORITHMIC_TANGENT || rThisVariable == ELASTIC_TANGENT || rThisVariable == THREED_ALGORITHMIC_TANGENT)
    {
        if(rValue.size1() != 6 || rValue.size2() != 6)
            rValue.resize(6, 6, false);
        CalculateElasticMatrix( rValue, mE, mNU );
    }
    else if(rThisVariable == ELASTIC_STRAIN_TENSOR)
    {
        VectorType StrainVector(6);
        this->CalculateStrain(mE, mNU, mCurrentStress + mPrestressFactor * mPrestress, StrainVector);
        SD_MathUtils<DataType>::StrainVectorToTensor(StrainVector, rValue);
    }
    else if(rThisVariable == CAUCHY_STRESS_TENSOR)
    {
        if (rValue.size1() != 3 || rValue.size2() != 3)
            rValue.resize(3, 3, false);
        SD_MathUtils<DataType>::StrainVectorToTensor(mCurrentStress, rValue);
    }

    return rValue ;
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::SetValue( const Variable<DataType>& rThisVariable, const DataType& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
    if ( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::SetValue( const Variable<array_1d<DataType, 3> >& rThisVariable,
                            const array_1d<DataType, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::SetValue( const Variable<VectorType>& rThisVariable, const VectorType& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        noalias(mPrestress) = rValue;
    }
    else if ( rThisVariable == STRESSES || rThisVariable == INITIAL_STRESS )
    {
        noalias(mCurrentStress) = rValue;
    }
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::SetValue( const Variable<MatrixType>& rThisVariable, const MatrixType& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mCurrentStress = ZeroVectorType( 6 );
    mPrestress = ZeroVectorType( 6 );
    mPrestressFactor = 1.0;
    mE = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mDE = props[DENSITY];
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    noalias(mCurrentStress) = mPrestressFactor*mPrestress;
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::CalculateMaterialResponseCauchy (typename BaseType::Parameters& rValues)
{
    if (rValues.IsSetConstitutiveMatrix())
    {
        MatrixType& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        if(AlgorithmicTangent.size1() != 6 || AlgorithmicTangent.size2() != 6)
            AlgorithmicTangent.resize(6, 6, false);
        CalculateElasticMatrix( AlgorithmicTangent, mE, mNU );
    }

    if (rValues.IsSetStressVector())
    {
        const VectorType& StrainVector = rValues.GetStrainVector();
        VectorType& StressVector = rValues.GetStressVector();
        if(StressVector.size() != 6)
            StressVector.resize(6, false);

        if (rValues.IsSetConstitutiveMatrix())
        {
            MatrixType& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
            CalculateStress( StrainVector, AlgorithmicTangent, StressVector );
        }
        else
        {
            MatrixType AlgorithmicTangent(6, 6);
            CalculateElasticMatrix( AlgorithmicTangent, mE, mNU );
            CalculateStress( StrainVector, AlgorithmicTangent, StressVector );
        }
    }
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::CalculateMaterialResponse( const VectorType& StrainVector,
        const MatrixType& DeformationGradient,
        VectorType& StressVector,
        MatrixType& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo,
        const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        bool CalculateStresses,
        int CalculateTangent,
        bool SaveInternalVariables )
{
    typename BaseType::Parameters const_params;
    VectorType ThisStrainVector = StrainVector;
    const_params.SetStrainVector(ThisStrainVector);
    if (CalculateStresses)
        const_params.SetStressVector(StressVector);
    if (CalculateTangent)
        const_params.SetConstitutiveMatrix(AlgorithmicTangent);

    this->CalculateMaterialResponseCauchy(const_params);
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::CalculateElasticMatrix( MatrixType& C, const DataType E, const DataType NU )
{
    //setting up material matrix
    DataType c1 = E / (( 1.00 + NU ) * ( 1.00 - 2.00 * NU ) );
    DataType c2 = c1 * ( 1.00 - NU );
    DataType c3 = c1 * NU;
    DataType c4 = c1 * 0.5 * ( 1.00 - 2.00 * NU );

    //filling material matrix
    noalias(C) = ZeroMatrix(6, 6);
    C( 0, 0 ) = c2;
    C( 0, 1 ) = c3;
    C( 0, 2 ) = c3;
    C( 1, 0 ) = c3;
    C( 1, 1 ) = c2;
    C( 1, 2 ) = c3;
    C( 2, 0 ) = c3;
    C( 2, 1 ) = c3;
    C( 2, 2 ) = c2;
    C( 3, 3 ) = c4;
    C( 4, 4 ) = c4;
    C( 5, 5 ) = c4;
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::CalculateStress( const VectorType& StrainVector, const MatrixType& AlgorithmicTangent, VectorType& StressVector )
{
    if ( StressVector.size() != 6 )
    {
        StressVector.resize( 6, false );
    }

    noalias( StressVector ) = prod( AlgorithmicTangent, StrainVector ) - mPrestressFactor * mPrestress;

    noalias(mCurrentStress) = StressVector;
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::CalculateStress( const DataType E, const DataType NU, const VectorType& StrainVector, VectorType& StressVector ) const
{
    if ( StressVector.size() != 6 )
    {
        StressVector.resize( 6, false );
    }

    MatrixType Ce( 6, 6 );

    this->CalculateElasticMatrix( Ce, E, NU );

    noalias( StressVector ) = prod( Ce, StrainVector );
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::CalculateStrain( const DataType E, const DataType NU, const VectorType& StressVector, VectorType& StrainVector ) const
{
    const DataType c1 = 1.0 / E;
    const DataType c2 = -NU * c1;
    const DataType c3 = 2 * (1.0 + NU) / E;

    StrainVector(0) = c1 * StressVector(0) + c2 * StressVector(1) + c2 * StressVector(2);
    StrainVector(1) = c2 * StressVector(0) + c1 * StressVector(1) + c2 * StressVector(2);
    StrainVector(2) = c2 * StressVector(0) + c2 * StressVector(1) + c1 * StressVector(2);
    StrainVector(3) = c3 * StressVector(3);
    StrainVector(4) = c3 * StressVector(4);
    StrainVector(5) = c3 * StressVector(5);
}

template<class TNodeType>
void Isotropic3DImpl<TNodeType>::CalculateCauchyStresses(
    VectorType& rCauchy_StressVector,
    const MatrixType& rF,
    const VectorType& rPK2_StressVector,
    const VectorType& rGreenLagrangeStrainVector )
{
    MatrixType S = MathUtils<DataType>::StressVectorToTensor( rPK2_StressVector );

    DataType J = MathUtils<DataType>::Det3( rF );
    boost::numeric::ublas::bounded_matrix<DataType, 3, 3> mstemp;
    boost::numeric::ublas::bounded_matrix<DataType, 3, 3> msaux;

    noalias( mstemp ) = prod( rF, S );
    noalias( msaux ) = prod( mstemp, trans( rF ) );
    msaux *= J;

    if ( rCauchy_StressVector.size() != 6 )
        rCauchy_StressVector.resize( 6 );

    rCauchy_StressVector[0] = msaux( 0, 0 );

    rCauchy_StressVector[1] = msaux( 1, 1 );

    rCauchy_StressVector[2] = msaux( 2, 2 );

    rCauchy_StressVector[3] = msaux( 0, 1 );

    rCauchy_StressVector[4] = msaux( 0, 2 );

    rCauchy_StressVector[5] = msaux( 1, 2 );
}

template<class TNodeType>
int Isotropic3DImpl<TNodeType>::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY

    if ( !props.Has( YOUNG_MODULUS ) || !props.Has(POISSON_RATIO) )
    {
        KRATOS_ERROR << "this constitutive law requires YOUNG_MODULUS and POISSON_RATIO given as KRATOS variables";
    }

    double nu = props[POISSON_RATIO];

    if ( nu > 0.499 && nu < 0.501 )
    {
        KRATOS_ERROR << "invalid poisson ratio in input, close to incompressibility";
        return -1;
    }

    return 0;

    KRATOS_CATCH( "" );
}

template class Isotropic3DImpl<RealNode>;
template class Isotropic3DImpl<ComplexNode>;
template class Isotropic3DImpl<GComplexNode>;

} // Namespace Kratos
