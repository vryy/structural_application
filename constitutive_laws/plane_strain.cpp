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
*   Last Modified by:    $Author: kazem $
*   Date:                $Date: 2009-01-16 10:50:04 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws/plane_strain.h"
#include "constitutive_laws/isotropic_3d.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"

namespace Kratos
{

template<class TNodeType>
PlaneStrainImpl<TNodeType>::PlaneStrainImpl()
    : BaseType()
{
}

template<class TNodeType>
PlaneStrainImpl<TNodeType>::~PlaneStrainImpl()
{
}

template<class TNodeType>
bool PlaneStrainImpl<TNodeType>::Has(const Variable<int>& rThisVariable)
{
    return false;
}

template<class TNodeType>
bool PlaneStrainImpl<TNodeType>::Has( const Variable<DataType>& rThisVariable )
{
    return false;
}

template<class TNodeType>
bool PlaneStrainImpl<TNodeType>::Has( const Variable<VectorType>& rThisVariable )
{
    if(rThisVariable == VARSEL(DataType, STRESSES))
        return true;

    return false;
}

template<class TNodeType>
bool PlaneStrainImpl<TNodeType>::Has( const Variable<MatrixType>& rThisVariable )
{
    if(rThisVariable == VARSEL(DataType, ALGORITHMIC_TANGENT)
    || rThisVariable == VARSEL(DataType, THREED_ALGORITHMIC_TANGENT))
        return true;

    return false;
}

template<class TNodeType>
int& PlaneStrainImpl<TNodeType>::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    if (rThisVariable == IS_SHAPE_FUNCTION_REQUIRED)
        rValue = 0;

    return rValue;
}

template<class TNodeType>
typename PlaneStrainImpl<TNodeType>::DataType& PlaneStrainImpl<TNodeType>::GetValue( const Variable<DataType>& rThisVariable, DataType& rValue )
{
    if(rThisVariable == VARSEL(DataType, DELTA_TIME))
        rValue = sqrt(mE/mDE);

    if(rThisVariable == VARSEL(DataType, PRESTRESS_FACTOR))
        rValue = mPrestressFactor;

    if(rThisVariable == VARSEL(DataType, YOUNG_MODULUS))
        rValue = mE;

    if(rThisVariable == VARSEL(DataType, POISSON_RATIO))
        rValue = mNU;

    if (rThisVariable == VARSEL(DataType, PRESSURE_P))
    {
        DataType o_zz = mNU * (mCurrentStress[0] + mCurrentStress[1]);
        rValue = -(mCurrentStress[0] + mCurrentStress[1] + o_zz) / 3;
        return rValue;
    }

    if (rThisVariable == VARSEL(DataType, PRESSURE_Q) || rThisVariable == VARSEL(DataType, VON_MISES_STRESS))
    {
        DataType o_zz = mNU * (mCurrentStress[0] + mCurrentStress[1]);
        DataType p = (mCurrentStress[0] + mCurrentStress[1] + o_zz) / 3;
        DataType sxx = mCurrentStress[0] - p;
        DataType syy = mCurrentStress[1] - p;
        DataType szz = o_zz - p;
        DataType sxy = mCurrentStress[2];
        DataType syz = 0.0;
        DataType sxz = 0.0;

        rValue = std::sqrt( 3.0 * ( 0.5*(sxx*sxx + syy*syy + szz*szz) + sxy*sxy + syz*syz + sxz*sxz ) );
        return rValue;
    }

    return rValue;
}

template<class TNodeType>
typename PlaneStrainImpl<TNodeType>::VectorType& PlaneStrainImpl<TNodeType>::GetValue( const Variable<VectorType>& rThisVariable, VectorType& rValue )
{
    if constexpr (std::is_arithmetic<DataType>::value)
    {
        if ( rThisVariable == STRESSES || rThisVariable == STRESSES_OLD )
        {
            rValue = mCurrentStress;
        }
        else if ( rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS )
        {
            rValue = mPreStress;
        }
        else if ( rThisVariable == THREED_STRESSES )
        {
            if(rValue.size() != 6)
                rValue.resize(6, false);
            rValue(0) = mCurrentStress(0);
            rValue(1) = mCurrentStress(1);
            rValue(2) = mNU * (mCurrentStress(0) + mCurrentStress(1));
            rValue(3) = mCurrentStress(2);
            rValue(4) = 0.0;
            rValue(5) = 0.0;
        }
        else if ( rThisVariable == THREED_PRESTRESS )
        {
            if(rValue.size() != 6)
                rValue.resize(6, false);
            rValue(0) = mPreStress(0);
            rValue(1) = mPreStress(1);
            rValue(2) = mNU * (mPreStress(0) + mPreStress(1));
            rValue(3) = mPreStress(2);
            rValue(4) = 0.0;
            rValue(5) = 0.0;
        }
        else if ( rThisVariable == THREED_STRAIN )
        {
            // REF: http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_strain.cfm
            if(rValue.size() != 6)
                rValue.resize(6, false);
            DataType aux = (1.0+mNU)/mE;
            rValue(0) = aux * ((1.0-mNU)*mCurrentStress(0) - mNU*mCurrentStress(1));
            rValue(1) = aux * ((1.0-mNU)*mCurrentStress(1) - mNU*mCurrentStress(0));
            rValue(2) = 0.0;
            rValue(3) = 2.0*aux*mCurrentStress(2);
            rValue(4) = 0.0;
            rValue(5) = 0.0;
        }
    }

    return rValue;
}

template<class TNodeType>
typename PlaneStrainImpl<TNodeType>::MatrixType& PlaneStrainImpl<TNodeType>::GetValue( const Variable<MatrixType>& rThisVariable, MatrixType& rValue )
{
    if constexpr (std::is_arithmetic<DataType>::value)
    {
        if(rThisVariable == ALGORITHMIC_TANGENT || rThisVariable == ELASTIC_TANGENT)
        {
            if(rValue.size1() != 3 || rValue.size2() != 3)
                rValue.resize(3, 3, false);
            CalculateElasticMatrix( rValue, mE, mNU );
        }
        else if( rThisVariable == THREED_ALGORITHMIC_TANGENT )
        {
            if (rValue.size1() != 6 || rValue.size2() != 6)
                rValue.resize(6, 6, false);
            Isotropic3D::CalculateElasticMatrix(rValue, mE, mNU);
        }
        else if(rThisVariable == ELASTIC_STRAIN_TENSOR)
        {
            VectorType StrainVector(3);
            this->CalculateStrain(mE, mNU, mCurrentStress + mPrestressFactor * mPreStress, StrainVector);
            SD_MathUtils<DataType>::StrainVectorToTensor(StrainVector, rValue);
        }
        else if(rThisVariable == CAUCHY_STRESS_TENSOR)
        {
            if (rValue.size1() != 3 || rValue.size2() != 3)
                rValue.resize(3, 3, false);
            SD_MathUtils<DataType>::StrainVectorToTensor(mCurrentStress, rValue);
        }
    }

    return rValue ;
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::SetValue( const Variable<DataType>& rThisVariable, const DataType& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if constexpr (std::is_arithmetic<DataType>::value)
    {
        if ( rThisVariable == PRESTRESS_FACTOR )
            mPrestressFactor = rValue;
        if ( rThisVariable == YOUNG_MODULUS )
            mE = rValue;
        if ( rThisVariable == POISSON_RATIO )
            mNU = rValue;
    }
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::SetValue( const Variable<VectorType>& rThisVariable, const VectorType& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if constexpr (std::is_arithmetic<DataType>::value)
    {
        if ( rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS )
        {
            if (rValue.size() == 3)
            {
                noalias(mPreStress) = rValue;
            }
            else if (rValue.size() == 4 || rValue.size() == 6)
            {
                mPreStress(0) = rValue(0);
                mPreStress(1) = rValue(1);
                mPreStress(2) = rValue(3);
            }
        }
        else if ( rThisVariable == STRESSES || rThisVariable == INITIAL_STRESS )
        {
            if(mCurrentStress.size() != rValue.size())
                mCurrentStress.resize(rValue.size(), false);
            noalias(mCurrentStress) = rValue;
        }
        else if ( rThisVariable == THREED_STRESSES )
        {
            if (rValue.size() == 3)
            {
                noalias(mCurrentStress) = rValue;
            }
            else if (rValue.size() == 4 || rValue.size() == 6)
            {
                mCurrentStress(0) = rValue(0);
                mCurrentStress(1) = rValue(1);
                mCurrentStress(2) = rValue(3);
            }
        }
    }
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::SetValue( const Variable<MatrixType>& rThisVariable, const MatrixType& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::CalculateMaterialResponseCauchy (typename BaseType::Parameters& rValues)
{
    if (rValues.IsSetStressVector())
    {
        const VectorType& StrainVector = rValues.GetStrainVector();
        VectorType& StressVector = rValues.GetStressVector();

        if(StressVector.size() != 3)
            StressVector.resize(3, false);
        CalculateStress( StrainVector, StressVector );
    }

    if (rValues.IsSetConstitutiveMatrix())
    {
        MatrixType& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        if(AlgorithmicTangent.size1() != 3 || AlgorithmicTangent.size2() != 3)
            AlgorithmicTangent.resize(3, 3, false);
        CalculateConstitutiveMatrix( AlgorithmicTangent );
    }
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::CalculateMaterialResponse( const VectorType& StrainVector,
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
    const_params.SetStressVector(StressVector);
    const_params.SetConstitutiveMatrix(AlgorithmicTangent);

    this->CalculateMaterialResponseCauchy(const_params);
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mCurrentStress = ZeroVectorType( 3 );
    mPreStress = ZeroVectorType( 3 );
    mPrestressFactor = 1.0;
    mE  = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mDE = props[DENSITY];
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    noalias(mCurrentStress) = ZeroVectorType(3);
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::InitializeNonLinearIteration( const Properties& rMaterialProperties,
                                                const GeometryType& rElementGeometry,
                                                const Vector& rShapeFunctionsValues,
                                                const ProcessInfo& rCurrentProcessInfo )
{
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::FinalizeNonLinearIteration( const Properties& rMaterialProperties,
                                              const GeometryType& rElementGeometry,
                                              const Vector& rShapeFunctionsValues,
                                              const ProcessInfo& rCurrentProcessInfo )
{
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::CalculateElasticMatrix( MatrixType& C, const DataType E, const DataType NU ) const
{
    DataType c1 = E * ( 1.00 - NU ) / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
    DataType c2 = E * NU / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
    DataType c3 = 0.5 * E / ( 1.00 + NU );

    C( 0, 0 ) = c1;
    C( 0, 1 ) = c2;
    C( 0, 2 ) = 0.0;
    C( 1, 0 ) = c2;
    C( 1, 1 ) = c1;
    C( 1, 2 ) = 0.0;
    C( 2, 0 ) = 0.0;
    C( 2, 1 ) = 0.0;
    C( 2, 2 ) = c3;
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::CalculateStress( const VectorType& StrainVector, VectorType& StressVector )
{
    this->CalculateStress(mE, mNU, StrainVector, StressVector);

    noalias(StressVector) -= mPrestressFactor*mPreStress;

    noalias( mCurrentStress ) = StressVector;
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::CalculateStress( const DataType E, const DataType NU, const VectorType& StrainVector, VectorType& StressVector ) const
{
    DataType c1 = E * ( 1.00 - NU ) / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
    DataType c2 = E * NU / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
    DataType c3 = 0.5 * E / ( 1.00 + NU );

    // compute the stress based on strain
    StressVector[0] = c1 * StrainVector[0] + c2 * StrainVector[1];
    StressVector[1] = c1 * StrainVector[1] + c2 * StrainVector[0];
    StressVector[2] = c3 * StrainVector[2];
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::CalculateStrain(const DataType E, const DataType NU, const VectorType& StressVector, VectorType& StrainVector ) const
{
    const DataType c1 = ( 1.00 - NU * NU ) / E;
    const DataType c2 = -NU * ( 1.00 + NU ) / E;
    const DataType c3 = 2 * ( 1.00 + NU ) / E;

    StrainVector(0) = c1 * StressVector(0) + c2 * StressVector(1);
    StrainVector(1) = c2 * StressVector(0) + c1 * StressVector(1);
    StrainVector(2) = c3 * StressVector(2);
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::CalculateConstitutiveMatrix( MatrixType& rResult ) const
{
    CalculateElasticMatrix( rResult, mE, mNU );
}

template<class TNodeType>
std::size_t PlaneStrainImpl<TNodeType>::GetStrainSize() const
{
    return 3;
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::CalculateCauchyStresses(
    VectorType& rCauchy_StressVector,
    const MatrixType& rF,
    const VectorType& rPK2_StressVector,
    const VectorType& rGreenLagrangeStrainVector ) const
{
    MatrixType S = MathUtils<DataType>::StressVectorToTensor( rPK2_StressVector );

    DataType J = MathUtils<DataType>::Det2( rF );

    boost::numeric::ublas::bounded_matrix<DataType, 2, 2> temp;
    boost::numeric::ublas::bounded_matrix<DataType, 2, 2> aux;

    noalias( temp ) = prod( rF, S );
    noalias( aux ) = prod( temp, trans( rF ) );
    aux *= J;

    if ( rCauchy_StressVector.size() != 3 )
        rCauchy_StressVector.resize( 3 );

    rCauchy_StressVector[0] = aux( 0, 0 );

    rCauchy_StressVector[1] = aux( 1, 1 );

    rCauchy_StressVector[2] = aux( 0,1 );
}

template<class TNodeType>
void PlaneStrainImpl<TNodeType>::Calculate( const Variable<MatrixType>& rVariable, MatrixType& rResult,
                             const ProcessInfo& rCurrentProcessInfo ) const
{
}

template<class TNodeType>
int PlaneStrainImpl<TNodeType>::Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo) const
{
    if(YOUNG_MODULUS.Key() == 0 || props[YOUNG_MODULUS]<= 0.00)
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value";

    const double nu = std::abs(props[POISSON_RATIO]);
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );
    if(POISSON_RATIO.Key() == 0 || check==true) // props[POISSON_RATIO] == 1.00 || props[POISSON_RATIO] == -1.00)
        KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value";

    if(DENSITY.Key() == 0 || props[DENSITY]<0.00)
        KRATOS_ERROR << "DENSITY has Key zero or invalid value";

    return 0;
}

template class PlaneStrainImpl<RealNode>;
template class PlaneStrainImpl<ComplexNode>;
template class PlaneStrainImpl<GComplexNode>;

} // Namespace Kratos
