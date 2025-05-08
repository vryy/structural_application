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
/* *********************************************************
 *
 *   Modified by:         $Author: hbui $
 *   Date:                $Date: 2 May 2025 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

// System includes


// External includes

// Project includes
#include "custom_elements/kinematic_linear_anti_plane.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"

//#define ENABLE_DEBUG_CONSTITUTIVE_LAW

namespace Kratos
{
    template<class TNodeType>
    BaseKinematicLinearAntiPlane<TNodeType>::BaseKinematicLinearAntiPlane( IndexType NewId, typename GeometryType::Pointer pGeometry )
        : BaseType( NewId, pGeometry )
    {}

    template<class TNodeType>
    BaseKinematicLinearAntiPlane<TNodeType>::BaseKinematicLinearAntiPlane( IndexType NewId,
            typename GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties )
        : BaseType( NewId, pGeometry, pProperties )
    {}

    template<class TNodeType>
    typename BaseKinematicLinearAntiPlane<TNodeType>::ElementType::Pointer BaseKinematicLinearAntiPlane<TNodeType>::Create( IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties ) const
    {
        return typename ElementType::Pointer( new BaseKinematicLinearAntiPlane( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    template<class TNodeType>
    typename BaseKinematicLinearAntiPlane<TNodeType>::ElementType::Pointer BaseKinematicLinearAntiPlane<TNodeType>::Create( IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties ) const
    {
        return typename ElementType::Pointer( new BaseKinematicLinearAntiPlane( NewId, pGeom, pProperties ) );
    }

    template<class TNodeType>
    BaseKinematicLinearAntiPlane<TNodeType>::~BaseKinematicLinearAntiPlane()
    {
    }

    template<class TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        // compute the full 3D matrix and vector
        MatrixType FullLeftHandSideMatrix;
        VectorType FullRightHandSideVector;

        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        BaseType::CalculateAll( FullLeftHandSideMatrix, FullRightHandSideVector, rCurrentProcessInfo,
                CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

        // extract the matrix and vector for the z-component
        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int mat_size = this->GetGeometry().size();

        if ( rLeftHandSideMatrix.size1() != mat_size
          || rLeftHandSideMatrix.size2() != mat_size )
        {
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );
            noalias( rLeftHandSideMatrix ) = ZeroMatrixType( mat_size, mat_size ); //resetting LHS
        }

        if ( rRightHandSideVector.size() != mat_size )
        {
            rRightHandSideVector.resize( mat_size, false );
            noalias( rRightHandSideVector ) = ZeroVectorType( mat_size ); //resetting RHS
        }

        for (unsigned int i = 0; i < mat_size; ++i)
        {
            rRightHandSideVector(i) = FullRightHandSideVector(dim*i + 2);
            for (unsigned int j = 0; j < mat_size; ++j)
            {
                rLeftHandSideMatrix(i, j) = FullLeftHandSideMatrix(dim*i + 2, dim*j + 2);
            }
        }
    }

    template<class TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        // compute the full 3D vector
        MatrixType FullLeftHandSideMatrix;
        VectorType FullRightHandSideVector;

        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        BaseType::CalculateAll( FullLeftHandSideMatrix, FullRightHandSideVector, rCurrentProcessInfo,
                CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

        // extract the vector for the z-component
        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int mat_size = this->GetGeometry().size();

        if ( rRightHandSideVector.size() != mat_size )
        {
            rRightHandSideVector.resize( mat_size, false );
        }

        for (unsigned int i = 0; i < mat_size; ++i)
        {
            rRightHandSideVector(i) = FullRightHandSideVector(dim*i + 2);
        }
    }

    template<class TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        MatrixType FullMassMatrix;
        BaseType::CalculateMassMatrix( FullMassMatrix, rCurrentProcessInfo );

        // extract the matrix for the z-component
        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int mat_size = this->GetGeometry().size();

        if ( rMassMatrix.size1() != mat_size
          || rMassMatrix.size2() != mat_size )
        {
            rMassMatrix.resize( mat_size, mat_size, false );
            noalias( rMassMatrix ) = ZeroMatrixType( mat_size, mat_size ); //resetting LHS
        }

        for (unsigned int i = 0; i < mat_size; ++i)
        {
            for (unsigned int j = 0; j < mat_size; ++j)
            {
                rMassMatrix(i, j) = FullMassMatrix(dim*i + 2, dim*j + 2);
            }
        }
    }

    template<class TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        // extract the matrix for the z-component
        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int mat_size = this->GetGeometry().size();

        if ( rDampMatrix.size1() != mat_size
          || rDampMatrix.size2() != mat_size )
        {
            rDampMatrix.resize( mat_size, mat_size, false );
            noalias( rDampMatrix ) = ZeroMatrixType( mat_size, mat_size ); //resetting LHS
        }

        // Rayleigh damping

        ValueType alpha = 0.0, beta = 0.0;

        if(this->GetProperties().Has(RAYLEIGH_DAMPING_ALPHA))
        {
            alpha = this->GetProperties()[RAYLEIGH_DAMPING_ALPHA];
        }

        if(this->GetProperties().Has(RAYLEIGH_DAMPING_BETA))
        {
            beta = this->GetProperties()[RAYLEIGH_DAMPING_BETA];
        }

        if (alpha > 0.0)
        {
            this->CalculateMassMatrix( rDampMatrix, rCurrentProcessInfo );

            rDampMatrix *= alpha;
        }

        if (beta > 0.0)
        {
            MatrixType StiffnessMatrix = ZeroMatrixType( mat_size, mat_size );

            MatrixType FullLeftHandSideMatrix;
            VectorType FullRightHandSideVector;

            bool CalculateStiffnessMatrixFlag = true;
            bool CalculateResidualVectorFlag = false;
            BaseType::CalculateAll( FullLeftHandSideMatrix, FullRightHandSideVector, rCurrentProcessInfo,
                    CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

            for (unsigned int i = 0; i < mat_size; ++i)
            {
                for (unsigned int j = 0; j < mat_size; ++j)
                {
                    StiffnessMatrix(i, j) = FullLeftHandSideMatrix(dim*i + 2, dim*j + 2);
                }
            }

            noalias( rDampMatrix ) += beta * StiffnessMatrix;
        }
    }

    template<class TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::CalculateBoperator( MatrixType& B_Operator, const Vector& N, const MatrixType& DN_DX ) const
    {
        KRATOS_TRY

        B_Operator.clear();

        const unsigned int number_of_nodes = this->GetGeometry().size();

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            B_Operator( 0, i*3     ) = DN_DX( i, 0 );
            B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
            B_Operator( 2, i*3 + 2 ) = 0.0;
            B_Operator( 3, i*3     ) = DN_DX( i, 1 );
            B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
            B_Operator( 4, i*3 + 1 ) = 0.0;
            B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
            B_Operator( 5, i*3     ) = 0.0;
            B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
        }

        KRATOS_CATCH( "" )
    }

    template<class TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::CalculateBBaroperator( MatrixType& B_Operator, const MatrixType& DN_DX, const MatrixType& Bdil_bar ) const
    {
        const unsigned int number_of_nodes = this->GetGeometry().size();

        B_Operator.clear();

        constexpr ValueType R1R3 = (1.0 / 3);
        DataType tmp1;
        DataType tmp2;
        DataType tmp3;
        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            tmp1 = R1R3 * (Bdil_bar( i, 0 ) - DN_DX( i, 0 ));
            tmp2 = R1R3 * (Bdil_bar( i, 1 ) - DN_DX( i, 1 ));
            tmp3 = R1R3 * Bdil_bar( i, 2 );

            B_Operator( 0, i*3 ) = DN_DX( i, 0 ) + tmp1;
            B_Operator( 1, i*3 ) = tmp1;
            B_Operator( 2, i*3 ) = tmp1;

            B_Operator( 0, i*3 + 1)  = tmp2;
            B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 ) + tmp2;
            B_Operator( 2, i*3 + 1 ) = tmp2;

            B_Operator( 0, i*3 + 2)  = tmp3;
            B_Operator( 1, i*3 + 2 ) = tmp3;
            B_Operator( 2, i*3 + 2 ) = tmp3;

            B_Operator( 3, i*3 )     = DN_DX( i, 1 );
            B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
            B_Operator( 4, i*3 + 1 ) = 0.0;
            B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
            B_Operator( 5, i*3 )     = 0.0;
            B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
        }
    }

    template<class TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
    {
        unsigned int mat_size = this->GetGeometry().size();

        if ( rResult.size() != mat_size )
            rResult.resize( mat_size, false );

        for ( unsigned int i = 0 ; i < this->GetGeometry().size() ; ++i )
        {
            rResult[i] = this->GetGeometry()[i].GetDof( VARSELC( DataType, DISPLACEMENT, Z ) ).EquationId();
        }
    }

    template<class TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo ) const
    {
        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0 ; i < this->GetGeometry().size() ; ++i )
        {
            ElementalDofList.push_back( this->GetGeometry()[i].pGetDof( VARSELC( DataType, DISPLACEMENT, Z ) ) );
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::GetValuesVector( VectorType& values, int Step ) const
    {
        const unsigned int number_of_nodes = this->GetGeometry().size();
        unsigned int mat_size = number_of_nodes;

        if ( values.size() != mat_size )
            values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            values[i] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, DISPLACEMENT, Z ), Step );
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::GetFirstDerivativesVector( VectorType& values, int Step ) const
    {
        const unsigned int number_of_nodes = this->GetGeometry().size();
        unsigned int mat_size = number_of_nodes;

        if ( values.size() != mat_size )
            values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            values[i] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, DISPLACEMENT_DT, Z ), Step );
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinearAntiPlane<TNodeType>::GetSecondDerivativesVector( VectorType& values, int Step ) const
    {
        const unsigned int number_of_nodes = this->GetGeometry().size();
        unsigned int mat_size = number_of_nodes;

        if ( values.size() != mat_size )
            values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            values[i] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, ACCELERATION, Z ), Step );
        }
    }

    template class BaseKinematicLinearAntiPlane<RealNode>;
    template class BaseKinematicLinearAntiPlane<ComplexNode>;
    template class BaseKinematicLinearAntiPlane<GComplexNode>;

} // Namespace Kratos

#undef ENABLE_DEBUG_CONSTITUTIVE_LAW
