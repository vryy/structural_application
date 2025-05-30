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
//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 11 Mar 2025 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes


// External includes


// Project includes
#include "geometries/point_3d.h"
#include "custom_conditions/elastic_line_springs.h"
#include "structural_application_variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
template<typename TNodeType>
BaseElasticLineSprings<TNodeType>::BaseElasticLineSprings(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

template<typename TNodeType>
BaseElasticLineSprings<TNodeType>::BaseElasticLineSprings(IndexType NewId, typename GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

template<typename TNodeType>
BaseElasticLineSprings<TNodeType>::BaseElasticLineSprings(IndexType NewId, typename TNodeType::Pointer const& pNode, typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, typename GeometryType::Pointer( new Point3D<TNodeType>( pNode ) ), pProperties)
{
}

template<typename TNodeType>
BaseElasticLineSprings<TNodeType>::~BaseElasticLineSprings()
{
}

//************************************************************************************
//************************************************************************************
template<typename TNodeType>
typename BaseElasticLineSprings<TNodeType>::BaseType::Pointer BaseElasticLineSprings<TNodeType>::Create(IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const
{
    return typename BaseType::Pointer(new BaseElasticLineSprings(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template<typename TNodeType>
typename BaseElasticLineSprings<TNodeType>::BaseType::Pointer BaseElasticLineSprings<TNodeType>::Create(IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const
{
    return typename BaseType::Pointer(new BaseElasticLineSprings(NewId, pGeom, pProperties));
}

//************************************************************************************
//************************************************************************************
template<typename TNodeType>
GeometryData::IntegrationMethod BaseElasticLineSprings<TNodeType>::GetIntegrationMethod() const
{
    if(this->Has( INTEGRATION_ORDER ))
    {
        if(this->GetValue(INTEGRATION_ORDER) == 1)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 2)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 3)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 4)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 5)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_ERROR << Info() << " does not support for integration order " << this->GetValue(INTEGRATION_ORDER);
    }
    else if(this->GetProperties().Has( INTEGRATION_ORDER ))
    {
        if(this->GetProperties()[INTEGRATION_ORDER] == 1)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if(this->GetProperties()[INTEGRATION_ORDER] == 2)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if(this->GetProperties()[INTEGRATION_ORDER] == 3)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if(this->GetProperties()[INTEGRATION_ORDER] == 4)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if(this->GetProperties()[INTEGRATION_ORDER] == 5)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_ERROR << Info() << " does not support for integration order " << this->GetProperties()[INTEGRATION_ORDER];
    }
    else
        return this->GetGeometry().GetDefaultIntegrationMethod(); // default method
}

//************************************************************************************
//************************************************************************************
template<typename TNodeType>
void BaseElasticLineSprings<TNodeType>::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    KRATOS_CATCH("")
}

template<typename TNodeType>
void BaseElasticLineSprings<TNodeType>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    unsigned int number_of_nodes = this->GetGeometry().size();
    unsigned int dim = 2;

    if(rMassMatrix.size1() != number_of_nodes*dim || rMassMatrix.size2() != number_of_nodes*dim)
        rMassMatrix.resize(number_of_nodes*dim, number_of_nodes*dim, false);
    noalias(rMassMatrix) = ZeroMatrix(number_of_nodes*dim, number_of_nodes*dim);
}

template<typename TNodeType>
void BaseElasticLineSprings<TNodeType>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    unsigned int number_of_nodes = this->GetGeometry().size();
    unsigned int dim = 2;

    if(rDampingMatrix.size1() != number_of_nodes*dim || rDampingMatrix.size2() != number_of_nodes*dim)
        rDampingMatrix.resize(number_of_nodes*dim, number_of_nodes*dim, false);
    noalias(rDampingMatrix) = ZeroMatrix(number_of_nodes*dim, number_of_nodes*dim);
}

//************************************************************************************
//************************************************************************************
template<typename TNodeType>
void BaseElasticLineSprings<TNodeType>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
template<typename TNodeType>
void BaseElasticLineSprings<TNodeType>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                                      const ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag,
                                                      bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    unsigned int number_of_nodes = this->GetGeometry().size();
    unsigned int dim = 2;

    //resizing LHS and RHS where needed
    if (CalculateStiffnessMatrixFlag)
    {
        if(rLeftHandSideMatrix.size1() != number_of_nodes*dim || rLeftHandSideMatrix.size2() != number_of_nodes*dim)
            rLeftHandSideMatrix.resize(number_of_nodes*dim, number_of_nodes*dim, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrixType(number_of_nodes*dim, number_of_nodes*dim);
    }

    if (CalculateResidualVectorFlag)
    {
        if(rRightHandSideVector.size() != number_of_nodes*dim)
            rRightHandSideVector.resize(number_of_nodes*dim, false);
        noalias(rRightHandSideVector) = ZeroVectorType(number_of_nodes*dim);
    }

    if( number_of_nodes == 1 )
    {
        const auto& bedding = this->GetGeometry()[0].GetSolutionStepValue(ELASTIC_BEDDING_STIFFNESS);
        const auto& displacement = this->GetGeometry()[0].GetSolutionStepValue(VARSEL(DataType, DISPLACEMENT));
        for( unsigned int i = 0; i < 3; i++ )
        {
            rLeftHandSideMatrix(i, i) = bedding[i];
            rRightHandSideVector[i] = -displacement[i]*bedding[i];
        }
    }
    else
    {
        const GeometryData::IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        this->GetGeometry().Initialize(ThisIntegrationMethod);
        #endif

//        std::cout << "Displacements:" << std::endl;
//        for ( unsigned int node = 0; node < this->GetGeometry().size(); node++ )
//            std::cout << " " << this->GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT ) << std::endl;

        //reading integration points and local gradients
        const typename GeometryType::IntegrationPointsArrayType& integration_points = this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );
        const typename GeometryType::ShapeFunctionsGradientsType& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients( ThisIntegrationMethod );
        const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

        bool use_distributed_properties = false;
        if (this->Has(USE_DISTRIBUTED_PROPERTIES))
            use_distributed_properties = this->GetValue(USE_DISTRIBUTED_PROPERTIES);

        array_1d<DataType, 3> bedding;
        if (use_distributed_properties)
        {
            if (!this->Has(ELASTIC_BEDDING_STIFFNESS))
                KRATOS_ERROR << "ELASTIC_BEDDING_STIFFNESS is not set, this is required since condition " << this->Id()
                             << " USE_DISTRIBUTED_PROPERTIES flag is set";
            noalias(bedding) = this->GetValue(ELASTIC_BEDDING_STIFFNESS);
        }
        else
        {
            noalias(bedding) = this->GetProperties()[ELASTIC_BEDDING_STIFFNESS];
        }

        // bedding[0]: normal springs component
        // bedding[1]: tangent springs component

        VectorType tangent(dim), normal_vector(dim);
        VectorType displacement(dim);
        DataType norm_t, dA;
        DataType du_n, du_t;
        Vector N(number_of_nodes);

//        KRATOS_WATCH(ThisIntegrationMethod)
//        KRATOS_WATCH(integration_points.size())

        for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            noalias(N) = row(Ncontainer, PointNumber);

            noalias(displacement) = ZeroVectorType(dim);

            noalias(tangent) = ZeroVectorType(dim);

            for ( unsigned int node = 0; node < number_of_nodes; ++node )
            {
                tangent += subrange(this->GetGeometry()[node].GetInitialPosition(), 0, dim) * DN_De[PointNumber]( node, 0 );

                noalias( displacement ) += subrange(this->GetGeometry()[node].GetSolutionStepValue( VARSEL(DataType, DISPLACEMENT) ), 0, dim) * N[node];
            }

            normal_vector[0] = -tangent[1];
            normal_vector[1] = tangent[0];

            DataType IntegrationWeight = integration_points[PointNumber].Weight();
            //modify integration weight in case of 2D
            IntegrationWeight *= this->GetProperties()[THICKNESS];

            // normalize
            norm_t = norm_2(tangent);
            dA = norm_t;

            tangent /= norm_t;
            normal_vector /= norm_t;

            if( CalculateResidualVectorFlag )
            {
                du_n = inner_prod(displacement, normal_vector);
                du_t = inner_prod(displacement, tangent);

                for ( unsigned int prim = 0; prim < this->GetGeometry().size(); ++prim )
                {
                    subrange(rRightHandSideVector, prim*dim, (prim+1)*dim) -=
                        (du_n*bedding[0]*normal_vector + du_t*bedding[1]*tangent)
                            * N[prim] * IntegrationWeight * dA;
                }
            }

            if( CalculateStiffnessMatrixFlag )
            {
                for ( unsigned int prim = 0; prim < this->GetGeometry().size(); ++prim )
                {
                    for ( unsigned int sec = 0; sec < this->GetGeometry().size(); ++sec )
                    {
                        subrange(rLeftHandSideMatrix, prim*dim, (prim+1)*dim, sec*dim, (sec+1)*dim) +=
                            (bedding[0]*outer_prod(normal_vector, normal_vector) + bedding[1]*outer_prod(tangent, tangent) )
                          * N[prim] * N[sec] * IntegrationWeight * dA;
                    }
                }
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //clean the internal data of the geometry
        this->GetGeometry().Clean();
        #endif

       // if( CalculateResidualVectorFlag )
       //     KRATOS_WATCH( rRightHandSideVector )

       // if( CalculateStiffnessMatrixFlag )
       //     KRATOS_WATCH( rLeftHandSideMatrix )
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
template<typename TNodeType>
void BaseElasticLineSprings<TNodeType>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int number_of_nodes = this->GetGeometry().size();
    unsigned int index;
    unsigned int dim = 2;

    if(rResult.size() != number_of_nodes*dim)
        rResult.resize(number_of_nodes*dim);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        index = i*dim;
        rResult[index] = (this->GetGeometry()[i].GetDof(VARSELC(DataType, DISPLACEMENT, X)).EquationId());
        rResult[index+1] = (this->GetGeometry()[i].GetDof(VARSELC(DataType, DISPLACEMENT, Y)).EquationId());
    }
}

//************************************************************************************
//************************************************************************************
template<typename TNodeType>
void BaseElasticLineSprings<TNodeType>::GetDofList(DofsVectorType& rConditionalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    unsigned int number_of_nodes = this->GetGeometry().size();
    unsigned int index;
    unsigned int dim = 2;

    if(rConditionalDofList.size() != number_of_nodes*dim)
        rConditionalDofList.resize(number_of_nodes*dim);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        index = i*dim;
        rConditionalDofList[index] = (this->GetGeometry()[i].pGetDof(VARSELC(DataType, DISPLACEMENT, X)));
        rConditionalDofList[index+1] = (this->GetGeometry()[i].pGetDof(VARSELC(DataType, DISPLACEMENT, Y)));
    }
}

//************************************************************************************
//************************************************************************************
template class BaseElasticLineSprings<RealNode>;
template class BaseElasticLineSprings<ComplexNode>;

} // Namespace Kratos
