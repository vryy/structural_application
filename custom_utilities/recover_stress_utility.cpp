//
//   Project Name:        KratosStructuralApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 12 May 2026 $
//
//


// Project includes
#include "utilities/math_utils.h"
#include "recover_stress_utility.h"
#include "structural_application_variables.h"


namespace Kratos
{

template<typename TEntityType>
typename TEntityType::DataType RecoverStressUtility::ComputeZZErrorEstimation(TEntityType& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    typedef typename TEntityType::DataType DataType;

    typedef typename TEntityType::ValueType ValueType;

    typedef typename TEntityType::VectorType VectorType;

    typedef typename TEntityType::MatrixType MatrixType;

    typedef typename TEntityType::GeometryType GeometryType;

    GeometryType& rGeometry = rElement.GetGeometry();

    const IntegrationMethod ThisIntegrationMethod = rElement.GetIntegrationMethod();

    #ifdef ENABLE_BEZIER_GEOMETRY
    rGeometry.Initialize(ThisIntegrationMethod);
    #endif

    const typename GeometryType::IntegrationPointsArrayType& integration_points =
                rGeometry.IntegrationPoints( ThisIntegrationMethod );

    std::vector<DataType> E;
    rElement.CalculateOnIntegrationPoints( VARSEL(DataType, YOUNG_MODULUS), E, rCurrentProcessInfo );

    std::vector<DataType> NU;
    rElement.CalculateOnIntegrationPoints( VARSEL(DataType, POISSON_RATIO), NU, rCurrentProcessInfo );

    std::vector<DataType> J0;
    rElement.CalculateOnIntegrationPoints( VARSEL(DataType, JACOBIAN_0), J0, rCurrentProcessInfo );

    std::vector<VectorType> stress;
    rElement.CalculateOnIntegrationPoints( VARSEL(DataType, STRESSES), stress, rCurrentProcessInfo );

    std::vector<VectorType> rstress;
    rElement.CalculateOnIntegrationPoints( VARSEL(DataType, RECOVERY_STRESSES), rstress, rCurrentProcessInfo );

    DataType result = 0.0;
    MatrixType stress_tensor(3, 3), rstress_tensor(3, 3), ds(3, 3), Eds(3, 3);
    typename SD_MathUtils<DataType>::Fourth_Order_Tensor Ci;

    for(std::size_t PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        ValueType weight = integration_points[PointNumber].Weight();

        SD_MathUtils<DataType>::StressVectorToTensor(stress[PointNumber], stress_tensor);
        SD_MathUtils<DataType>::StressVectorToTensor(rstress[PointNumber], rstress_tensor);

        SD_MathUtils<DataType>::CalculateInversedElasticTensor(Ci, E[PointNumber], NU[PointNumber]);

        noalias(ds) = (rstress_tensor - stress_tensor);
        Eds.clear();
        SD_MathUtils<DataType>::ContractFourthOrderTensor(1.0, Ci, ds, Eds);

        result += SD_MathUtils<DataType>::mat_inner_prod(ds, Eds) * weight * J0[PointNumber];
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    rGeometry.Clean();
    #endif

    DataType error = std::sqrt(result);
    rElement.SetValue(VARSEL( DataType, LOCAL_ERROR ), error);

    return error;
}

template<typename TVariableType>
double RecoverStressUtility::ComputeKellyErrorEstimation(const InterfaceContainer<ModelPart>& icon,
        const TVariableType& rVariable)
{
    const auto& interfaces = icon.GetInterfaces();

    #ifdef DUMP_INTERFACE
    KRATOS_WATCH(icon)
    #endif

    // evaluate the jump in each half face

    double sum = 0.0;

    for (auto& hf : interfaces)
    {
        // integrate the jump on the surface
        if ((hf.first != nullptr) && (hf.second != nullptr))
        {
            auto& face_geom = *(hf.left);

            const unsigned int local_dim = face_geom.LocalSpaceDimension();
            const unsigned int dim = face_geom.WorkingSpaceDimension();

            const auto& geom1 = hf.first->GetGeometry();
            const auto& geom2 = hf.second->GetGeometry();

            const auto ThisIntegrationMethod = face_geom.GetDefaultIntegrationMethod();

            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            face_geom.Initialize(ThisIntegrationMethod);
            #endif

            const auto& integration_points = face_geom.IntegrationPoints( ThisIntegrationMethod );

            const ShapeFunctionsGradientsType& DN_De = face_geom.ShapeFunctionsLocalGradients( ThisIntegrationMethod );

            // const Matrix& Ncontainer = face_geom.ShapeFunctionsValues( ThisIntegrationMethod );

            JacobiansType J;
            J = face_geom.Jacobian( J, ThisIntegrationMethod );

            double DetJ;

            CoordinatesArrayType P;
            LocalCoordinatesArrayType p1, p2;
            Vector t1(dim), t2(dim), n(dim);
            double jump = 0.;
            for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
            {
                // compute the determinant of Jacobian
                DetJ = MathUtils<double>::Det(prod(trans(J[PointNumber]), J[PointNumber]));

                // obtain the physical coordinates of the integration point

                P = face_geom.GlobalCoordinates( P, integration_points[PointNumber] );

                // compute the tangents
                if (local_dim > 0)
                {
                    noalias(t1) = ZeroVector(dim);
                    for ( unsigned int i = 0; i < face_geom.size(); ++i )
                    {
                        t1[0] += face_geom[i].X0() * DN_De[PointNumber]( i, 0 );
                        t1[1] += face_geom[i].Y0() * DN_De[PointNumber]( i, 0 );
                        if(dim == 3)
                            t1[2] += face_geom[i].Z0() * DN_De[PointNumber]( i, 0 );
                    }
                }

                if (local_dim > 1)
                {
                    noalias(t2) = ZeroVector(dim);
                    for ( unsigned int i = 0; i < face_geom.size(); ++i )
                    {
                        t2[0] += face_geom[i].X0() * DN_De[PointNumber]( i, 1 );
                        t2[1] += face_geom[i].Y0() * DN_De[PointNumber]( i, 1 );
                        if(dim == 3)
                            t2[2] += face_geom[i].Z0() * DN_De[PointNumber]( i, 1 );
                    }
                }

                // compute the normal
                if (local_dim == 1)
                {
                    n[0] = -t1[1];
                    n[1] = t1[0];
                    n[2] = 0.;
                }
                else if (local_dim == 2)
                {
                    noalias(n) = MathUtils<double>::CrossProduct(t1, t2);
                }

                n /= norm_2(n);

                // compute the local coordinates in both side
                p1 = geom1.PointLocalCoordinates(p1, P, true, 1e-10);
                p2 = geom2.PointLocalCoordinates(p2, P, true, 1e-10);

                // contribute to the jump
                double tmp = ComputeJump(geom1, geom2, p1, p2, n, rVariable);
                jump += tmp * tmp * DetJ * integration_points[PointNumber].Weight();
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            face_geom.Clean();
            #endif

            // hf.jump = jump;
            sum += jump;

            // sum up the jump over the element
            hf.first->GetValue(LOCAL_ERROR) += jump;
            hf.second->GetValue(LOCAL_ERROR) += jump;
        }
    }

    return sum;
}

Vector RecoverStressUtility::ComputeGradient(const GeometryType& rGeometry, const LocalCoordinatesArrayType& rPoint, const Variable<double>& rVariable)
{
    const unsigned int local_dim = rGeometry.LocalSpaceDimension();

    Matrix DN_De;
    DN_De = rGeometry.ShapeFunctionsLocalGradients(DN_De, rPoint);

    Matrix J;
    J = rGeometry.Jacobian(J, rPoint);

    Matrix InvJ;
    double DetJ;
    MathUtils<double>::InvertMatrix( J, InvJ, DetJ );

    Matrix DN_DX(DN_De.size1(), DN_De.size2());
    noalias(DN_DX) = prod(DN_De, InvJ);

    Vector grad(local_dim);
    noalias(grad) = ZeroVector(local_dim);
    for (unsigned int i = 0; i < rGeometry.size(); ++i)
    {
        for (unsigned int j = 0; j < local_dim; ++j)
        {
            grad(j) += DN_DX(i, j) * rGeometry[i].GetSolutionStepValue(rVariable);
        }
    }

    return grad;
}

Matrix RecoverStressUtility::ComputeGradient(const GeometryType& rGeometry, const LocalCoordinatesArrayType& rPoint, const Variable<array_1d<double, 3> >& rVariable)
{
    const unsigned int local_dim = rGeometry.LocalSpaceDimension();

    Matrix DN_De;
    DN_De = rGeometry.ShapeFunctionsLocalGradients(DN_De, rPoint);

    Matrix J;
    J = rGeometry.Jacobian(J, rPoint);

    Matrix InvJ;
    double DetJ;
    MathUtils<double>::InvertMatrix( J, InvJ, DetJ );

    Matrix DN_DX(DN_De.size1(), DN_De.size2());
    noalias(DN_DX) = prod(DN_De, InvJ);

    Matrix grad(3, local_dim);
    noalias(grad) = ZeroMatrix(3, local_dim);
    for (unsigned int i = 0; i < rGeometry.size(); ++i)
    {
        for (unsigned int j = 0; j < local_dim; ++j)
        {
            noalias( column(grad, j) ) += DN_DX(i, j) * rGeometry[i].GetSolutionStepValue(rVariable);
        }
    }

    return grad;
}

double RecoverStressUtility::ComputeJump(const GeometryType& rGeometry1, const GeometryType& rGeometry2,
            const LocalCoordinatesArrayType& p1, const LocalCoordinatesArrayType& p2,
            const Vector& n, const Variable<double>& rVariable)
{
    auto du1 = ComputeGradient(rGeometry1, p1, rVariable);
    auto du2 = ComputeGradient(rGeometry2, p2, rVariable);
    return inner_prod(du1 - du2, n);
}

double RecoverStressUtility::ComputeJump(const GeometryType& rGeometry1, const GeometryType& rGeometry2,
            const LocalCoordinatesArrayType& p1, const LocalCoordinatesArrayType& p2,
            const Vector& n, const Variable<array_1d<double, 3> >& rVariable)
{
    auto du1 = ComputeGradient(rGeometry1, p1, rVariable);
    auto du2 = ComputeGradient(rGeometry2, p2, rVariable);
    return norm_2(prod(du1 - du2, n));
}

////////////////////

// function template specialization

template Element::DataType RecoverStressUtility::ComputeZZErrorEstimation(Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double RecoverStressUtility::ComputeKellyErrorEstimation(const InterfaceContainer<ModelPart>&, const Variable<double>&);
template double RecoverStressUtility::ComputeKellyErrorEstimation(const InterfaceContainer<ModelPart>&, const Variable<array_1d<double, 3> >&);

#undef DUMP_INTERFACE

}  // namespace Kratos.
