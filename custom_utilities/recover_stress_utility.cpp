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

typedef RecoverStressUtility::HalfFace HalfFace;
typedef RecoverStressUtility::Comparator Comparator;

template<typename TVariableType>
double RecoverStressUtility::ComputeKellyErrorEstimation(ModelPart& rModelPart,
        const TVariableType& rVariable)
{
    KRATOS_TRY

    // construct the half face structure
    auto half_face_set = ConstructHalfFaceStructure(rModelPart.Elements());

    // evaluate the jump in each half face

    double sum = 0.0;

    for (auto& hf : half_face_set)
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

            const Matrix& Ncontainer = face_geom.ShapeFunctionsValues( ThisIntegrationMethod );

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

            hf.jump = jump;
            sum += jump;
        }
    }

    // sum up the jump over the element

    for (auto it = rModelPart.Elements().begin(); it != rModelPart.Elements().end(); ++it)
    {
        auto faces = it->GetGeometry().Faces();

        double ejump = 0.0;
        for (std::size_t i = 0; i < faces.size(); ++i)
        {
            HalfFace hf;
            hf.face = faces[i];

            auto itf = half_face_set.find(hf);
            if (itf != half_face_set.end())
            {
                ejump += itf->jump;
            }
        }

        it->SetValue(LOCAL_ERROR, ejump);
    }

    return sum;

    KRATOS_CATCH("")
}

std::set<HalfFace, Comparator> RecoverStressUtility::ConstructHalfFaceStructure(const ElementsContainerType& rElements)
{
    std::set<HalfFace, Comparator> half_face_set;

    for (auto it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
    {
        auto faces = (*it)->GetGeometry().Faces();

        for (std::size_t i = 0; i < faces.size(); ++i)
        {
            HalfFace hf;
            hf.face = faces[i];

            auto itf = half_face_set.find(hf);
            if (itf == half_face_set.end())
            {
                hf.left = faces(i);
                hf.first = *it;
                half_face_set.insert(hf);
            }
            else
            {
                hf.left = itf->left;
                hf.first = itf->first;
                hf.right = faces(i);
                hf.second = *it;
                half_face_set.erase(itf);
                half_face_set.insert(hf);
            }
        }
    }

    return half_face_set;
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

template double RecoverStressUtility::ComputeKellyErrorEstimation(ModelPart&, const Variable<double>&);
template double RecoverStressUtility::ComputeKellyErrorEstimation(ModelPart&, const Variable<array_1d<double, 3> >&);

}  // namespace Kratos.
