/*
==============================================================================
see soil_mechanics_appplication/LICENSE.txt
==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Giang Bui-Hoang $
//   Date:                $Date: 22 Feb 2024 $
//   Revision:            $Revision: 1.0 $
//
//



#if !defined(KRATOS_SOIL_MECHANICS_MULTI_SURFACE_ELASTOPLASTICITY_LAW_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_MULTI_SURFACE_ELASTOPLASTICITY_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "plasticity_laws/multi_surface_plasticity_law.h"


namespace Kratos
{

/**
 * Implementation of the multi yield surface constitutive law
 * REF: gen_plas.pdf, hs_small_dev.pdf
 */
class KRATOS_API(STRUCTURAL_APPLICATION) MultiSurfaceElastoplasticityLaw : public MultiSurfacePlasticityLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(MultiSurfaceElastoplasticityLaw);
    typedef MultiSurfacePlasticityLaw BaseType;
    typedef BaseType::Third_Order_Tensor Third_Order_Tensor;
    typedef BaseType::Fourth_Order_Tensor Fourth_Order_Tensor;

    /**
     * Constructor.
     */
    MultiSurfaceElastoplasticityLaw() : BaseType()
    {}

    /**
     * Destructor.
     */
    virtual ~MultiSurfaceElastoplasticityLaw()
    {}

    /**
     * Operations
     */

    ///////////////////////////////////////////////////////////
    ///////////// PLASTIC CALCULATION SUBROUTINES /////////////
    ///////////////////////////////////////////////////////////

    /// Perform the plastic integration
    int PlasticIntegration(const std::vector<int>& active_surfaces,
        Matrix& stress, std::vector<Vector>& q, std::vector<Vector>& alpha,
        std::vector<double>& dlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props,
        const int debug_level = 0) const;

    /// Compute the (consistent) plastic tangent according to the procedure in plastic integration
    void ComputeConsistentPlasticTangent(const std::vector<int>& active_surfaces,
        Fourth_Order_Tensor& Cep, const Fourth_Order_Tensor& Ce,
        const Matrix& stress, const std::vector<Vector>& q,
        const std::vector<Vector>& alpha, const std::vector<double>& dlambda,
        const ProcessInfo& CurrentProcessInfo, const Properties& props,
        const int debug_level = 0) const;

    ///////////////////////////////////////////////////////////

    /**
     * Turn back information as a string.
     */
    std::string Name() const override
    {
        std::stringstream ss;
        ss << "MultiSurfaceElastoplasticityLaw<";
        for (std::size_t i = 0; i < mpPlasticityLaws.size(); ++i)
            ss << mpPlasticityLaws[i]->Name() << ",";
        ss << ">";
        return ss.str();
    }

private:

    void PlasticIntegration_ComputeRHS(Vector& rhs, const std::vector<int>& active_surfaces,
        const Matrix& stress, const std::vector<Vector>& q, const std::vector<Vector>& alpha,
        const std::vector<double>& lambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    void PlasticIntegration_ComputeNumLHS(Matrix& lhs, const std::vector<int>& active_surfaces,
        const Matrix& stress, const std::vector<Vector>& q, const std::vector<Vector>& alpha,
        const std::vector<double> dlambda, const double ddlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

}; /* Class MultiSurfaceElastoplasticityLaw */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_MULTI_SURFACE_ELASTOPLASTICITY_LAW_H_INCLUDED  defined */
