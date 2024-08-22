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



#if !defined(KRATOS_SOIL_MECHANICS_MULTI_SURFACE_PLASTICITY_LAW_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_MULTI_SURFACE_PLASTICITY_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "plasticity_laws/general_plasticity_law.h"


namespace Kratos
{

/**
 * Implementation of the multi yield surface constitutive law
 * REF: gen_plas.pdf, hs_small_dev.pdf
 */
class MultiSurfacePlasticityLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(MultiSurfacePlasticityLaw);
    typedef GeneralPlasticityLaw::Third_Order_Tensor Third_Order_Tensor;
    typedef GeneralPlasticityLaw::Fourth_Order_Tensor Fourth_Order_Tensor;

    /**
     * Constructor.
     */
    MultiSurfacePlasticityLaw()
    {}

    /**
     * Destructor.
     */
    virtual ~MultiSurfacePlasticityLaw()
    {}

    /**
     * Operations
     */

    void SetLaw(const int i, GeneralPlasticityLaw::Pointer pLaw)
    {
        if (i >= mpPlasticityLaws.size() || i < 0)
            KRATOS_ERROR << "Invalid plastic surface index " << i;
        mpPlasticityLaws[i] = pLaw;
    }

    void AddLaw(GeneralPlasticityLaw::Pointer pLaw)
    {
        mpPlasticityLaws.push_back(pLaw);
    }

    GeneralPlasticityLaw::Pointer pGetLaw(const unsigned int i) const
    {
        return mpPlasticityLaws[i];
    }

    /// Check which yield surfaces are active
    virtual std::vector<int> CheckActiveSurfaces(const Matrix& stress, const std::vector<Vector>& q,
        const std::vector<Vector>& alpha, const double FTOL,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /**
     * Turn back information as a string.
     */
    virtual std::string Name() const
    {
        std::stringstream ss;
        ss << "MultiSurfacePlasticityLaw<";
        for (std::size_t i = 0; i < mpPlasticityLaws.size(); ++i)
            ss << mpPlasticityLaws[i]->Name() << ",";
        ss << ">";
        return ss.str();
    }

    /**
     * Print information about this object.
     */
    virtual void Print(std::ostream& rOStream) const
    {
        rOStream << Name();
    }

protected:

    std::vector<GeneralPlasticityLaw::Pointer> mpPlasticityLaws;

private:

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        // TODO save the array of plasticity laws
        // for (auto pLaw : mpPlasticityLaws)
        //     pLaw->save(rSerializer);
    }

    virtual void load(Serializer& rSerializer)
    {
        // TODO load the array of plasticity laws
        // for (auto pLaw : mpPlasticityLaws)
        //     pLaw->load(rSerializer);
    }

    ///@}

}; /* Class MultiSurfacePlasticityLaw */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_MULTI_SURFACE_PLASTICITY_LAW_H_INCLUDED  defined */
