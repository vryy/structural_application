/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 11 Mar 2021 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_STRUCTURAL_APPLICATION_CALCULATE_STRAIN_ENERGY_PROCESS_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APPLICATION_CALCULATE_STRAIN_ENERGY_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "structural_application_variables.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** A class to compute the strain energy of the model_part
*/
class CalculateStrainEnergyProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    typedef Element::GeometryType GeometryType;
    typedef GeometryType::PointType NodeType;
    typedef GeometryType::PointType::PointType PointType;
    typedef ModelPart::ElementsContainerType ElementsContainerType;

    /// Pointer definition of CalculateStrainEnergyProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateStrainEnergyProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// This constructor will take all elements of the model_part for topology optimization
    CalculateStrainEnergyProcess(const ModelPart& r_model_part)
    : mr_model_part(r_model_part), mEnergy(0.0)
    {}

    /// Destructor.
    virtual ~CalculateStrainEnergyProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {

        const ElementsContainerType& rElements = mr_model_part.Elements();

        mEnergy = 0.0;
        for(ElementsContainerType::const_iterator i_element = rElements.begin() ; i_element != rElements.end(); ++i_element)
        {
            // get the integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = i_element->GetGeometry().IntegrationPoints( i_element->GetIntegrationMethod() );

            // compute strain energy (i.e compliance) at the integration points
            std::vector<double> StrainEnergyAtIntegrationPoints;
            i_element->CalculateOnIntegrationPoints(STRAIN_ENERGY, StrainEnergyAtIntegrationPoints, mr_model_part.GetProcessInfo());

            // get the Jacobian at integration points
            std::vector<double> Jacobian;
            i_element->CalculateOnIntegrationPoints(JACOBIAN_0, Jacobian, mr_model_part.GetProcessInfo());

            // compute strain energy at the element
            double StrainEnergy = 0.0;
            for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
            {
                StrainEnergy += integration_points[PointNumber].Weight() * StrainEnergyAtIntegrationPoints[PointNumber] * Jacobian[PointNumber];
            }

            // add to system energy
            mEnergy += StrainEnergy;
        }

        #ifdef SD_APP_FORWARD_COMPATIBILITY
        mr_model_part.GetCommunicator().GetDataCommunicator().SumAll(mEnergy);
        #else
        mr_model_part.GetCommunicator().SumAll(mEnergy);
        #endif

    }

    ///@}
    ///@name Access
    ///@{

    /// Get the computed strain energy
    double GetEnergy() const
    {
        return mEnergy;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CalculateStrainEnergyProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateStrainEnergyProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Total strain energy: " << mEnergy;
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    const ModelPart& mr_model_part;
    double mEnergy;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    CalculateStrainEnergyProcess& operator=(CalculateStrainEnergyProcess const& rOther);

    /// Copy constructor.
    //CalculateStrainEnergyProcess(FindConditionsNeighboursProcess const& rOther);


    ///@}

}; // Class FindConditionsNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_APPLICATION_CALCULATE_STRAIN_ENERGY_PROCESS_H_INCLUDED  defined
