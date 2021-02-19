/*
==============================================================================
KratosR1StructuralApplication
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
//   Last Modified by:    $Author: jelena $
//   Date:                $Date: 2015-09-17 07:11:02 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_TIP_UTILITY_INCLUDED )
#define  KRATOS_TIP_UTILITY_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "geometries/geometry_data.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "custom_conditions/tip_condition.h"


namespace Kratos
{
    /**
     * Tip Utility
     */
    class TipUtility
    {
        public:
            typedef ModelPart::ElementsContainerType ElementsContainerType;
            typedef ModelPart::ConditionsContainerType ConditionsContainerType;
            typedef Element::GeometryType GeometryType;
            typedef GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
            typedef GeometryType::PointType NodeType;
            typedef NodeType::PointType PointType;
            typedef Properties PropertiesType;
            typedef GeometryData::IntegrationMethod IntegrationMethod;
            typedef std::size_t IndexType;

            /**
             * class pointer definition
             */
            KRATOS_CLASS_POINTER_DEFINITION( TipUtility );

            /**
             * Constructor.
             */
            TipUtility()
            {
            }

            /**
             * Destructor.
             */
            virtual ~TipUtility()
            {
            }

            /// Get the last condition id of the model part
            static std::size_t GetLastConditionId(ModelPart& r_model_part)
            {
                std::size_t lastCondId = 0;
                for(typename ModelPart::ConditionsContainerType::ptr_iterator it = r_model_part.Conditions().ptr_begin();
                        it != r_model_part.Conditions().ptr_end(); ++it)
                {
                    if((*it)->Id() > lastCondId)
                        lastCondId = (*it)->Id();
                }

                return lastCondId;
            }

            /**
             * Initializes mesh tying by means of lagrange multipliers
             * tip_elements: list of beam elements
             * tip_soil_elements: list of solid elements tied to the beam
             */
            void InitializeTipUtility( ModelPart& model_part,
                std::vector<unsigned int>& tip_elements,
                std::vector<unsigned int>& tip_soil_elements,
                Properties::Pointer linkProperties )
            {
                ElementsContainerType::Pointer tips( new ElementsContainerType() );
                ElementsContainerType::Pointer tip_soil_elems( new ElementsContainerType() );

                IntegrationMethod ThisIntegrationMethod;
                GeometryType::Pointer tempGeometry;

                std::cout << "Initializing TipUtility..." << std::endl;
                for( unsigned int it = 0; it != tip_elements.size(); it++ )
                {
                    tips->push_back( model_part.pGetElement( tip_elements[it]) );
                }

                for( unsigned int it = 0; it != tip_soil_elements.size(); it++ )
                {
                    tip_soil_elems->push_back( model_part.pGetElement( tip_soil_elements[it]) );
                }

                // loop through all conditions to compute the lastCondId
                std::size_t lastCondId = GetLastConditionId(model_part);
                KRATOS_WATCH(lastCondId)

                PointType TipPoint;
                PointType TipLocalPoint;
                PointType TipSoilLocalPoint;
                PointType TipSoilGlobalPoint;
                Element::Pointer TargetElement;
                std::size_t number_of_tip_conditions = 0;
                for( ElementsContainerType::ptr_iterator it = tips->ptr_begin(); it != tips->ptr_end(); ++it )
                {
                    /******KRATOS_WATCH(it);*/
//                    KRATOS_WATCH(*it);
                    ThisIntegrationMethod = (*it)->GetGeometry().GetDefaultIntegrationMethod();
                    const IntegrationPointsArrayType& integration_points = (*it)->GetGeometry().IntegrationPoints(ThisIntegrationMethod);

                    for( IndexType i = 0; i < integration_points.size(); ++i )
                    {
                        noalias(TipLocalPoint) = integration_points[i];
                        (*it)->GetGeometry().GlobalCoordinates( TipPoint, TipLocalPoint );
                        if( FindPartnerElement( TipPoint, tip_soil_elems, TargetElement, TipSoilLocalPoint ) )
                        {
                            TargetElement->GetGeometry().GlobalCoordinates( TipSoilGlobalPoint, TipSoilLocalPoint );

                            tempGeometry = GeometryType::Pointer( new GeometryType() );
                            Condition::Pointer newLink = Condition::Pointer( new TipCondition( ++lastCondId, tempGeometry, linkProperties, TargetElement, *it, TipSoilLocalPoint, TipLocalPoint ) );

                            model_part.Conditions().push_back( newLink );
                            ++number_of_tip_conditions;
                        }
                    }
                }

                std::cout << "Setup tip conditions completed, " << number_of_tip_conditions << " conditions are added to model_part" << std::endl;
            }//InitializeTipUtility

            /**
             * calculates for a point given with the physical coords newNode
             * the element oldElement where it lays in and the natural coords
             * localPoint within this element
             * @return whether a corresponding element and natural coords could be found
             * @param newNode physical coordinates of given point
             * @param OldMeshElementsArray Array of elements wherein the search should be performed
             * @param oldElement corresponding element for newNode
             * @param rResult corresponding natural coords for newNode
             * TODO: find a faster method for outside search (hextree? etc.), maybe outside this
             * function by restriction of OldMeshElementsArray
             */
            bool FindPartnerElement( PointType& sourcePoint,
                                     const ElementsContainerType::Pointer& TipSoilElements,
                                     Element::Pointer& TargetElement, PointType& rResult)
            {
                bool partner_found= false;
                ElementsContainerType::Pointer TipSoilElementsCandidates( new ElementsContainerType() );
                std::vector<double > OldMinDist;
                bool newMinDistFound= false;

                int counter = 0;
                do
                {
                    double minDist = 1.0e120;
                    newMinDistFound= false;
                    TipSoilElementsCandidates->clear();
                    // (global search)
                    for( ElementsContainerType::ptr_iterator it = TipSoilElements->ptr_begin();
                         it != TipSoilElements->ptr_end(); ++it )
                    {
                        //loop over all nodes in tested element
                        for( unsigned int n=0; n<(*it)->GetGeometry().size(); ++n )
                        {
                            double dist = ((*it)->GetGeometry().GetPoint(n).X0()-sourcePoint[0])
                                        *((*it)->GetGeometry().GetPoint(n).X0()-sourcePoint[0])
                                        +((*it)->GetGeometry().GetPoint(n).Y0()-sourcePoint[1])
                                        *((*it)->GetGeometry().GetPoint(n).Y0()-sourcePoint[1])
                                        +((*it)->GetGeometry().GetPoint(n).Z0()-sourcePoint[2])
                                        *((*it)->GetGeometry().GetPoint(n).Z0()-sourcePoint[2]);
                            if( fabs(dist-minDist) < 1e-7 )
                            {
                                TipSoilElementsCandidates->push_back(*it);
                            }
                            else if( dist < minDist )
                            {
                                bool alreadyUsed= false;
                                for(unsigned int old_dist= 0; old_dist<OldMinDist.size(); ++old_dist)
                                {
                                    if(fabs(dist- OldMinDist[old_dist])< 1e-7 )
                                        alreadyUsed= true;
                                }
                                if(!alreadyUsed)
                                {
                                    TipSoilElementsCandidates->clear();
                                    minDist = dist;
                                    TipSoilElementsCandidates->push_back(*it);
                                    newMinDistFound= true;
                                }
                            }
                        }
                    }

                    OldMinDist.push_back(minDist);

                    for( ElementsContainerType::ptr_iterator it = TipSoilElementsCandidates->ptr_begin();
                         it != TipSoilElementsCandidates->ptr_end(); ++it )
                    {
//                         std::cout << "checking elements list" << std::endl;
                        if( (*it)->GetGeometry().IsInside( sourcePoint, rResult ) )
                        {
//                             std::cout << "isInside" << std::endl;
                            TargetElement=*(it);
                            partner_found=true;
                            return partner_found;
                        }
                    }
//                     std::cout << counter << std::endl;
                    ++counter;
                    if( counter > 27 )
                        break;
                }while(newMinDistFound);

                if(!partner_found)
                    std::cout<<" !!!! NO PARTNER FOUND !!!! "<<std::endl;
                return partner_found;
            }
   /**
    * Calculates for given Loacal coordinates the global coordinates
    * @param Surface surface
    * @param rResult global coordinates
    * @param LocalCoordinates local coordinates
    * @return global coordinates
    */

         GeometryType::CoordinatesArrayType& GlobalCoordinatesTipSoil(Element::Pointer TipSoilElements, GeometryType::CoordinatesArrayType& rResult, GeometryType::CoordinatesArrayType const& LocalCoordinates)
        {
        noalias(rResult)= ZeroVector(3);

        for(IndexType i = 0 ; i < TipSoilElements->GetGeometry().size() ; i++)
        {
            double shape_func= TipSoilElements->GetGeometry().ShapeFunctionValue(i,LocalCoordinates);

            rResult(0) += shape_func*
                ((TipSoilElements->GetGeometry()[i]).X0()
                +(TipSoilElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_X));

            rResult(1) += shape_func*
                ((TipSoilElements->GetGeometry()[i]).Y0()
                +(TipSoilElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Y));

            rResult(2) += shape_func*
                ((TipSoilElements->GetGeometry()[i]).Z0()
                +(TipSoilElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Z));
        }
        return rResult;
        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         GeometryType::CoordinatesArrayType& GlobalCoordinatesTip(const Element::Pointer TipElements, GeometryType::CoordinatesArrayType& rResult, GeometryType::CoordinatesArrayType const& LocalCoordinates)
        {
        noalias(rResult)= ZeroVector(3);

        for(IndexType i = 0 ; i < TipElements->GetGeometry().size() ; i++)
        {
            double shape_func= TipElements->GetGeometry().ShapeFunctionValue(i,LocalCoordinates);

            rResult(0) += shape_func*
                ((TipElements->GetGeometry()[i]).X0()
                +(TipElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_X));

            rResult(1) += shape_func*
                ((TipElements->GetGeometry()[i]).Y0()
                +(TipElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Y));

            rResult(2) += shape_func*
                ((TipElements->GetGeometry()[i]).Z0()
                +(TipElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Z));
        }
        return rResult;
        }
    };//class TipUtility
}  // namespace Kratos.

#endif // KRATOS_PILE_UTILITY_INCLUDED defined
