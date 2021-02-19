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
//   Last Modified by:    $Author: hurga $
//   Date:                $Date: 2008-09-17 07:11:02 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_FOUNDATION_UTILITY_INCLUDED )
#define  KRATOS_FOUNDATION_UTILITY_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "structural_application_variables.h"
//#include "custom_conditions/tip_condition.h"
#include "custom_conditions/foundation_condition.h"
#include "utilities/math_utils.h"
/////#include "structural_application_variables.h"
#include "geometries/point_3d.h"


namespace Kratos
{
    /**
     * Foundation Utility
     */
    class FoundationUtility
    {
        public:
            typedef ModelPart::ElementsContainerType ElementsArrayType;
            typedef ModelPart::ConditionsContainerType ConditionsArrayType;
            typedef Element::GeometryType GeometryType;
            typedef GeometryType::PointType NodeType;
            typedef NodeType::PointType PointType;
            typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;
            typedef GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
            typedef Properties PropertiesType;
            typedef std::size_t IndexType;

            /**
             * class pointer definition
             */
            KRATOS_CLASS_POINTER_DEFINITION( FoundationUtility );

            /**
             * Constructor.
             */
            FoundationUtility()
            {
            }

            /**
             * Destructor.
             */
            virtual ~FoundationUtility()
            {
            }

            /**
             * Initializes mesh tying by means of lagrange multipliers
             */
            void InitializeFoundationUtility( ModelPart& model_part,
                std::vector<unsigned int>& foundation_elements,
                std::vector<unsigned int>& soil_elements,
                Properties::Pointer linkProperties )
            {
                ElementsArrayType::Pointer foundations( new ElementsArrayType() );
                ElementsArrayType::Pointer soil_elems( new ElementsArrayType() );
//        KRATOS_WATCH("foundation_utility, line 117");
//                 Geometry<NodeType >::Pointer tempGeometry;
//                 tempGeometry = geometry_object.Create( temp_points );
                // NodeType point( 0.0, 0.0, 0.0 );
                // GeometryType::Pointer tempGeometry = GeometryType::Pointer(new Point3D<NodeType>(point) );
                GeometryType::Pointer tempGeometry = GeometryType::Pointer(new GeometryType());
       //         KRATOS_WATCH( tempGeometry );
        //        KRATOS_WATCH( *tempGeometry );
//        KRATOS_WATCH("foundation_utiliy, line 120");

                std::cout << "Initializing FoundationUtility..." << std::endl;
                std::size_t nconds = 0;
                for( unsigned int it = 0; it != foundation_elements.size(); it++ )
                {
                    foundations->push_back( model_part.pGetElement( foundation_elements[it]) );
                }

                for( unsigned int it = 0; it != soil_elements.size(); it++ )
                {
                    soil_elems->push_back( model_part.pGetElement( soil_elements[it]) );
                }
                for( ElementsArrayType::ptr_iterator it = foundations->ptr_begin();
                                 it != foundations->ptr_end(); ++it )
                {
                    /******KRATOS_WATCH(it);*/
                    for( IndexType i = 0; i < (*it)->GetGeometry().IntegrationPoints().size(); i++ )
                    {
                        PointType FoundationPoint;
                        PointType FoundationLocalPoint = (*it)->GetGeometry().IntegrationPoints()[i];
                        (*it)->GetGeometry().GlobalCoordinates( FoundationPoint, FoundationLocalPoint );
                        PointType SoilLocalPoint;
                        Element::Pointer TargetElement;
                        if( FindPartnerElement( FoundationPoint, soil_elems, TargetElement, SoilLocalPoint ) )
                        {
                        ///    KRATOS_WATCH(i);
                        ///    KRATOS_WATCH("Soil element found:");
                        ///    KRATOS_WATCH(TargetElement->Id());
                        ///    KRATOS_WATCH(FoundationLocalPoint);
                        ///    KRATOS_WATCH(SoilLocalPoint);
                        ///    KRATOS_WATCH(FoundationPoint);
                            PointType SoilGlobalPoint;
                            TargetElement->GetGeometry().GlobalCoordinates( SoilGlobalPoint, SoilLocalPoint );
                            //j +=  TargetElement->GetGeometry().IntegrationPoints().size();
                    //        KRATOS_WATCH( SoilGlobalPoint );
                            int a = (model_part.Conditions().end()-1)->Id();
                        //    KRATOS_WATCH (a);
                            IndexType newId = a + 1;
                    //        KRATOS_WATCH( newId );
        //                    unsigned int j;
        //                    for( IndexType j = 0; j < ( TargetElement->GetGeometry().IntegrationPoints().size()); j++ )
        //                    {
                            Condition::Pointer newLink = Condition::Pointer( new FoundationCondition( newId, tempGeometry, linkProperties, TargetElement, *it, SoilLocalPoint, FoundationLocalPoint) );
                            model_part.Conditions().push_back( newLink );
                            ++nconds;
                    //        KRATOS_WATCH (i);
        //                    }
                        }
                    }
                }

                std::cout << "FoundationUtility completed, " << nconds << " FoundationCondition is added to the model_part" << std::endl;

                return;
            } // InitializeFoundationUtility

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
                                     const ElementsArrayType::Pointer& FoundationSoilElements,
                                     Element::Pointer& TargetElement, PointType& rResult)
            {
                bool partner_found= false;
                ElementsArrayType::Pointer SoilElementsCandidates( new ElementsArrayType() );
                std::vector<double > OldMinDist;
                bool newMinDistFound= false;

////std::cout << "line 187" << std::endl;
                int counter = 0;
                do
                {
                    double minDist = 1.0e120;
                    newMinDistFound= false;
                    SoilElementsCandidates->clear();
                // (global search)
                    for( ElementsArrayType::ptr_iterator it = FoundationSoilElements->ptr_begin();
                         it != FoundationSoilElements->ptr_end(); ++it )
                    {
            //loop over all nodes in tested element
                        for( unsigned int n=0; n<(*it)->GetGeometry().size(); n++ )
                        {
                            double dist = ((*it)->GetGeometry().GetPoint(n).X0()-sourcePoint[0])
                                        *((*it)->GetGeometry().GetPoint(n).X0()-sourcePoint[0])
                                        +((*it)->GetGeometry().GetPoint(n).Y0()-sourcePoint[1])
                                        *((*it)->GetGeometry().GetPoint(n).Y0()-sourcePoint[1])
                                        +((*it)->GetGeometry().GetPoint(n).Z0()-sourcePoint[2])
                                        *((*it)->GetGeometry().GetPoint(n).Z0()-sourcePoint[2]);
                            if( fabs(dist-minDist) < 1e-7 )
                            {
                                SoilElementsCandidates->push_back(*it);
                            }
                            else if( dist < minDist )
                            {
                                bool alreadyUsed= false;
                                for(unsigned int old_dist= 0; old_dist<OldMinDist.size(); old_dist++)
                                {
                                    if(fabs(dist- OldMinDist[old_dist])< 1e-7 )
                                        alreadyUsed= true;
                                }
                                if(!alreadyUsed)
                                {
                                    SoilElementsCandidates->clear();
                                    minDist = dist;
                                    SoilElementsCandidates->push_back(*it);
                                    newMinDistFound= true;
                                }
                            }
                        }
                    }

                    OldMinDist.push_back(minDist);
//                     KRATOS_WATCH(OldElementsSet->size());

                    for( ElementsArrayType::ptr_iterator it = SoilElementsCandidates->ptr_begin();
                         it != SoilElementsCandidates->ptr_end(); ++it )
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
                    counter++;
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

         GeometryType::CoordinatesArrayType& GlobalCoordinatesFoundationSoil(Element::Pointer FoundationSoilElements, GeometryType::CoordinatesArrayType& rResult, GeometryType::CoordinatesArrayType const& LocalCoordinates)
        {
        noalias(rResult)= ZeroVector(3);

        for(IndexType i = 0 ; i < FoundationSoilElements->GetGeometry().size() ; i++)
        {
            double shape_func= FoundationSoilElements->GetGeometry().ShapeFunctionValue(i,LocalCoordinates);

            rResult(0) += shape_func*
                ((FoundationSoilElements->GetGeometry()[i]).X0()
                +(FoundationSoilElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_X));

            rResult(1) += shape_func*
                ((FoundationSoilElements->GetGeometry()[i]).Y0()
                +(FoundationSoilElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Y));

            rResult(2) += shape_func*
                ((FoundationSoilElements->GetGeometry()[i]).Z0()
                +(FoundationSoilElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Z));
        }
        return rResult;
        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         GeometryType::CoordinatesArrayType& GlobalCoordinatesFoundation(const Element::Pointer FoundationElements, GeometryType::CoordinatesArrayType& rResult, GeometryType::CoordinatesArrayType const& LocalCoordinates)
        {
        noalias(rResult)= ZeroVector(3);

        for(IndexType i = 0 ; i < FoundationElements->GetGeometry().size() ; i++)
        {
            double shape_func= FoundationElements->GetGeometry().ShapeFunctionValue(i,LocalCoordinates);

            rResult(0) += shape_func*
                ((FoundationElements->GetGeometry()[i]).X0()
                +(FoundationElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_X));

            rResult(1) += shape_func*
                ((FoundationElements->GetGeometry()[i]).Y0()
                +(FoundationElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Y));

            rResult(2) += shape_func*
                ((FoundationElements->GetGeometry()[i]).Z0()
                +(FoundationElements->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Z));
        }
        return rResult;
        }
    };//class TipUtility
}  // namespace Kratos.

#endif // KRATOS_PILE_UTILITY_INCLUDED defined
