/*
see license.txt
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 6 Jul 2016 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_EMBEDDED_POINT_TYING_UTILITY_INCLUDED )
#define  KRATOS_EMBEDDED_POINT_TYING_UTILITY_INCLUDED

// System includes

// External includes 
#include "boost/smart_ptr.hpp"
#include "boost/progress.hpp"

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
//#include "spatial_containers/bounding_volume_tree.h"
#include "structural_application_variables.h"


namespace Kratos
{
    template<class TyingLinkType>
    class EmbeddedPointTyingUtility
    {
        public:

            typedef Element::GeometryType   GeometryType;
            typedef Element::NodeType       NodeType;
            typedef NodeType::PointType     PointType;

            /**
             * class pointer definition
             */
            KRATOS_CLASS_POINTER_DEFINITION( EmbeddedPointTyingUtility );

            /**
             * Constructor.
             */
            EmbeddedPointTyingUtility()
            {}

            /**
             * Destructor.
             */
            virtual ~EmbeddedPointTyingUtility()
            {}

            ModelPart::ConditionsContainerType SetUpTyingLinks1( ModelPart& r_model_part, boost::python::list list_slave_elements,
                    boost::python::list list_master_elements )
            {
                ModelPart::ElementsContainerType pSlaveElements;
                ModelPart::ElementsContainerType pMasterElements;

                typedef boost::python::stl_input_iterator<int> iterator_type;
                BOOST_FOREACH(const iterator_type::value_type& id,
                              std::make_pair(iterator_type(list_slave_elements), // begin
                                iterator_type() ) ) // end
                {
                    pSlaveElements.push_back(r_model_part.pGetElement(id));
                }

                typedef boost::python::stl_input_iterator<int> iterator_type;
                BOOST_FOREACH(const iterator_type::value_type& id,
                              std::make_pair(iterator_type(list_master_elements), // begin
                                iterator_type() ) ) // end
                {
                    pMasterElements.push_back(r_model_part.pGetElement(id));
                }

                return this->SetUpTyingLinks2(r_model_part, pSlaveElements, pMasterElements);
            }

            /**
             * Initializes tying conditions
             */
            ModelPart::ConditionsContainerType SetUpTyingLinks2( ModelPart& r_model_part, ModelPart::ElementsContainerType& rpSlaveElements,
                    ModelPart::ElementsContainerType& rpMasterElements )
            {
                std::cout << "Initializing EmbeddedPointTyingUtility..." << std::endl;

                // get the last condition id
                std::size_t lastCondId = 0;
                for( typename ModelPart::ConditionsContainerType::ptr_iterator it = r_model_part.Conditions().ptr_begin();
                        it != r_model_part.Conditions().ptr_end(); ++it)
                {
                    if((*it)->Id() > lastCondId)
                        lastCondId = (*it)->Id();
                }
                KRATOS_WATCH(lastCondId)

                // setup the tying links
                TyingLinkType SampleLink;
                ModelPart::ConditionsContainerType LinkingConditions;
                boost::progress_display show_progress( rpSlaveElements.size() );
                GeometryData::IntegrationMethod ThisIntegrationMethod;
                PointType LocalPoint;
                PointType GlobalPoint;
                for( typename ModelPart::ElementsContainerType::ptr_iterator it = rpSlaveElements.ptr_begin();
                        it != rpSlaveElements.ptr_end(); ++it )
                {
                    // define the list of potential master elements
                    //TODO

                    ThisIntegrationMethod = (*it)->GetGeometry().GetDefaultIntegrationMethod();
                    const GeometryType::IntegrationPointsArrayType& integration_points = (*it)->GetGeometry().IntegrationPoints(ThisIntegrationMethod);

                    for( std::size_t i = 0; i < integration_points.size(); ++i )
                    {
                        noalias(LocalPoint) = integration_points[i];
                        noalias(GlobalPoint) = (*it)->GetGeometry().GlobalCoordinates(GlobalPoint, LocalPoint);

                        Element::Pointer pTargetElement;
                        PointType TargetLocalPoint;
                        if( FindPartnerElement( GlobalPoint, rpMasterElements, pTargetElement, TargetLocalPoint ) )
                        {
                            // Create the linking condition
                            GeometryType::Pointer pTempGeometry = GeometryType::Pointer( new GeometryType() );
                            Condition::Pointer pNewLink = SampleLink.Create(++lastCondId, pTempGeometry, pTargetElement, *it, TargetLocalPoint, LocalPoint);
                            pNewLink->Set(ACTIVE, true);
                            LinkingConditions.push_back( pNewLink );
                        }
                    }
                    ++show_progress;
                }

                for(ConditionsContainerType::ptr_iterator it = LinkingConditions.ptr_begin(); it != LinkingConditions.ptr_end(); ++it)
                    r_model_part.Conditions().push_back(*it);

                std::cout << "Setup point tying links completed, "
                    << LinkingConditions.size() << " linking conditions of " << typeid(SampleLink).name()
                    << " was added to model_part" << std::endl;

                return LinkingConditions;
            }

        private:
            /**
             * calculates ,for a given node with the physical coords, a projected newNode
             * within the solid element where it lays in global/natural coords
             * @return whether a node lays within corresponding element
             * @param newNode physical coordinates of given point
             * @param OldMeshElementsArray Array of elements wherein the search should be performed
             * @param oldElement corresponding element for newNode
             * @param rResult corresponding natural coords for newNode
             * TODO: find a faster method for outside search (hextree? etc.), maybe outside this
             * function by restriction of OldMeshElementsArray
             */
            bool FindPartnerElement( const PointType& rSourcePoint, ModelPart::ElementsContainerType& pMasterElements,
                                     Element::Pointer& pTargetElement, PointType& rLocalTargetPoint )
            {
                ModelPart::ElementsContainerType pMasterElementsCandidates;

                // find the potential master elements
                pMasterElementsCandidates.clear();
                for( typename ModelPart::ElementsContainerType::ptr_iterator it = pMasterElements.ptr_begin();
                        it != pMasterElements.ptr_end(); ++it )
                {
                    // compute the center of the element
                    PointType Center(0.0, 0.0, 0.0);
                    for( unsigned int node = 0; node < (*it)->GetGeometry().size(); ++node )
                        noalias(Center) += (*it)->GetGeometry()[node];
                    Center /= (*it)->GetGeometry().size();

                    // compute the maximum distance from center to vertices
                    double max_dist = 0.0;
                    for( unsigned int node = 0; node < (*it)->GetGeometry().size(); ++node )
                    {
                        double dist = norm_2((*it)->GetGeometry()[node] - Center);
                        if(dist > max_dist)
                            max_dist = dist;
                    }

                    // compute the distance of the source point to the center of the element
                    double sdist = norm_2(rSourcePoint - Center);
                    if(sdist < max_dist)
                        pMasterElementsCandidates.push_back(*it);
                }

                for( typename ModelPart::ElementsContainerType::ptr_iterator it = pMasterElementsCandidates.ptr_begin();
                        it != pMasterElementsCandidates.ptr_end(); ++it )
                {
                    bool is_inside = (*it)->GetGeometry().IsInside( rSourcePoint, rLocalTargetPoint );
                    if( is_inside )
                    {
                        pTargetElement = *it;
                        return true;
                    }
                }

                std::cout << " !!!! WARNING: NO ELEMENT FOUND TO CONTAIN " << rSourcePoint << " !!!! " << std::endl;
                return false;
            }

    };//class EmbeddedPointTyingUtility
}  // namespace Kratos.
#endif // KRATOS_EMBEDDED_NODE_TYING_UTILITY_INCLUDED

