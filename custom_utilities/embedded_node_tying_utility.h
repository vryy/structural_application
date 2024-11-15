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


#if !defined(KRATOS_EMBEDDED_NODE_TYING_UTILITY_INCLUDED )
#define  KRATOS_EMBEDDED_NODE_TYING_UTILITY_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/legacy_structural_app_vars.h"
#include "utilities/math_utils.h"
#include "utilities/progress.h"
#include "geometries/point_3d.h"
//#include "spatial_containers/bounding_volume_tree.h"
#include "structural_application_variables.h"


namespace Kratos
{
    template<class TyingLinkType>
    class EmbeddedNodeTyingUtility
    {
        public:

            typedef Element::GeometryType   GeometryType;
            typedef Element::NodeType       NodeType;
            typedef NodeType::PointType     PointType;

            /**
             * class pointer definition
             */
            KRATOS_CLASS_POINTER_DEFINITION( EmbeddedNodeTyingUtility );

            /**
             * Constructor.
             */
            EmbeddedNodeTyingUtility()
            {}

            /**
             * Destructor.
             */
            virtual ~EmbeddedNodeTyingUtility()
            {}

            ModelPart::ConditionsContainerType SetUpTyingLinks1( ModelPart& r_model_part, boost::python::list list_nodes,
                    boost::python::list list_elements )
            {
                ModelPart::NodesContainerType pNodes;
                ModelPart::ElementsContainerType pElements;

                typedef boost::python::stl_input_iterator<int> iterator_type;
                BOOST_FOREACH(const iterator_type::value_type& id,
                              std::make_pair(iterator_type(list_nodes), // begin
                                iterator_type() ) ) // end
                {
                    pNodes.push_back(r_model_part.pGetNode(id));
                }

                typedef boost::python::stl_input_iterator<int> iterator_type;
                BOOST_FOREACH(const iterator_type::value_type& id,
                              std::make_pair(iterator_type(list_elements), // begin
                                iterator_type() ) ) // end
                {
                    pElements.push_back(r_model_part.pGetElement(id));
                }

                return this->SetUpTyingLinks2(r_model_part, pNodes, pElements);
            }

            /**
             * Initializes tying conditions
             */
            ModelPart::ConditionsContainerType SetUpTyingLinks2( ModelPart& r_model_part, ModelPart::NodesContainerType& rpNodes,
                    ModelPart::ElementsContainerType& rpElements )
            {
                std::cout << "Initializing EmbeddedNodeTyingUtility..." << std::endl;

                // get the last condition id
                std::size_t lastCondId = 0;
                for( typename ModelPart::ConditionsContainerType::ptr_iterator it = r_model_part.Conditions().ptr_begin();
                        it != r_model_part.Conditions().ptr_end(); ++it)
                {
                    if((*it)->Id() > lastCondId)
                        lastCondId = (*it)->Id();
                }
                KRATOS_WATCH(lastCondId)

//                // build the Bounding Volume Tree to dicretize the elements
//                BoundingVolumePartitioner<ModelPart::ElementsContainerType>::Pointer pPartitioner
//                    = BoundingVolumePartitioner<ModelPart::ElementsContainerType>::Pointer(new SimpleBoundingVolumePartitioner<ModelPart::ElementsContainerType>());
//                BoundingVolumeTree<ModelPart::ElementsContainerType>::Pointer pTree
//                    = BoundingVolumeTree<ModelPart::ElementsContainerType>::Pointer(new BoundingVolumeTree<ModelPart::ElementsContainerType>(26));

//                pTree->BuildTreeTopDown(r_model_part.Elements(), *pPartitioner);

                // setup the tying links
                TyingLinkType SampleLink;
                ModelPart::ConditionsContainerType LinkingConditions;
//                Kratos::progress_display show_progress( rpNodes.size() );
                for( typename ModelPart::NodesContainerType::ptr_iterator it = rpNodes.ptr_begin();
                        it != rpNodes.ptr_end(); ++it )
                {
                    // define the list of potential elements
                    //TODO

                    Element::Pointer pTargetElement;
                    PointType LocalPoint;
                    if( FindPartnerElement( *(*it), rpElements, pTargetElement, LocalPoint ) )
                    {
                        // Create the linking condition
                        GeometryType::Pointer pTempGeometry = GeometryType::Pointer( new GeometryType() );
                        Condition::Pointer pNewLink = SampleLink.Create(++lastCondId, pTempGeometry, *it, pTargetElement, LocalPoint);
                        pNewLink->Set(ACTIVE, true);
                        pNewLink->SetValue(ASSOCIATED_NODE, *it);
                        pNewLink->SetValue(ASSOCIATED_ELEMENT, pTargetElement);
                        LinkingConditions.push_back( pNewLink );
                    }
//                    ++show_progress;
                }

                for( ModelPart::ConditionsContainerType::ptr_iterator it = LinkingConditions.ptr_begin();
                        it != LinkingConditions.ptr_end(); ++it )
                {
                    r_model_part.Conditions().push_back(*it);
                }

                std::cout << "Setup node tying links completed, "
                    << LinkingConditions.size() << " linking conditions of " << typeid(SampleLink).name()
                    << " was added to model_part" << std::endl;

                return LinkingConditions;
            }

            /********************************************************************************
            INTERFACE FOR TYING USING 2-NODE ELEMENT
            ********************************************************************************/

            ModelPart::ConditionsContainerType SetUpTyingLinks3( ModelPart& r_model_part, boost::python::list list_trusses,
                    boost::python::list list_elements )
            {
                ModelPart::ElementsContainerType pTrusses;
                ModelPart::ElementsContainerType pElements;

                typedef boost::python::stl_input_iterator<int> iterator_type;
                BOOST_FOREACH(const iterator_type::value_type& id,
                              std::make_pair(iterator_type(list_trusses), // begin
                                iterator_type() ) ) // end
                {
                    Element::Pointer pElem = r_model_part.pGetElement(id);
                    pTrusses.push_back(pElem);
                    if(pElem->GetGeometry().GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Line3D2)
                        KRATOS_THROW_ERROR(std::runtime_error, "The truss element is invalid, its Id is", pElem->Id())
                }

                typedef boost::python::stl_input_iterator<int> iterator_type;
                BOOST_FOREACH(const iterator_type::value_type& id,
                              std::make_pair(iterator_type(list_elements), // begin
                                iterator_type() ) ) // end
                {
                    pElements.push_back(r_model_part.pGetElement(id));
                }

                return this->SetUpTyingLinks4(r_model_part, pTrusses, pElements);
            }

            /**
             * Initializes tying conditions
             */
            ModelPart::ConditionsContainerType SetUpTyingLinks4( ModelPart& r_model_part, ModelPart::ElementsContainerType& rpTrusses,
                    ModelPart::ElementsContainerType& rpElements )
            {

                std::cout << "Initializing EmbeddedNodeTyingUtility..." << std::endl;
            }


            /********************************************************************************
            UTILITY FUNCTION
            ********************************************************************************/
            ModelPart::ConditionsContainerType Combine(ModelPart::ConditionsContainerType& List1, ModelPart::ConditionsContainerType& List2)
            {
                ModelPart::ConditionsContainerType Conds;

                for(ModelPart::ConditionsContainerType::ptr_iterator it = List1.ptr_begin(); it != List1.ptr_end(); ++it)
                {
                    Conds.push_back(*it);
                }
                for(ModelPart::ConditionsContainerType::ptr_iterator it = List2.ptr_begin(); it != List2.ptr_end(); ++it)
                {
                    Conds.push_back(*it);
                }

                return Conds;
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

    };//class EmbeddedNodeTyingUtility
}  // namespace Kratos.
#endif // KRATOS_EMBEDDED_NODE_TYING_UTILITY_INCLUDED

