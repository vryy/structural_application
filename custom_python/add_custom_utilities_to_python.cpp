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
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-23 14:27:01 $
//   Revision:            $Revision: 1.20 $
//
//


#if !defined(KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED


// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/foreach.hpp>

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/define.h"
#include "custom_utilities/deactivation_utility.h"
#include "custom_utilities/variable_transfer_utility.h"
#include "custom_utilities/variable_projection_utility.h"
#include "custom_utilities/variable_advanced_transfer_utility.h"

#ifdef _OPENMP
#include "custom_utilities/parallel_variable_transfer_utility.h"
#endif

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/contact_utility.h"
#include "custom_utilities/volume_utility.h"
#include "custom_utilities/restart_utility.h"
#include "custom_utilities/node_snapping_utility.h"
#include "custom_elements/rigid_body_3D.h"
#include "custom_utilities/output_utility.h"
#include "custom_utilities/dof_utility.h"
#include "custom_utilities/smoothing_utility.h"
//#include "custom_utilities/tip_utility.h"
#include "custom_utilities/pile_utility.h"
#include "custom_utilities/foundation_utility.h"

//#include "custom_utilities/detect_elements_utility.h"
#include "custom_utilities/intra_fracture_triangle_utility.h"
#include "custom_utilities/inter_fracture_triangle_utility.h"
#include "custom_utilities/inter_fracture_tetrahedra_utility.h"
//#include "custom_utilities/mark_element_for_refinement.h"
#include "custom_utilities/disconnect_utility.h"

#include "custom_utilities/embedded_node_tying_utility.h"
#include "custom_conditions/embedded_node_lagrange_tying_condition.h"
#include "custom_conditions/embedded_node_penalty_tying_condition.h"
#include "custom_utilities/embedded_point_tying_utility.h"
#include "custom_conditions/embedded_point_lagrange_tying_condition.h"
#include "custom_conditions/embedded_point_penalty_tying_condition.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void AddNewRigidBody3D( ModelPart& structural_model_part,
                        ModelPart& skin_model_part,
                        Variable<double>& rSelectionVariable,
                        double selection_value,
                        Node<3>::Pointer CenterNode,
                        Element::PropertiesType::Pointer pProperties,
                        double nodal_mass,
                        Matrix& Inertia
                      )
{
    Geometry<Node<3> >::Pointer skin_nodes_geometry( new Geometry<Node<3> > ); ;

    //selecting the nodes in the model part having rSelectionVariable==selection_value

    for ( ModelPart::NodesContainerType::iterator it = skin_model_part.NodesBegin(); it != skin_model_part.NodesEnd(); it++ )
    {
        if ( it->FastGetSolutionStepValue( rSelectionVariable ) == selection_value )
            skin_nodes_geometry->push_back( *( it.base() ) );
    }

    //creating a geometry containing the center node
    Geometry<Node<3> >::Pointer center_node_geometry( new Geometry<Node<3> > ) ;

    center_node_geometry->push_back( Node<3>::Pointer( CenterNode ) );

    unsigned int last_id = 1;

    if ( structural_model_part.Elements().size() != 0 )
        last_id = ( structural_model_part.ElementsEnd() - 1 )->Id() + 1;

    array_1d<double, 3> zero = ZeroVector( 3 );

    Element::Pointer new_el = RigidBody3D::Pointer( new  RigidBody3D( last_id,
                              center_node_geometry,
                              pProperties,
                              skin_nodes_geometry,
                              nodal_mass,
                              Inertia, zero, zero ) );

    structural_model_part.Elements().push_back(
        new_el
    );
}

void AddNewRigidBodyAndSpring3D( ModelPart& structural_model_part,
                                 ModelPart& skin_model_part,
                                 Variable<double>& rSelectionVariable,
                                 double selection_value,
                                 Node<3>::Pointer CenterNode,
                                 Element::PropertiesType::Pointer pProperties,
                                 double nodal_mass,
                                 Matrix& Inertia,
                                 array_1d<double, 3>& translational_stiffness,
                                 array_1d<double, 3>& rotational_stiffness
                               )
{
    Geometry<Node<3> >::Pointer skin_nodes_geometry( new Geometry<Node<3> > ); ;

    //selecting the nodes in the model part having rSelectionVariable==selection_value

    for ( ModelPart::NodesContainerType::iterator it = skin_model_part.NodesBegin(); it != skin_model_part.NodesEnd(); it++ )
    {
        if ( it->FastGetSolutionStepValue( rSelectionVariable ) == selection_value )
            skin_nodes_geometry->push_back( *( it.base() ) );
    }

    //creating a geometry containing the center node
    Geometry<Node<3> >::Pointer center_node_geometry( new Geometry<Node<3> > ) ;

    center_node_geometry->push_back( Node<3>::Pointer( CenterNode ) );

    unsigned int last_id = 1;

    if ( structural_model_part.Elements().size() != 0 )
        last_id = ( structural_model_part.ElementsEnd() - 1 )->Id() + 1;

    array_1d<double, 3> zero = ZeroVector( 3 );

    Element::Pointer new_el = RigidBody3D::Pointer( new  RigidBody3D( last_id,
                              center_node_geometry,
                              pProperties,
                              skin_nodes_geometry,
                              nodal_mass,
                              Inertia,
                              translational_stiffness,
                              rotational_stiffness ) );

    structural_model_part.Elements().push_back(
        new_el
    );
}

void DoubleTransferVariablesToNodes(VariableTransferUtility& dummy,
        ModelPart& model_part, Variable<double>& rThisVariable)
{
    dummy.TransferVariablesToNodes(model_part, rThisVariable);
}

void Array1DTransferVariablesToNodes(VariableTransferUtility& dummy,
        ModelPart& model_part, Variable<array_1d<double, 3> >& rThisVariable)
{
    dummy.TransferVariablesToNodes(model_part, rThisVariable);
}

void Array1DTransferVariablesToNodesForElements(VariableTransferUtility& dummy,
        ModelPart& model_part, ModelPart::ElementsContainerType& rElements, Variable<array_1d<double, 3> >& rThisVariable)
{
    dummy.TransferVariablesToNodes(model_part, rElements, rThisVariable);
}

void Array1DTransferVariablesToNodesForElementsAsList(VariableTransferUtility& dummy,
        ModelPart& model_part, boost::python::list& listElements, Variable<array_1d<double, 3> >& rThisVariable)
{
    ModelPart::ElementsContainerType rElements;
    typedef boost::python::stl_input_iterator<Element::Pointer> iterator_type;
    BOOST_FOREACH(const iterator_type::value_type& v,
                  std::make_pair(iterator_type(listElements), // begin
                    iterator_type() ) ) // end
        rElements.push_back(v);
    dummy.TransferVariablesToNodes(model_part, rElements, rThisVariable);
}

void Array1DTransferVariablesToNodesForConditions(VariableTransferUtility& dummy,
        ModelPart& model_part, ModelPart::ConditionsContainerType& rConditions, Variable<array_1d<double, 3> >& rThisVariable)
{
    dummy.TransferVariablesToNodes(model_part, rConditions, rThisVariable);
}

void Array1DTransferVariablesToNodesForConditionsAsList(VariableTransferUtility& dummy,
        ModelPart& model_part, boost::python::list& listConditions, Variable<array_1d<double, 3> >& rThisVariable)
{
    ModelPart::ConditionsContainerType rConditions;
    typedef boost::python::stl_input_iterator<Condition::Pointer> iterator_type;
    BOOST_FOREACH(const iterator_type::value_type& v,
                  std::make_pair(iterator_type(listConditions), // begin
                    iterator_type() ) ) // end
        rConditions.push_back(v);
    dummy.TransferVariablesToNodes(model_part, rConditions, rThisVariable);
}

void VectorTransferVariablesToNodes(VariableTransferUtility& dummy,
        ModelPart& model_part, Variable<Vector>& rThisVariable)
{
    dummy.TransferVariablesToNodes(model_part, rThisVariable);
}

void VectorTransferVariablesToNodesComponents(VariableTransferUtility& dummy,
        ModelPart& model_part, Variable<Vector>& rThisVariable, const std::size_t& ncomponents)
{
    dummy.TransferVariablesToNodes(model_part, rThisVariable, ncomponents);
}

void VectorTransferVariablesToNodesComponentsForElements(VariableTransferUtility& dummy,
        ModelPart& model_part, ModelPart::ElementsContainerType& rElements, Variable<Vector>& rThisVariable, const std::size_t& ncomponents)
{
    dummy.TransferVariablesToNodes(model_part, rElements, rThisVariable, ncomponents);
}

void VectorTransferVariablesToNodesComponentsForElementsAsList(VariableTransferUtility& dummy,
        ModelPart& model_part, boost::python::list& listElements, Variable<Vector>& rThisVariable, const std::size_t& ncomponents)
{
    ModelPart::ElementsContainerType rElements;
    typedef boost::python::stl_input_iterator<Element::Pointer> iterator_type;
    BOOST_FOREACH(const iterator_type::value_type& v,
                  std::make_pair(iterator_type(listElements), // begin
                    iterator_type() ) ) // end
        rElements.push_back(v);
    dummy.TransferVariablesToNodes(model_part, rElements, rThisVariable, ncomponents);
}

void VectorTransferVariablesToNodesComponentsForConditions(VariableTransferUtility& dummy,
        ModelPart& model_part, ModelPart::ConditionsContainerType& rConditions, Variable<Vector>& rThisVariable, const std::size_t& ncomponents)
{
    dummy.TransferVariablesToNodes(model_part, rConditions, rThisVariable, ncomponents);
}

void VectorTransferVariablesToNodesComponentsForConditionsAsList(VariableTransferUtility& dummy,
        ModelPart& model_part, boost::python::list& listConditions, Variable<Vector>& rThisVariable, const std::size_t& ncomponents)
{
    ModelPart::ConditionsContainerType rConditions;
    typedef boost::python::stl_input_iterator<Condition::Pointer> iterator_type;
    BOOST_FOREACH(const iterator_type::value_type& v,
                  std::make_pair(iterator_type(listConditions), // begin
                    iterator_type() ) ) // end
        rConditions.push_back(v);
    dummy.TransferVariablesToNodes(model_part, rConditions, rThisVariable, ncomponents);
}

void DoubleTransferVariablesToGaussPoints(VariableTransferUtility& dummy,
        ModelPart& source_model_part, ModelPart& target_model_part, Variable<double>& rThisVariable)
{
    dummy.TransferVariablesToGaussPoints(source_model_part, target_model_part, rThisVariable);
}

void VectorTransferVariablesToGaussPoints(VariableTransferUtility& dummy,
        ModelPart& source_model_part, ModelPart& target_model_part, Variable<Vector>& rThisVariable)
{
    dummy.TransferVariablesToGaussPoints(source_model_part, target_model_part, rThisVariable);
}

void VectorTransferVariablesToGaussPointsElementComponents(VariableTransferUtility& dummy,
        ModelPart& source_model_part, Element::Pointer pTargetElement, Variable<Vector>& rThisVariable, std::size_t ncomponents)
{
    dummy.TransferVariablesToGaussPoints(source_model_part, pTargetElement, rThisVariable, ncomponents);
}

void VectorTransferVariablesToGaussPointsElementComponentsIdentically(VariableTransferUtility& dummy,
        ModelPart& source_model_part, Element::Pointer pTargetElement, Variable<Vector>& rThisVariable, std::size_t ncomponents)
{
    dummy.TransferVariablesToGaussPointsIdentically(source_model_part, pTargetElement, rThisVariable, ncomponents);
}

void DoubleTransferVariablesBetweenMeshes(VariableTransferUtility& dummy,
        ModelPart& rSource, ModelPart& rTarget, Variable<double>& rThisVariable)
{
    dummy.TransferVariablesBetweenMeshes(rSource, rTarget, rThisVariable);
}

void VectorTransferVariablesBetweenMeshes(VariableTransferUtility& dummy,
        ModelPart& rSource, ModelPart& rTarget, Variable<Kratos::Vector>& rThisVariable)
{
    dummy.TransferVariablesBetweenMeshes(rSource, rTarget, rThisVariable);
}

void ListDofs(DofUtility& dummy, ModelPart::DofsArrayType& rDofSet, std::size_t EquationSystemSize)
{
    dummy.ListDofs(rDofSet, EquationSystemSize);
}

template<class TVariableType>
void PrintKey(DofUtility& dummy, const TVariableType& rThisVariable)
{
    dummy.PrintKey(rThisVariable);
}

void SetAssociatedElement(DeactivationUtility& rDummy, Condition::Pointer pCond, Element::Pointer pElem)
{
    pCond->SetValue(ASSOCIATED_ELEMENT, pElem);
}

///////////////////////////////////////////////////////////////////////
// Auxilliary Utilities for connection of the building to the ground //
///////////////////////////////////////////////////////////////////////
void InitializePileUtility( PileUtility& dummy, ModelPart& model_part,
                            boost::python::list pile_elements, int len_pile_elements,
                            boost::python::list soil_elements, int len_soil_elements )
{
    std::vector<unsigned int> vec_pile_elements;
    std::vector<unsigned int> vec_soil_elements;

    for ( int it = 0; it < len_pile_elements; it++ )
    {
       boost::python::extract<int> x( pile_elements[it] );
 
       if ( x.check() )
           vec_pile_elements.push_back(( unsigned int )x );
       else break;
    }
 
    for ( int it = 0; it < len_soil_elements; it++ )
    {
        boost::python::extract<int> x( soil_elements[it] );

        if ( x.check() )
           vec_soil_elements.push_back(( unsigned int )x );
        else break;
    }
 
    dummy.InitializePileUtility( model_part, vec_pile_elements, vec_soil_elements );
}

///////////////////////////////////////////////////////////////////////
// Auxilliary Utilities for connection of the building to the ground //
///////////////////////////////////////////////////////////////////////
void InitializeFoundationUtility( FoundationUtility& dummy, ModelPart& model_part,
                            boost::python::list foundation_elements, int len_foundation_elements,
                            boost::python::list soil_elements, int len_soil_elements )
{
    std::vector<unsigned int> vec_foundation_elements;
    std::vector<unsigned int> vec_soil_elements;

    for ( int it = 0; it < len_foundation_elements; it++ )
    {
       boost::python::extract<int> x( foundation_elements[it] );
 
       if ( x.check() )
           vec_foundation_elements.push_back(( unsigned int )x );
       else break;
    }
 
    for ( int it = 0; it < len_soil_elements; it++ )
    {
        boost::python::extract<int> x( soil_elements[it] );

        if ( x.check() )
           vec_soil_elements.push_back(( unsigned int )x );
        else break;
    }
 
    dummy.InitializeFoundationUtility( model_part, vec_foundation_elements, vec_soil_elements );
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

void  AddCustomUtilitiesToPython()
{
    class_<DeactivationUtility, boost::noncopyable >
    ( "DeactivationUtility", init<>() )
    .def( init<int>() )
    .def( "Deactivate", &DeactivationUtility::Deactivate )
    .def( "Reactivate", &DeactivationUtility::Reactivate )
    .def( "ReactivateStressFree", &DeactivationUtility::ReactivateStressFree )
    .def( "ReactivateAll", &DeactivationUtility::ReactivateAll )
    .def( "Initialize", &DeactivationUtility::Initialize )
    .def( "GetName", &DeactivationUtility::GetName<Element> )
    .def( "GetName", &DeactivationUtility::GetName<Condition> )
    .def( "SetAssociatedElement", &SetAssociatedElement )
    ;

    class_<VariableTransferUtility, boost::noncopyable >
    ( "VariableTransferUtility", init<>() )
    .def(init<VariableTransferUtility::LinearSolverType::Pointer>())
    .def( "TransferNodalVariables", &VariableTransferUtility::TransferNodalVariables )
    .def( "TransferNodalVariables", &VariableTransferUtility::TransferGeneralNodalVariables<Variable<double> > )
    .def( "TransferNodalVariables", &VariableTransferUtility::TransferGeneralNodalVariables<Variable<Vector> > )
    .def( "TransferConstitutiveLawVariables", &VariableTransferUtility::TransferConstitutiveLawVariables )
    .def( "TransferInSituStress", &VariableTransferUtility::TransferInSituStress )
    .def( "TransferPrestress", &VariableTransferUtility::TransferPrestress )
    .def( "TransferPrestressIdentically", &VariableTransferUtility::TransferPrestressIdentically )
    .def( "TransferSpecificVariable", &VariableTransferUtility::TransferSpecificVariable )
    .def( "TransferSpecificVariableWithComponents", &VariableTransferUtility::TransferSpecificVariableWithComponents )
    .def( "InitializeModelPart", &VariableTransferUtility::InitializeModelPart )
    .def("TransferVariablesToNodes", &DoubleTransferVariablesToNodes)
    .def("TransferVariablesToNodes", &Array1DTransferVariablesToNodes)
    .def("TransferVariablesToNodes", &Array1DTransferVariablesToNodesForElements)
    .def("TransferVariablesToNodesForElements", &Array1DTransferVariablesToNodesForElementsAsList)
    .def("TransferVariablesToNodes", &Array1DTransferVariablesToNodesForConditions)
    .def("TransferVariablesToNodesForConditions", &Array1DTransferVariablesToNodesForConditionsAsList)
    .def("TransferVariablesToNodes", &VectorTransferVariablesToNodes)
    .def("TransferVariablesToNodes", &VectorTransferVariablesToNodesComponents)
    .def("TransferVariablesToNodes", &VectorTransferVariablesToNodesComponentsForElements)
    .def("TransferVariablesToNodesForElements", &VectorTransferVariablesToNodesComponentsForElementsAsList)
    .def("TransferVariablesToNodes", &VectorTransferVariablesToNodesComponentsForConditions)
    .def("TransferVariablesToNodesForConditions", &VectorTransferVariablesToNodesComponentsForConditionsAsList)
    .def("TransferVariablesToGaussPoints", &DoubleTransferVariablesToGaussPoints)
    .def("TransferVariablesToGaussPoints", &VectorTransferVariablesToGaussPoints)
    .def("TransferVariablesToGaussPoints", &VectorTransferVariablesToGaussPointsElementComponents)
    .def("TransferVariablesToGaussPointsIdentically", &VectorTransferVariablesToGaussPointsElementComponentsIdentically)
    .def("TransferVariablesFromNodeToNode", &VariableTransferUtility::TransferVariablesFromNodeToNode<Variable<Vector> >)
    .def("TransferVariablesFromNodeToNode", &VariableTransferUtility::TransferVariablesFromNodeToNode<Variable<array_1d<double, 3> > >)
    .def("TransferVariablesFromNodeToNode", &VariableTransferUtility::TransferVariablesFromNodeToNode<Variable<double> >)
    .def("TransferVariablesBetweenMeshes", &DoubleTransferVariablesBetweenMeshes)
    .def("TransferVariablesBetweenMeshes", &VectorTransferVariablesBetweenMeshes)
    ;

    class_<VariableProjectionUtility, boost::noncopyable >
    ( "VariableProjectionUtility", init<VariableProjectionUtility::LinearSolverType::Pointer>() )
    .def("SetEchoLevel", &VariableProjectionUtility::SetEchoLevel)
    .def("BeginProjection", &VariableProjectionUtility::BeginProjection)
    .def("TransferVariablesToNodes", &VariableProjectionUtility::TransferVariablesToNodes)
    .def("EndProjection", &VariableProjectionUtility::EndProjection)
    ;

    class_<VariableAdvancedTransferUtility, boost::noncopyable >
    ( "VariableAdvancedTransferUtility", init<ModelPart::ElementsContainerType&, const double&, const double&, const double&>() )
    .def("SetEchoLevel", &VariableAdvancedTransferUtility::SetEchoLevel)
    .def("TransferVariablesFromNodeToNode", &VariableAdvancedTransferUtility::TransferVariablesFromNodeToNode<Variable<double> >)
    .def("TransferVariablesFromNodeToNode", &VariableAdvancedTransferUtility::TransferVariablesFromNodeToNode<Variable<array_1d<double, 3> > >)
    .def("TransferVariablesFromNodeToNode", &VariableAdvancedTransferUtility::TransferVariablesFromNodeToNode<Variable<Vector> >)
    .def("TransferVariablesFromNodeToNode", &VariableAdvancedTransferUtility::TransferVariablesFromNodeToNode<Variable<Matrix> >)
    ;

#ifdef _OPENMP
    class_<ParallelVariableTransferUtility, boost::noncopyable >
    ( "ParallelVariableTransferUtility",
      init<>() )
    .def( "TransferNodalVariables", &ParallelVariableTransferUtility::TransferNodalVariables )
    .def( "TransferConstitutiveLawVariables", &ParallelVariableTransferUtility::TransferConstitutiveLawVariables )
    .def( "TransferInSituStress", &ParallelVariableTransferUtility::TransferInSituStress )
    .def( "InitializeModelPart", &ParallelVariableTransferUtility::InitializeModelPart )
    ;
#endif

    class_<ContactUtility, boost::noncopyable >
    ( "ContactUtility",
      init<int>() )
    .def( "SetUpContactConditions", &ContactUtility::SetUpContactConditions )
    .def( "SetUpContactConditionsLagrangeTying", &ContactUtility::SetUpContactConditionsLagrangeTying )
    .def( "Update", &ContactUtility::Update )
    .def( "IsConverged", &ContactUtility::IsConverged )
    .def( "Clean", &ContactUtility::Clean )
    .def( "CleanLagrangeTying", &ContactUtility::CleanLagrangeTying )
    ;
// VM
    class_<VolumeUtility, boost::noncopyable >
    ( "VolumeUtility",
      init<int>() )
    .def( "Calculate_this_Volume", &VolumeUtility::CalculateVolume ) // VM
    ;
//VM

    class_<RestartUtility, boost::noncopyable >
    ( "RestartUtility",
      init< std::string const& >() )
    .def( "ChangeFileName", &RestartUtility::ChangeFileName )
    .def( "StoreNodalVariables", &RestartUtility::StoreNodalVariables )
    .def( "WriteNodalVariables", &RestartUtility::WriteNodalVariables )
    .def( "StoreInSituStress", &RestartUtility::StoreInSituStress )
    .def( "WriteConstitutiveLawVariables", &RestartUtility::WriteConstitutiveLawVariables )
    .def( "StoreConstitutiveLawVariables", &RestartUtility::StoreConstitutiveLawVariables )
    .def( "WriteInSituStress", &RestartUtility::WriteInSituStress )
    ;

    class_<NodeSnappingUtility, boost::noncopyable >
    ( "NodeSnappingUtility",
      init<>() )
    .def( "MoveNode", &NodeSnappingUtility::MoveNode )
    .def( "AdjustNodes", &NodeSnappingUtility::AdjustNodes )
    .def( "AdjustToCircle", &NodeSnappingUtility::AdjustToCircle )
    .def( "AdjustToCylinder", &NodeSnappingUtility::AdjustToCylinder )
    .def( "AdjustToClosedCylinder", &NodeSnappingUtility::AdjustToClosedCylinder )
    .def( "IdentifyInsideElements", &NodeSnappingUtility::IdentifyInsideElements )
    .def( "SetInsituStress", &NodeSnappingUtility::SetInsituStress )
    .def( "ExtractCapNodes", &NodeSnappingUtility::ExtractCapNodes )
    .def( "TestElements", &NodeSnappingUtility::TestElements )
    ;

    class_<OutputUtility, boost::noncopyable >
    ( "OutputUtility",
      init<>() )
    .def( "GetStrain", &OutputUtility::GetStrain )
    .def( "GetStress", &OutputUtility::GetStress )
    .def( "GetInternalVariables", &OutputUtility::GetInternalVariables )
    ;


    def( "AddNewRigidBody3D", AddNewRigidBody3D );
    def( "AddNewRigidBodyAndSpring3D", AddNewRigidBodyAndSpring3D );
    ;

    /*
                class_<Detect_Elements_And_Nodes, boost::noncopyable >
                        ("DetectElementsAndNodes", init<ModelPart&, int >() )
          .def("DetectNode",              &Detect_Elements_And_Nodes::Detect_Node_To_Be_Splitted)
                        .def("DetectElements",          &Detect_Elements_And_Nodes::Detect_Elements_To_Be_Splitted)
                        .def("CalculateMapFailure",     &Detect_Elements_And_Nodes::Calculate_Map_Failure)
                        .def("Finalize",                &Detect_Elements_And_Nodes::Finalize)
                        ;
    */
    class_<Smoothing_Utility, boost::noncopyable >
    ( "SmoothingUtility", init<ModelPart&, int >() )
    .def( "WeightedRecoveryGradients", &Smoothing_Utility::WeightedRecoveryGradients<double> )
    .def( "WeightedRecoveryGradients", &Smoothing_Utility::WeightedRecoveryGradients<Matrix> ) // for matrices
    .def( "InterpolatedRecoveryGradients", &Smoothing_Utility::InterpolatedRecoveryGradients<Matrix> )
    .def( "SettingNodalValues", &Smoothing_Utility::SettingNodalValues )
    .def( "RecomputeValuesForNewMesh", &Smoothing_Utility::Recompute_Values_For_New_Mesh )
    .def( "Finalize", &Smoothing_Utility::Finalize )
    .def( "SettingNodalValues", &Smoothing_Utility::SettingNodalValues )
    ;




    class_<Disconnect_Triangle_Utilities, boost::noncopyable >
    ( "DisconnectTriangle", init<ModelPart&>() )
    .def( "DisconnectElements", &Disconnect_Triangle_Utilities::Disconnect_Elements )
    ;


    class_<Intra_Fracture_Triangle, boost::noncopyable >
    ( "IntraFractureTriangle", init<ModelPart&, int >() )
    .def( "DetectAndSplitElements",              &Intra_Fracture_Triangle::Detect_And_Split_Elements )
    ;

    class_<Inter_Fracture_Triangle, boost::noncopyable >
    ( "InterFractureTriangle", init<ModelPart&, int >() )
    .def( "DetectAndSplitElementsHeuristicFormula", &Inter_Fracture_Triangle::Detect_And_Split_Elements_Heuristic_Formula )
    .def( "DetectAndSplitElements",              &Inter_Fracture_Triangle::Detect_And_Split_Elements )
    .def( "Finalize",                            &Inter_Fracture_Triangle::Finalize )
    ;

    class_<Inter_Fracture_Tetrahedra, boost::noncopyable >
    ( "InterFractureTetrahedra", init<ModelPart&, int >() )
    .def( "DetectAndSplitElements",              &Inter_Fracture_Tetrahedra::Detect_And_Split_Elements )
    ;

    class_<DofUtility, boost::noncopyable >
    ( "DofUtility", init<>() )
    .def( "ListDofs", &ListDofs )
    .def( "PrintKey", &PrintKey<Variable<double> > )
    .def( "PrintKey", &PrintKey<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > > )
    ;

    typedef EmbeddedNodeTyingUtility<EmbeddedNodeLagrangeTyingCondition> EmbeddedNodeLagrangeTyingUtilityType;
    class_<EmbeddedNodeLagrangeTyingUtilityType, EmbeddedNodeLagrangeTyingUtilityType::Pointer, boost::noncopyable >
    ( "EmbeddedNodeLagrangeTyingUtility", init<>() )
    .def( "SetUpTyingLinks", &EmbeddedNodeLagrangeTyingUtilityType::SetUpTyingLinks1 )
    .def( "SetUpTyingLinks", &EmbeddedNodeLagrangeTyingUtilityType::SetUpTyingLinks2 )
    ;

    typedef EmbeddedNodeTyingUtility<EmbeddedNodePenaltyTyingCondition> EmbeddedNodePenaltyTyingUtilityType;
    class_<EmbeddedNodePenaltyTyingUtilityType, EmbeddedNodePenaltyTyingUtilityType::Pointer, boost::noncopyable >
    ( "EmbeddedNodePenaltyTyingUtility", init<>() )
    .def( "SetUpTyingLinks", &EmbeddedNodePenaltyTyingUtilityType::SetUpTyingLinks1 )
    .def( "SetUpTyingLinks", &EmbeddedNodePenaltyTyingUtilityType::SetUpTyingLinks2 )
    .def( "Combine", &EmbeddedNodePenaltyTyingUtilityType::Combine )
    ;

//    typedef EmbeddedNodeTyingUtility<YourCondition> EmbeddedNodeFrictionalTrussTyingUtilityType;
//    class_<EmbeddedNodeFrictionalTrussTyingUtilityType, EmbeddedNodeFrictionalTrussTyingUtilityType::Pointer, boost::noncopyable >
//    ( "EmbeddedNodeFrictionalTrussTyingUtility", init<>() )
//    .def( "SetUpTyingLinks", &EmbeddedNodeFrictionalTrussTyingUtilityType::SetUpTyingLinks3 )
//    .def( "SetUpTyingLinks", &EmbeddedNodeFrictionalTrussTyingUtilityType::SetUpTyingLinks4 )
//    ;

    typedef EmbeddedPointTyingUtility<EmbeddedPointLagrangeTyingCondition> EmbeddedPointLagrangeTyingUtilityType;
    class_<EmbeddedPointLagrangeTyingUtilityType, EmbeddedPointLagrangeTyingUtilityType::Pointer, boost::noncopyable >
    ( "EmbeddedPointLagrangeTyingUtility", init<>() )
    .def( "SetUpTyingLinks", &EmbeddedPointLagrangeTyingUtilityType::SetUpTyingLinks1 )
    .def( "SetUpTyingLinks", &EmbeddedPointLagrangeTyingUtilityType::SetUpTyingLinks2 )
    ;

    typedef EmbeddedPointTyingUtility<EmbeddedPointPenaltyTyingCondition> EmbeddedPointPenaltyTyingUtilityType;
    class_<EmbeddedPointPenaltyTyingUtilityType, EmbeddedPointPenaltyTyingUtilityType::Pointer, boost::noncopyable >
    ( "EmbeddedPointPenaltyTyingUtility", init<>() )
    .def( "SetUpTyingLinks", &EmbeddedPointPenaltyTyingUtilityType::SetUpTyingLinks1 )
    .def( "SetUpTyingLinks", &EmbeddedPointPenaltyTyingUtilityType::SetUpTyingLinks2 )
    ;

    class_<PileUtility, boost::noncopyable >
    ( "PileUtility", init<>() )
    .def( "InitializePileUtility", &InitializePileUtility )
    ;

    class_<FoundationUtility, boost::noncopyable >
    ( "FoundationUtility", init<>() )
    .def( "InitializeFoundationUtility", &InitializeFoundationUtility )
    ;

}
}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED  defined 
