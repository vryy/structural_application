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

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/foreach.hpp>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/deactivation_utility.h"
#include "custom_utilities/variable_transfer_utility.h"
#include "custom_utilities/variable_utility.h"
#include "custom_utilities/variable_projection_utility.h"
#include "custom_utilities/variable_interpolation_utility.h"
#include "custom_utilities/variable_binning_interpolation_utility.h"
#include "custom_utilities/variable_bvh_interpolation_utility.h"

#ifdef _OPENMP
#include "custom_utilities/parallel_variable_transfer_utility.h"
#endif

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/contact_utility.h"
#include "custom_utilities/restart_utility.h"
#include "custom_elements/rigid_body_3D.h"
#include "custom_utilities/output_utility.h"
#include "custom_utilities/dof_utility.h"
#include "custom_utilities/smoothing_utility.h"
#include "custom_utilities/tip_utility.h"
#include "custom_utilities/pile_utility.h"
#include "custom_utilities/foundation_utility.h"

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

void DoubleTransferVariablesToNodesForElements(VariableTransferUtility& dummy,
        ModelPart& model_part, ModelPart::ElementsContainerType& rElements, Variable<double>& rThisVariable)
{
    dummy.TransferVariablesToNodes(model_part, rElements, rThisVariable);
}

void DoubleTransferVariablesToNodesForElementsAsList(VariableTransferUtility& dummy,
        ModelPart& model_part, boost::python::list& listElements, Variable<double>& rThisVariable)
{
    ModelPart::ElementsContainerType rElements;
    typedef boost::python::stl_input_iterator<Element::Pointer> iterator_type;
    BOOST_FOREACH(const iterator_type::value_type& v,
                  std::make_pair(iterator_type(listElements), // begin
                    iterator_type() ) ) // end
        rElements.push_back(v);
    dummy.TransferVariablesToNodes(model_part, rElements, rThisVariable);
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

void DoubleTransferVariablesToNode(VariableTransferUtility& dummy,
        ElementsContainerType& rElements, Variable<double>& rThisVariable, Element::GeometryType::PointType& rNode)
{
    dummy.TransferVariablesToNode(rElements, rThisVariable, rNode);
}

void Array1DTransferVariablesToNode(VariableTransferUtility& dummy,
        ElementsContainerType& rElements, Variable<array_1d<double, 3> >& rThisVariable, Element::GeometryType::PointType& rNode)
{
    dummy.TransferVariablesToNode(rElements, rThisVariable, rNode);
}

boost::python::list DoubleComputeExtrapolatedNodalValues(VariableTransferUtility& dummy,
        Element& source_element, Variable<double>& rThisVariable,
        const ProcessInfo& CurrentProcessInfo)
{
    std::vector<double> values;
    dummy.ComputeExtrapolatedNodalValues(values, source_element, rThisVariable, CurrentProcessInfo);
    boost::python::list output;
    for (std::size_t i = 0; i < values.size(); ++i)
        output.append(values[i]);
    return output;
}

boost::python::list Array1DComputeExtrapolatedNodalValues(VariableTransferUtility& dummy,
        Element& source_element, Variable<array_1d<double, 3> >& rThisVariable,
        const ProcessInfo& CurrentProcessInfo)
{
    std::vector<array_1d<double, 3> > values;
    dummy.ComputeExtrapolatedNodalValues(values, source_element, rThisVariable, CurrentProcessInfo);
    boost::python::list output;
    for (std::size_t i = 0; i < values.size(); ++i)
        output.append(values[i]);
    return output;
}

boost::python::list VectorComputeExtrapolatedNodalValues(VariableTransferUtility& dummy,
        Element& source_element, Variable<Vector>& rThisVariable,
        const ProcessInfo& CurrentProcessInfo, const std::size_t& ncomponents)
{
    std::vector<Vector> values;
    dummy.ComputeExtrapolatedNodalValues(values, source_element, rThisVariable, CurrentProcessInfo, ncomponents);
    boost::python::list output;
    for (std::size_t i = 0; i < values.size(); ++i)
        output.append(values[i]);
    return output;
}

void DoubleTransferVariablesToGaussPointsFromNodalValues(VariableTransferUtility& dummy,
        boost::python::list& list_values, Element& target_element, Variable<double>& rThisVariable,
        const ProcessInfo& CurrentProcessInfo)
{
    std::vector<double> values;
    typedef boost::python::stl_input_iterator<double> iterator_type;
    BOOST_FOREACH(const iterator_type::value_type& v,
            std::make_pair(iterator_type(list_values), iterator_type() ) )
        values.push_back(v);
    dummy.TransferVariablesToGaussPoints(values, target_element, rThisVariable, CurrentProcessInfo);
}

void Array1DTransferVariablesToGaussPointsFromNodalValues(VariableTransferUtility& dummy,
        boost::python::list& list_values, Element& target_element, Variable<array_1d<double, 3> >& rThisVariable,
        const ProcessInfo& CurrentProcessInfo)
{
    std::vector<array_1d<double, 3> > values;
    typedef boost::python::stl_input_iterator<array_1d<double, 3> > iterator_type;
    BOOST_FOREACH(const iterator_type::value_type& v,
            std::make_pair(iterator_type(list_values), iterator_type() ) )
        values.push_back(v);
    dummy.TransferVariablesToGaussPoints(values, target_element, rThisVariable, CurrentProcessInfo);
}

void VectorTransferVariablesToGaussPointsFromNodalValues(VariableTransferUtility& dummy,
        boost::python::list& list_values, Element& target_element, Variable<Vector>& rThisVariable,
        const ProcessInfo& CurrentProcessInfo, const std::size_t& ncomponents)
{
    std::vector<Vector> values;
    typedef boost::python::stl_input_iterator<Vector> iterator_type;
    BOOST_FOREACH(const iterator_type::value_type& v,
            std::make_pair(iterator_type(list_values), iterator_type() ) )
        values.push_back(v);
    dummy.TransferVariablesToGaussPoints(values, target_element, rThisVariable, CurrentProcessInfo, ncomponents);
}

void DoubleTransferVariablesToGaussPoints(VariableTransferUtility& dummy,
        ModelPart& source_model_part, ModelPart& target_model_part, Variable<double>& rThisVariable)
{
    dummy.TransferVariablesToGaussPoints(source_model_part, target_model_part, rThisVariable);
}

void DoubleTransferVariablesToGaussPointsLocal(VariableTransferUtility& dummy,
        Element& source_element, Element& target_element, Variable<double>& rThisVariable,
        const ProcessInfo& CurrentProcessInfo)
{
    dummy.TransferVariablesToGaussPoints(source_element, target_element, rThisVariable, CurrentProcessInfo);
}

void Array1DTransferVariablesToGaussPointsElements(VariableTransferUtility& dummy,
        ModelPart::ElementsContainerType& source_elements,
        ModelPart::ElementsContainerType& target_elements,
        Variable<array_1d<double, 3> >& rThisVariable,
        const ProcessInfo& CurrentProcessInfo)
{
    dummy.TransferVariablesToGaussPoints(source_elements, target_elements, rThisVariable, CurrentProcessInfo);
}

void Array1DTransferVariablesToGaussPoints(VariableTransferUtility& dummy,
        ModelPart& source_model_part, ModelPart& target_model_part, Variable<array_1d<double, 3> >& rThisVariable)
{
    dummy.TransferVariablesToGaussPoints(source_model_part, target_model_part, rThisVariable);
}

void Array1DTransferVariablesToGaussPointsLocal(VariableTransferUtility& dummy,
        Element& source_element, Element& target_element, Variable<array_1d<double, 3> >& rThisVariable,
        const ProcessInfo& CurrentProcessInfo)
{
    dummy.TransferVariablesToGaussPoints(source_element, target_element, rThisVariable, CurrentProcessInfo);
}

void VectorTransferVariablesToGaussPoints(VariableTransferUtility& dummy,
        ModelPart& source_model_part, ModelPart& target_model_part, Variable<Vector>& rThisVariable)
{
    dummy.TransferVariablesToGaussPoints(source_model_part, target_model_part, rThisVariable);
}

void VectorTransferVariablesToGaussPointsFromElementsToElements(VariableTransferUtility& dummy,
        ModelPart::ElementsContainerType& source_elements, ModelPart::ElementsContainerType& target_elements,
        const ProcessInfo& CurrentProcessInfo, Variable<Vector>& rThisVariable)
{
    dummy.TransferVariablesToGaussPoints(source_elements, target_elements, CurrentProcessInfo, rThisVariable);
}

void VectorTransferVariablesToGaussPointsComponents(VariableTransferUtility& dummy,
        ModelPart& source_model_part, ModelPart& target_model_part, Variable<Vector>& rThisVariable, const std::size_t& ncomponents)
{
    dummy.TransferVariablesToGaussPoints(source_model_part, target_model_part, rThisVariable, ncomponents);
}

void VectorTransferVariablesToGaussPointsFromElementsToElementsComponents(VariableTransferUtility& dummy,
        ModelPart::ElementsContainerType& source_elements, ModelPart::ElementsContainerType& target_elements,
        const ProcessInfo& CurrentProcessInfo, Variable<Vector>& rThisVariable, const std::size_t& ncomponents)
{
    dummy.TransferVariablesToGaussPoints(source_elements, target_elements, CurrentProcessInfo, rThisVariable, ncomponents);
}

void VectorTransferVariablesToGaussPointsLocal(VariableTransferUtility& dummy,
        Element& source_element, Element& target_element, Variable<Vector>& rThisVariable,
        const ProcessInfo& CurrentProcessInfo, const std::size_t& ncomponents)
{
    dummy.TransferVariablesToGaussPoints(source_element, target_element, rThisVariable, CurrentProcessInfo, ncomponents);
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

void DoubleTransferVariablesToGaussPointsElementComponentsIdentically(VariableTransferUtility& dummy,
        ModelPart::ElementsContainerType& source_elements, ModelPart::ElementsContainerType& target_elements,
        Variable<double>& rThisVariable, const ProcessInfo& CurrentProcessInfo)
{
    dummy.TransferVariablesToGaussPointsIdentically(source_elements, target_elements, rThisVariable, CurrentProcessInfo);
}

void VectorTransferVariablesToGaussPointsElementComponentsIdentically2(VariableTransferUtility& dummy,
        ModelPart::ElementsContainerType& source_elements, ModelPart::ElementsContainerType& target_elements,
        Variable<Vector>& rThisVariable, const ProcessInfo& CurrentProcessInfo)
{
    dummy.TransferVariablesToGaussPointsIdentically(source_elements, target_elements, rThisVariable, CurrentProcessInfo);
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
void RemoveDof(DofUtility& dummy, ModelPart::DofsArrayType& rDofSet, const TVariableType& rThisVariable)
{
    dummy.RemoveDof(rDofSet, rThisVariable);
}

template<class TVariableType>
void PrintKey(DofUtility& dummy, const TVariableType& rThisVariable)
{
    dummy.PrintKey(rThisVariable);
}

void PrintMaxUnbalancedForce(DofUtility& dummy, ModelPart::DofsArrayType& rDofSet, const Vector& forces)
{
    dummy.PrintMaxUnbalancedForce(rDofSet, forces);
}

void PrintMinUnbalancedForce(DofUtility& dummy, ModelPart::DofsArrayType& rDofSet, const Vector& forces)
{
    dummy.PrintMinUnbalancedForce(rDofSet, forces);
}

void SetAssociatedElement(DeactivationUtility& rDummy, Condition::Pointer pCond, Element::Pointer pElem)
{
    pCond->SetValue(ASSOCIATED_ELEMENT, pElem);
}

///////////////////////////////////////////////////////////////////////
// Auxilliary Utilities for connection of the pile to the ground //
///////////////////////////////////////////////////////////////////////
void InitializePileUtility( PileUtility& dummy, ModelPart& model_part,
                            boost::python::list pile_elements,
                            boost::python::list soil_elements,
                            Properties::Pointer linkProperties )
{
    std::vector<unsigned int> vec_pile_elements;
    std::vector<unsigned int> vec_soil_elements;

    int len_pile_elements = boost::python::len(pile_elements);
    for ( int it = 0; it < len_pile_elements; it++ )
    {
       boost::python::extract<int> x( pile_elements[it] );

       if ( x.check() )
           vec_pile_elements.push_back(( unsigned int )x );
       else break;
    }

    int len_soil_elements = boost::python::len(soil_elements);
    for ( int it = 0; it < len_soil_elements; it++ )
    {
        boost::python::extract<int> x( soil_elements[it] );

        if ( x.check() )
           vec_soil_elements.push_back(( unsigned int )x );
        else break;
    }

    dummy.InitializePileUtility( model_part, vec_pile_elements, vec_soil_elements, linkProperties );
}

///////////////////////////////////////////////////////////////////////
// Auxilliary Utilities for connection of the pile to the ground //
///////////////////////////////////////////////////////////////////////
void InitializeTipUtility( TipUtility& dummy, ModelPart& model_part,
                           boost::python::list tip_elements,
                           boost::python::list tip_soil_elements,
                           Properties::Pointer linkProperties )
{
    std::vector<unsigned int> vec_tip_elements;
    std::vector<unsigned int> vec_tip_soil_elements;

    int len_tip_elements = boost::python::len(tip_elements);
    for ( int it = 0; it < len_tip_elements; it++ )
    {
       boost::python::extract<int> x( tip_elements[it] );

       if ( x.check() )
           vec_tip_elements.push_back(( unsigned int )x );
       else break;
    }

    int len_tip_soil_elements = boost::python::len(tip_soil_elements);
    for ( int it = 0; it < len_tip_soil_elements; it++ )
    {
        boost::python::extract<int> x( tip_soil_elements[it] );

        if ( x.check() )
           vec_tip_soil_elements.push_back(( unsigned int )x );
        else break;
    }

    dummy.InitializeTipUtility( model_part, vec_tip_elements, vec_tip_soil_elements, linkProperties );
}

///////////////////////////////////////////////////////////////////////
// Auxilliary Utilities for connection of the building to the ground //
///////////////////////////////////////////////////////////////////////
void InitializeFoundationUtility( FoundationUtility& dummy, ModelPart& model_part,
                            boost::python::list foundation_elements,
                            boost::python::list soil_elements,
                            Properties::Pointer linkProperties )
{
    std::vector<unsigned int> vec_foundation_elements;
    std::vector<unsigned int> vec_soil_elements;

    int len_foundation_elements = boost::python::len(foundation_elements);
    for ( int it = 0; it < len_foundation_elements; it++ )
    {
       boost::python::extract<int> x( foundation_elements[it] );

       if ( x.check() )
           vec_foundation_elements.push_back(( unsigned int )x );
       else break;
    }

    int len_soil_elements = boost::python::len(soil_elements);
    for ( int it = 0; it < len_soil_elements; it++ )
    {
        boost::python::extract<int> x( soil_elements[it] );

        if ( x.check() )
           vec_soil_elements.push_back(( unsigned int )x );
        else break;
    }

    dummy.InitializeFoundationUtility( model_part, vec_foundation_elements, vec_soil_elements, linkProperties );
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
    .def( "InitializeWithThreads", &DeactivationUtility::InitializeWithThreads )
    .def( "GetName", &DeactivationUtility::GetName<Element> )
    .def( "GetName", &DeactivationUtility::GetName<Condition> )
    .def( "SetAssociatedElement", &SetAssociatedElement )
    ;

    void(VariableTransferUtility::*pointer_to_TransferPrestressIdentically)(ModelPart&, ModelPart&) = &VariableTransferUtility::TransferPrestressIdentically;
    void(VariableTransferUtility::*pointer_to_TransferPrestressIdenticallyWithCheck)(ModelPart&, ModelPart&) = &VariableTransferUtility::TransferPrestressIdenticallyWithCheck;
    void(VariableTransferUtility::*pointer_to_TransferPrestressIdenticallyForElement)(Element&, Element&, const ProcessInfo&) = &VariableTransferUtility::TransferPrestressIdentically;

    class_<VariableTransferUtility, boost::noncopyable >
    ( "VariableTransferUtility", init<>() )
    .def(init<VariableTransferUtility::LinearSolverType::Pointer>())
    .def( "TransferNodalVariables", &VariableTransferUtility::TransferNodalVariables )
    .def( "TransferNodalVariables", &VariableTransferUtility::TransferGeneralNodalVariables<Variable<double> > )
    .def( "TransferNodalVariables", &VariableTransferUtility::TransferGeneralNodalVariables<Variable<Vector> > )
    .def( "TransferConstitutiveLawVariables", &VariableTransferUtility::TransferConstitutiveLawVariables )
    .def( "TransferInSituStress", &VariableTransferUtility::TransferInSituStress )
    .def( "TransferPrestress", &VariableTransferUtility::TransferPrestress )
    .def( "TransferPrestressIdentically", pointer_to_TransferPrestressIdentically )
    .def( "TransferPrestressIdenticallyWithCheck", pointer_to_TransferPrestressIdenticallyWithCheck )
    .def( "TransferPrestressIdentically", pointer_to_TransferPrestressIdenticallyForElement )
    .def( "TransferSpecificVariable", &VariableTransferUtility::TransferSpecificVariable )
    .def( "TransferSpecificVariableWithComponents", &VariableTransferUtility::TransferSpecificVariableWithComponents )
    .def( "InitializeModelPart", &VariableTransferUtility::InitializeModelPart )
    .def("TransferVariablesToNodes", &DoubleTransferVariablesToNodes)
    .def("TransferVariablesToNodes", &DoubleTransferVariablesToNodesForElements)
    .def("TransferVariablesToNodesForElements", &DoubleTransferVariablesToNodesForElementsAsList)
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
    .def("ComputeExtrapolatedNodalValues", &DoubleComputeExtrapolatedNodalValues)
    .def("ComputeExtrapolatedNodalValues", &Array1DComputeExtrapolatedNodalValues)
    .def("ComputeExtrapolatedNodalValues", &VectorComputeExtrapolatedNodalValues)
    .def("TransferVariablesToGaussPoints", &DoubleTransferVariablesToGaussPoints)
    .def("TransferVariablesToGaussPoints", &DoubleTransferVariablesToGaussPointsLocal)
    .def("TransferVariablesToGaussPoints", &Array1DTransferVariablesToGaussPointsElements)
    .def("TransferVariablesToGaussPoints", &Array1DTransferVariablesToGaussPoints)
    .def("TransferVariablesToGaussPoints", &Array1DTransferVariablesToGaussPointsLocal)
    .def("TransferVariablesToGaussPoints", &VectorTransferVariablesToGaussPoints)
    .def("TransferVariablesToGaussPoints", &VectorTransferVariablesToGaussPointsFromElementsToElements)
    .def("TransferVariablesToGaussPoints", &VectorTransferVariablesToGaussPointsFromElementsToElementsComponents)
    .def("TransferVariablesToGaussPoints", &VectorTransferVariablesToGaussPointsComponents)
    .def("TransferVariablesToGaussPoints", &VectorTransferVariablesToGaussPointsLocal)
    .def("TransferVariablesToGaussPoints", &VectorTransferVariablesToGaussPointsElementComponents)
    .def("TransferVariablesToGaussPoints", &DoubleTransferVariablesToGaussPointsFromNodalValues)
    .def("TransferVariablesToGaussPoints", &Array1DTransferVariablesToGaussPointsFromNodalValues)
    .def("TransferVariablesToGaussPoints", &VectorTransferVariablesToGaussPointsFromNodalValues)
    .def("TransferVariablesToGaussPointsIdentically", &DoubleTransferVariablesToGaussPointsElementComponentsIdentically)
    .def("TransferVariablesToGaussPointsIdentically", &VectorTransferVariablesToGaussPointsElementComponentsIdentically)
    .def("TransferVariablesToGaussPointsIdentically", &VectorTransferVariablesToGaussPointsElementComponentsIdentically2)
    .def("TransferVariablesFromNodeToNode", &VariableTransferUtility::TransferVariablesFromNodeToNode<Variable<Vector> >)
    .def("TransferVariablesFromNodeToNode", &VariableTransferUtility::TransferVariablesFromNodeToNode<Variable<array_1d<double, 3> > >)
    .def("TransferVariablesFromNodeToNode", &VariableTransferUtility::TransferVariablesFromNodeToNode<Variable<double> >)
    .def("TransferVariablesBetweenMeshes", &DoubleTransferVariablesBetweenMeshes)
    .def("TransferVariablesBetweenMeshes", &VectorTransferVariablesBetweenMeshes)
    .def("TransferVariablesToNode", &DoubleTransferVariablesToNode)
    .def("TransferVariablesToNode", &Array1DTransferVariablesToNode)
    ;

    class_<VariableUtility, boost::noncopyable >
    ( "VariableUtility", init<ModelPart::ElementsContainerType&>() )
    .add_property("EchoLevel", &VariableProjectionUtility::GetEchoLevel, &VariableProjectionUtility::SetEchoLevel)
    ;

    class_<VariableProjectionUtility, bases<VariableUtility>, boost::noncopyable >
    ( "VariableProjectionUtility", init<ModelPart::ElementsContainerType&, VariableProjectionUtility::LinearSolverType::Pointer>() )
    .def("TransferVariablesToNodes", &VariableProjectionUtility::TransferVariablesToNodes)
    ;

    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToNodesDouble1)(ModelPart&, const Variable<double>&) = &VariableInterpolationUtility::TransferVariablesToNodes;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToNodesDouble2)(ModelPart::NodesContainerType&, const Variable<double>&) = &VariableInterpolationUtility::TransferVariablesToNodes;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToNodesArray1D1)(ModelPart&, const Variable<array_1d<double, 3> >&) = &VariableInterpolationUtility::TransferVariablesToNodes;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToNodesArray1D2)(ModelPart::NodesContainerType&, const Variable<array_1d<double, 3> >&) = &VariableInterpolationUtility::TransferVariablesToNodes;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToNodesVector1)(ModelPart&, const Variable<Vector>&, const std::size_t&) = &VariableInterpolationUtility::TransferVariablesToNodes;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToNodesVector2)(ModelPart::NodesContainerType&, const Variable<Vector>&, const std::size_t&) = &VariableInterpolationUtility::TransferVariablesToNodes;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToGaussPointsDouble1)(ModelPart&, const Variable<double>&) = &VariableInterpolationUtility::TransferVariablesToGaussPoints;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToGaussPointsDouble2)(ModelPart::ElementsContainerType&, const Variable<double>&, const ProcessInfo&) = &VariableInterpolationUtility::TransferVariablesToGaussPoints;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToGaussPointsArray1D1)(ModelPart&, const Variable<array_1d<double, 3> >&) = &VariableInterpolationUtility::TransferVariablesToGaussPoints;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToGaussPointsArray1D2)(ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3> >&, const ProcessInfo&) = &VariableInterpolationUtility::TransferVariablesToGaussPoints;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToGaussPointsVector1)(ModelPart&, const Variable<Vector>&, std::size_t) = &VariableInterpolationUtility::TransferVariablesToGaussPoints;
    void(VariableInterpolationUtility::*pointer_to_TransferVariablesToGaussPointsVector2)(ModelPart::ElementsContainerType&, const Variable<Vector>&, const ProcessInfo&, std::size_t) = &VariableInterpolationUtility::TransferVariablesToGaussPoints;
    ModelPart::ElementsContainerType(VariableInterpolationUtility::*pointer_to_FindPotentialPartners)(const Element::GeometryType::PointType::PointType&) const = &VariableInterpolationUtility::FindPotentialPartners;
    Element::Pointer(VariableInterpolationUtility::*pointer_to_SearchPartner)(const Element::GeometryType::PointType::PointType&, ModelPart::ElementsContainerType&) const = &VariableInterpolationUtility::SearchPartner;

    class_<VariableInterpolationUtility, bases<VariableUtility>, boost::noncopyable >
    ( "VariableInterpolationUtility", init<ModelPart::ElementsContainerType&>() )
    .def("FindPotentialPartners", pointer_to_FindPotentialPartners)
    .def("SearchPartner", pointer_to_SearchPartner)
    .def("TransferVariablesToNodes", pointer_to_TransferVariablesToNodesDouble1)
    .def("TransferVariablesToNodes", pointer_to_TransferVariablesToNodesDouble2)
    .def("TransferVariablesToNodes", pointer_to_TransferVariablesToNodesArray1D1)
    .def("TransferVariablesToNodes", pointer_to_TransferVariablesToNodesArray1D2)
    .def("TransferVariablesToNodes", pointer_to_TransferVariablesToNodesVector1)
    .def("TransferVariablesToNodes", pointer_to_TransferVariablesToNodesVector2)
    .def("TransferVariablesToGaussPoints", pointer_to_TransferVariablesToGaussPointsDouble1)
    .def("TransferVariablesToGaussPoints", pointer_to_TransferVariablesToGaussPointsDouble2)
    .def("TransferVariablesToGaussPoints", pointer_to_TransferVariablesToGaussPointsArray1D1)
    .def("TransferVariablesToGaussPoints", pointer_to_TransferVariablesToGaussPointsArray1D2)
    .def("TransferVariablesToGaussPoints", pointer_to_TransferVariablesToGaussPointsVector1)
    .def("TransferVariablesToGaussPoints", pointer_to_TransferVariablesToGaussPointsVector2)
    ;

    class_<VariableBinningInterpolationUtility, bases<VariableInterpolationUtility>, boost::noncopyable >
    ( "VariableBinningInterpolationUtility", init<ModelPart::ElementsContainerType&, const double&, const double&, const double&>() )
    ;

    class_<VariableBVHInterpolationUtility, bases<VariableInterpolationUtility>, boost::noncopyable >
    ( "VariableBVHInterpolationUtility", init<ModelPart::ElementsContainerType&, const int&>() )
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

    class_<OutputUtility, boost::noncopyable >
    ( "OutputUtility",
      init<>() )
    .def( "GetStrain", &OutputUtility::GetStrain )
    .def( "GetStress", &OutputUtility::GetStress )
    .def( "GetInternalVariables", &OutputUtility::GetInternalVariables )
    .def( "GetNumberOfPlasticPoints", &OutputUtility::GetNumberOfPlasticPoints )
    .def( "ListPlasticPoints", &OutputUtility::ListPlasticPoints )
    ;

    def( "AddNewRigidBody3D", AddNewRigidBody3D );
    def( "AddNewRigidBodyAndSpring3D", AddNewRigidBodyAndSpring3D );
    ;

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


    class_<DofUtility, boost::noncopyable >
    ( "DofUtility", init<>() )
    .def( "ListDofs", &ListDofs )
    .def( "RemoveDof", &RemoveDof<Variable<double> > )
    .def( "RemoveDof", &RemoveDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > > )
    .def( "PrintKey", &PrintKey<Variable<double> > )
    .def( "PrintKey", &PrintKey<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > > )
    .def( "PrintMaxUnbalancedForce", &PrintMaxUnbalancedForce )
    .def( "PrintMinUnbalancedForce", &PrintMinUnbalancedForce )
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

    class_<TipUtility, boost::noncopyable >
    ( "TipUtility", init<>() )
    .def( "InitializeTipUtility", &InitializeTipUtility )
    ;

    class_<FoundationUtility, boost::noncopyable >
    ( "FoundationUtility", init<>() )
    .def( "InitializeFoundationUtility", &InitializeFoundationUtility )
    ;

}
}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED  defined
