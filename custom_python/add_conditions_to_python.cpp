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
//   Last Modified by:    $Author: nelson $
//   Date:                $Date: 2008-12-09 15:23:36 $
//   Revision:            $Revision: 1.2 $
//
//



// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <cstring>
// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "custom_python/add_conditions_to_python.h"
#include "includes/define.h"
#include "includes/condition.h"
#include "custom_conditions/pointforce3D.h"
#include "custom_conditions/faceforce3D.h"
#include "custom_conditions/point_point_joint_condition.h"
#include "custom_conditions/face_face_joint_condition.h"
#include "custom_conditions/point_point_lagrange_condition.h"
#include "custom_conditions/elastic_constraint.h"
#include "custom_conditions/collinear_constraint_2d.h"
#include "custom_conditions/inclined_constraint_2d.h"
#include "custom_conditions/inclined_constraint_3d.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "python/add_mesh_to_python.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"

namespace Kratos
{

namespace Python
{
using namespace boost::python;
typedef Condition ConditionBaseType;
typedef Geometry<Node<3> > GeometryType;
typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
typedef GeometryType::PointsArrayType NodesArrayType;

void  AddCustomConditionsToPython()
{

    class_< PointForce3D, PointForce3D::Pointer, bases< ConditionBaseType > >
    ("PointForce3D",
     init<int, Node<3>::Pointer, Properties::Pointer>() )
    ;

    class_< FaceForce3D, FaceForce3D::Pointer, bases< ConditionBaseType > >
    ("FaceForce3D",
     init<int, GeometryType::Pointer, Properties::Pointer>() )
    ;

    class_< PointPointJointCondition, PointPointJointCondition::Pointer, bases< ConditionBaseType > >
    ("PointPointJointCondition",
     init<int, Node<3>::Pointer, Node<3>::Pointer, Properties::Pointer>() )
    ;

    class_< FaceFaceJointCondition, FaceFaceJointCondition::Pointer, bases< ConditionBaseType > >
    ("FaceFaceJointCondition",
     init<int, Condition::Pointer, Condition::Pointer, Properties::Pointer>() )
    ;

    class_< PointPointLagrangeCondition, PointPointLagrangeCondition::Pointer, bases< ConditionBaseType > >
    ("PointPointLagrangeCondition",
     init<int, Node<3>::Pointer, Node<3>::Pointer, Properties::Pointer>() )
    ;

    class_< ElasticConstraint, ElasticConstraint::Pointer, bases< ConditionBaseType > >
    ("ElasticConstraint",
     init<int, Node<3>::Pointer, Properties::Pointer>() )
    ;

    class_< CollinearConstraint2D, CollinearConstraint2D::Pointer, bases< ConditionBaseType > >
    ("CollinearConstraint2D",
     init<int, Node<3>::Pointer, Node<3>::Pointer, Node<3>::Pointer, Properties::Pointer>() )
    ;

    class_< InclinedConstraint2D, InclinedConstraint2D::Pointer, bases< ConditionBaseType > >
    ("InclinedConstraint2D",
     init<int, Node<3>::Pointer, const double&, const double&, const double&, Properties::Pointer>() )
    ;

    class_< InclinedConstraint3D, InclinedConstraint3D::Pointer, bases< ConditionBaseType > >
    ("InclinedConstraint3D",
     init<int, Node<3>::Pointer, const double&, const double&, const double&, const double&, Properties::Pointer>() )
    ;

}

}  // namespace Python.

}  // namespace Kratos.
