
#if !defined(KRATOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

void StructuralApplication_AddCustomProcessesToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED  defined 
