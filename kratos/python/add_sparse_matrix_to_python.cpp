//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//



// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "python/matrix_python_interface.h"
#include "python/matrix_scalar_operator_python.h"
#include "python/matrix_scalar_assignment_operator_python.h"
#include "python/matrix_vector_operator_python.h"

namespace Kratos
{
namespace Python
{

using namespace boost::python;


void  AddSparseMatrixToPython()
{
    /*	  MatrixPythonInterface<mapped_matrix<double> >::CreateInterface("SparseMatrix")
    		  .def(init<mapped_matrix<double>::size_type, mapped_matrix<double>::size_type>())
    		  .def(MatrixScalarOperatorPython<mapped_matrix<double>, double>())
    		  .def(MatrixScalarAssignmentOperatorPython<mapped_matrix<double>, double>())
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, zero_matrix<double>, mapped_matrix<double> >())
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, identity_matrix<double>, mapped_matrix<double> >())
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, scalar_matrix<double>, mapped_matrix<double> >())
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, banded_matrix<double>, banded_matrix<double> >())
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, triangular_matrix<double, upper>, matrix<double> >())
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, triangular_matrix<double, lower>, matrix<double> >())
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, symmetric_matrix<double, upper>, matrix<double> >())
    #if defined KRATOS_ADD_HERMITIAN_MATRIX_INTERFACE
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, hermitian_matrix<double, upper>, matrix<double> >())
    #endif
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, matrix<double>, matrix<double> >())
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, compressed_matrix<double>, compressed_matrix<double> >())
    #if defined KRATOS_ADD_COORDINATE_MATRIX_INTERFACE
    	          .def(MatrixMatrixOperatorPython<mapped_matrix<double>, coordinate_matrix<double>, mapped_matrix<double> >())
    #endif
    		  ;*/
}

}  // namespace Python.

} // Namespace Kratos
