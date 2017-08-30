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
//                   Riccardo Rossi
//



// System includes

// External includes
#include <boost/python.hpp>
#include "python/add_parallel_linear_solvers_to_python.h"

//nothing will be compiled if an openmp compiler is not found
#ifdef _OPENMP

// Project includes
#include "includes/define.h"
#include "python/add_equation_systems_to_python.h"
#include "linear_solvers/cg_solver.h"
#include "linear_solvers/cg_solver.h"
#include "linear_solvers/bicgstab_solver.h"
#include "linear_solvers/tfqmr_solver.h"
#include "includes/dof.h"
#include "spaces/parallel_ublas_space.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
//#include "linear_solvers/superlu_solver.h"


#endif


namespace Kratos
{

namespace Python
{

void  AddParallelLinearSolversToPython()
{
//nothing will be compiled if an openmp compiler is not found
#ifdef _OPENMP

    typedef ParallelUblasSpace<double, CompressedMatrix, Vector> ParallelSpaceType;
    typedef UblasSpace<double, Matrix, Vector> ParallelLocalSpaceType;
    typedef LinearSolver<ParallelSpaceType,  ParallelLocalSpaceType> ParallelLinearSolverType;
    typedef IterativeSolver<ParallelSpaceType,  ParallelLocalSpaceType> ParallelIterativeSolverType;
    typedef CGSolver<ParallelSpaceType,  ParallelLocalSpaceType> ParallelCGSolverType;
    //typedef BICGSTABSolver<SpaceType,  LocalSpaceType> BICGSTABSolverType;
    //typedef TFQMRSolver<SpaceType,  LocalSpaceType> TFQMRSolverType;

    bool (ParallelLinearSolverType::*pointer_to_solve)(ParallelLinearSolverType::SparseMatrixType& rA, ParallelLinearSolverType::VectorType& rX, ParallelLinearSolverType::VectorType& rB) = &ParallelLinearSolverType::Solve;

    using namespace boost::python;

    //****************************************************************************************************
    //preconditioners
    //****************************************************************************************************
    typedef Preconditioner<ParallelSpaceType,  ParallelLocalSpaceType> ParallelPreconditionerType;

    class_<ParallelPreconditionerType, ParallelPreconditionerType::Pointer>("ParallelPreconditioner")
    .def(self_ns::str(self))
    ;

    typedef DiagonalPreconditioner<ParallelSpaceType,  ParallelLocalSpaceType> ParallelDiagonalPreconditionerType;
    class_<ParallelDiagonalPreconditionerType, ParallelDiagonalPreconditionerType::Pointer, bases<ParallelPreconditionerType> >("ParallelDiagonalPreconditioner")
    .def(self_ns::str(self))
    ;


    //****************************************************************************************************
    //linear solvers
    //****************************************************************************************************
    class_<ParallelLinearSolverType, ParallelLinearSolverType::Pointer>("ParallelLinearSolver")
    .def("Initialize",&ParallelLinearSolverType::Initialize)
    .def("Solve",pointer_to_solve)
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<ParallelIterativeSolverType, ParallelIterativeSolverType::Pointer, bases<ParallelLinearSolverType> >("ParallelIterativeSolver")
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<ParallelCGSolverType, ParallelCGSolverType::Pointer, bases<ParallelIterativeSolverType> >("ParallelCGSolver")
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(init<double, unsigned int,  ParallelPreconditionerType::Pointer>())
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    //  class_<BICGSTABSolverType, BICGSTABSolverType::Pointer, bases<IterativeSolverType> >("ParallelBICGSTABSolver")
    //		  .def(init<double>())
    //		  .def(init<double, unsigned int>())
    //		  .def(self_ns::str(self))
    //		  .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    //		  ;
    //
    //  class_<TFQMRSolverType, TFQMRSolverType::Pointer, bases<IterativeSolverType> >("ParallelTFQMRSolver")
    //		  .def(init<double>())
    //		  .def(init<double, unsigned int>())
    //		  .def(self_ns::str(self))
    //		  .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    //		  ;

    typedef Reorderer<ParallelSpaceType,  ParallelLocalSpaceType > ParallelReordererType;
    typedef DirectSolver<ParallelSpaceType,  ParallelLocalSpaceType, ParallelReordererType > ParallelDirectSolverType;
    typedef SkylineLUFactorizationSolver<ParallelSpaceType,  ParallelLocalSpaceType, ParallelReordererType > ParallelSkylineLUFactorizationSolverType;

    class_<ParallelReordererType, ParallelReordererType::Pointer >("ParallelReorderer")
    .def( init< >() )
    .def(self_ns::str(self))
    .def( "Initialize",&ParallelReordererType::Initialize)
    .def( "Reorder",&ParallelReordererType::Reorder)
    .def( "InverseReorder",&ParallelReordererType::InverseReorder)
    ;

    class_<ParallelDirectSolverType, ParallelDirectSolverType::Pointer, bases<ParallelLinearSolverType> >("ParallelDirectSolver")
    .def( init< >() )
    .def(self_ns::str(self))
    ;

    class_<ParallelSkylineLUFactorizationSolverType, ParallelSkylineLUFactorizationSolverType::Pointer, bases<ParallelDirectSolverType> >("ParallelSkylineLUFactorizationSolver")
    .def(init< >())
    .def(self_ns::str(self))
    ;

    //	typedef SuperLUSolver<SparseSpaceType, LocalSpaceType> SuperLUSolverType;
    //	class_<SuperLUSolverType, bases<DirsectSolverType>, boost::noncopyable >
    //		( "SuperLUSolver",
    //			init<>() )


//it will compile nothing if openmp is not found
#endif
}

}  // namespace Python.

} // Namespace Kratos


