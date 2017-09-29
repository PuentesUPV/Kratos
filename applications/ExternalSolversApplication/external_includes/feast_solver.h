//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:		 BSD License
//   Kratos default license: kratos/license.txt
//
//   Main authors:    Michael Andre
//     Co-authors:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <iostream>
#include <complex>
#include <vector>
#include <unordered_set>
#include <algorithm>

// External includes
#include <boost/smart_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/value_type/complex.hpp>
#include <amgcl/solver/skyline_lu.hpp>
extern "C" 
{
    #include "feast.h"
}

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/direct_solver.h"
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"

#if !defined(KRATOS_FEAST_SOLVER)
#define  KRATOS_FEAST_SOLVER

namespace Kratos 
{

///@name Kratos Classes
///@{

template <class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class SkylineLUSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SkylineLUSolver);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TSparseSpaceType::DataType DataType;
    
    typedef typename amgcl::backend::builtin<DataType>::matrix BuiltinMatrixType;
    
    typedef amgcl::solver::skyline_lu<DataType> SolverType;

    ~SkylineLUSolver() override
    {
        Clear();
    }

    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        Clear();

        pBuiltinMatrix = amgcl::adapter::zero_copy(
                rA.size1(),
                rA.index1_data().begin(),
                rA.index2_data().begin(),
                rA.value_data().begin());

        pSolver = boost::make_shared<SolverType>(*pBuiltinMatrix);
    }

    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        std::vector<DataType> x(rX.size());
        std::vector<DataType> b(rB.size());

        std::copy(std::begin(rB), std::end(rB), std::begin(b));

        (*pSolver)(b, x);

        std::copy(std::begin(x), std::end(x), std::begin(rX));

        return true;
    }

    void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        Clear();
    }

    void Clear() override
    {
        pSolver.reset();
        pBuiltinMatrix.reset();
    }

private:

    boost::shared_ptr<BuiltinMatrixType> pBuiltinMatrix;

    boost::shared_ptr<SolverType> pSolver;
};

/// Adapter to FEAST eigenvalue problem solver.
template<class TSparseSpaceType, class TDenseSpaceType,
        class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class FEASTSolver: public LinearSolver<TSparseSpaceType, TDenseSpaceType,
        TReordererType> {

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( FEASTSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType SparseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef std::complex<double> ComplexType;

    typedef boost::numeric::ublas::compressed_matrix<ComplexType> ComplexSparseMatrixType;

    typedef boost::numeric::ublas::matrix<ComplexType> ComplexDenseMatrixType;

    typedef boost::numeric::ublas::vector<ComplexType> ComplexVectorType;

    typedef UblasSpace<ComplexType, ComplexSparseMatrixType, ComplexVectorType> ComplexSparseSpaceType;

    typedef UblasSpace<ComplexType, ComplexDenseMatrixType, ComplexVectorType> ComplexDenseSpaceType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for built-in linear solver.
    /**
     * Parameters let the user control the settings of the FEAST library.
     */
    FEASTSolver(Parameters::Pointer pParam) : mpParam(pParam)
    {
        Parameters default_params(R"(
        {
            "solver_type": "FEAST",
            "print_feast_output": false,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem": true,
            "lambda_min": 0.0,
            "lambda_max": 1.0,
            "echo_level": 1,
            "number_of_eigenvalues": 0,
            "search_dimension": 10,
            "linear_solver_settings": {
                "solver_type": "skyline_lu"
            }
        })");

        mpParam->RecursivelyValidateAndAssignDefaults(default_params);

        if (mpParam->GetValue("linear_solver_settings")["solver_type"].GetString() != "skyline_lu")
        {
            KRATOS_ERROR << "Built-in solver type must be used with this constructor" << std::endl;
        }

        mpLinearSolver = boost::make_shared<SkylineLUSolver<ComplexSparseSpaceType, ComplexDenseSpaceType>>();
    }

    /// Constructor for externally provided linear solver.
    /**
     * Parameters let the user control the settings of the FEAST library.
     * Warning: For iterative solvers, very small tolerances (~1e-15)
     *          may be needed for FEAST to work properly. Common iterative 
     *          solvers normally don't perform efficiently with FEAST 
     *          (M. Galgon et al., Parallel Computing (49) 2015 153-163).
     */
    FEASTSolver(
        Parameters::Pointer pParam, 
        ComplexLinearSolverType::Pointer pLinearSolver
        ): mpParam(pParam)
    {
        Parameters default_params(R"(
        {
            "solver_type": "FEAST",
            "print_feast_output": false,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem": true,
            "lambda_min": 0.0,
            "lambda_max": 1.0,
            "echo_level": 1,
            "number_of_eigenvalues": 0,
            "search_dimension": 10,
            "linear_solver_settings": {}
        })");

        // Don't validate linear_solver_settings here
        mpParam->ValidateAndAssignDefaults(default_params);
        
        if (pLinearSolver != nullptr)
        {
            mpLinearSolver = pLinearSolver;
        }
        else
        {
            mpLinearSolver = boost::make_shared<SkylineLUSolver<ComplexSparseSpaceType, ComplexDenseSpaceType>>();
        }
        
    }

    /// Deleted copy constructor.
    FEASTSolver(const FEASTSolver& Other) = delete;

    /// Destructor.
    ~FEASTSolver() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    FEASTSolver& operator=(const FEASTSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * Solve the generalized eigenvalue problem.
     * @param K: is a symmetric matrix. 
     * @param M: Is a symmetric positive-definite matrix.
     * @param Eigenvalues: The vector containing all the eigenvalues
     * @param Eigenvector: The matrix containing the eigen vectors
     */
    void Solve(
        SparseMatrixType& K,
        SparseMatrixType& M,
        DenseVectorType& Eigenvalues,
        DenseMatrixType& Eigenvectors
        ) override
    {
        const auto system_size = K.size1();

        Parameters& feast_settings = *mpParam;
        const double eigen_value_range_min = feast_settings["lambda_min"].GetDouble();
        const double eigen_value_range_max = feast_settings["lambda_max"].GetDouble();

        int search_dimension = feast_settings["search_dimension"].GetInt();
        int num_eigenvalues = feast_settings["number_of_eigenvalues"].GetInt();
        const int echo_level = feast_settings["echo_level"].GetInt();

        Eigenvalues.resize(search_dimension,false);
        Eigenvectors.resize(search_dimension,system_size,false);

        if (feast_settings["perform_stochastic_estimate"].GetBool())
        {
            // This estimates the number of eigenvalues in the interval [lambda_min, lambda_max]
            Calculate(M,K,eigen_value_range_min,eigen_value_range_max,search_dimension,
                    num_eigenvalues,Eigenvalues,Eigenvectors,true);

            if (echo_level > 0)
            {
                std::cout << "Estimated number of eigenvalues = " << num_eigenvalues << std::endl;
            }

            // Recommended estimate of search dimension from FEAST documentation
            search_dimension = num_eigenvalues + num_eigenvalues/2 + 1;
            feast_settings["search_dimension"].SetInt(search_dimension);
        }
        if (feast_settings["solve_eigenvalue_problem"].GetBool())
        {
            // This attempts to solve the generalized eigenvalue problem
            Calculate(M,K,eigen_value_range_min,eigen_value_range_max,search_dimension,
                    num_eigenvalues,Eigenvalues,Eigenvectors,false);

            Eigenvalues.resize(num_eigenvalues,true);
            Eigenvectors.resize(num_eigenvalues,system_size,true);
        }
        feast_settings["number_of_eigenvalues"].SetInt(num_eigenvalues);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FEAST solver.";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Parameters::Pointer mpParam;

    ComplexLinearSolverType::Pointer mpLinearSolver;

    ///@}
    ///@name Private Operations
    ///@{

    /// Wrapper for FEAST library.
    void Calculate(
            SparseMatrixType& rMassMatrix,
            SparseMatrixType& rStiffnessMatrix,
            double EigenvalueRangeMin,
            double EigenvalueRangeMax,
            int SearchDimension,
            int& rNumEigenvalues,
            DenseVectorType& rEigenvalues,
            DenseMatrixType& rEigenvectors,
            bool PerformStochasticEstimate
            )
    {
        KRATOS_TRY

        int feast_params[64] = {};
        int num_iter, info, system_size;
        double eps_out;
        DenseVectorType residual(SearchDimension);
        std::vector<std::complex<double> > integration_nodes, integration_weights;
        system_size = static_cast<int>(rMassMatrix.size1());
        matrix<double,column_major> work(system_size,SearchDimension);
        matrix<std::complex<double>,column_major> zwork(system_size,SearchDimension);
        matrix<double,column_major> Aq(SearchDimension,SearchDimension);
        matrix<double,column_major> Bq(SearchDimension,SearchDimension);
        std::complex<double> Ze;
        ComplexSparseMatrixType Az;
        ComplexVectorType b(system_size);
        ComplexVectorType x(system_size);

        this->InitializeFEASTSystemMatrix(rMassMatrix, rStiffnessMatrix, Az);

        Parameters& feast_settings = *mpParam;

        // initialize FEAST eigenvalue solver (see FEAST documentation for details)
        feastinit(feast_params);
        if (feast_settings["print_feast_output"].GetBool())
            feast_params[0] = 1;
        feast_params[2] = 8; // stopping convergence criteria 10^-feast_params[2]
        feast_params[28] = 1;// not sure if this is needed
        if (PerformStochasticEstimate)
        {
            feast_params[1] = 4; // number of quadrature points (default: 8)
            feast_params[13] = 2;
        }

        integration_nodes.resize(feast_params[1]);
        integration_weights.resize(feast_params[1]);

        // get quadrature nodes and weights
        zfeast_contour(&EigenvalueRangeMin,
                &EigenvalueRangeMax,
                &feast_params[1],
                &feast_params[15],
                &feast_params[17],
                (double *)integration_nodes.data(),
                (double *)integration_weights.data());

        int ijob = -1;
        // solve the eigenvalue problem
        while (ijob != 0)
        {
            // FEAST's reverse communication interface
            dfeast_srcix(&ijob,&system_size,(double *)&Ze,(double *)work.data().begin(),
                    (double *)zwork.data().begin(),(double *)Aq.data().begin(),
                    (double *)Bq.data().begin(),feast_params,&eps_out,&num_iter,
                    &EigenvalueRangeMin,&EigenvalueRangeMax,&SearchDimension,
                    (double *)rEigenvalues.data().begin(),
                    (double *)rEigenvectors.data().begin(),
                    &rNumEigenvalues,(double *)residual.data().begin(),&info,
                    (double *)integration_nodes.data(),
                    (double *)integration_weights.data());

            switch (ijob)
            {
                case 10:
                {
                    // set up quadrature matrix (ZeM-K) and solver
                    this->CalculateFEASTSystemMatrix(Ze, rMassMatrix, rStiffnessMatrix, Az);
                    mpLinearSolver->Clear();
                    mpLinearSolver->Initialize(Az,x,b);
                    mpLinearSolver->InitializeSolutionStep(Az, x, b);
                } break;
                case 11:
                {
                    // solve the linear system for one quadrature point
                    for (int j=0; j < feast_params[22]; j++)
                    {
                        for (int i=0; i < system_size; i++)
                            b[i] = zwork(i,j);
                        mpLinearSolver->Solve(Az,x,b);
                        for (int i=0; i < system_size; i++)
                            zwork(i,j) = x[i];
                    }
                } break;
                case 30:
                {
                    // multiply Kx
                    for (int i=0; i < feast_params[24]; i++)
                    {
                        int k = feast_params[23]-1+i;
                        noalias(column(work,k)) = prod(rStiffnessMatrix,row(rEigenvectors,k));
                    }
                } break;
                case 40:
                {
                    // multiply Mx
                    for (int i=0; i < feast_params[24]; i++)
                    {
                        int k = feast_params[23]-1+i;
                        noalias(column(work,k)) = prod(rMassMatrix,row(rEigenvectors,k));
                    }
                }
            } // switch
        } // while

        KRATOS_CATCH("")
    }

    /**
     * Initialize CSR matrix structure for FEAST system matrix: C = z * B - A.
     */
    void InitializeFEASTSystemMatrix(
        const SparseMatrixType& B,
        const SparseMatrixType& A,
        ComplexSparseMatrixType& C
        )
    {
        C.resize(B.size1(), B.size2(), false);

        std::vector<std::unordered_set<std::size_t> > indices(C.size1());

        // indices for row begin / end
        C.index1_data()[0] = 0;
        for (std::size_t i = 0; i < C.size1(); ++i)
        {
            std::size_t row_begin, row_end;
            indices[i].reserve(40); // initialize C's indices

            row_begin = B.index1_data()[i];
            row_end = B.index1_data()[i + 1];
            indices[i].insert(B.index2_data().begin() + row_begin,
                    B.index2_data().begin() + row_end); // insert B's column indices for row i

            row_begin = A.index1_data()[i];
            row_end = A.index1_data()[i + 1];
            indices[i].insert(A.index2_data().begin() + row_begin,
                    A.index2_data().begin() + row_end); // insert A's column indices for row i

            // C.index1_data()[i+1] = number of non-zeros in rows <= i
            C.index1_data()[i + 1] = C.index1_data()[i] + indices[i].size();
        }

        // C.index1_data()[C.size1()] = number of non-zeros
        C.reserve(C.index1_data()[C.size1()]);

        // column indices
        std::size_t k = 0;
        for (std::size_t i = 0; i < C.size1(); ++i)
        {
            for (std::size_t j : indices[i])
                C.index2_data()[k++] = j; // fill C's column indices

            indices[i].clear();

            std::sort(C.index2_data().begin() + C.index1_data()[i],
                    C.index2_data().begin() + C.index1_data()[i + 1]);
        }

        C.set_filled(C.size1() + 1, C.index1_data()[C.size1()]);
    }

    /**
     * Calculate FEAST system matrix: C = z * B - A. Similar to FEAST's zdaddcsr subroutine.
     */
    void CalculateFEASTSystemMatrix(
        std::complex<double> z,
        SparseMatrixType& B,
        SparseMatrixType& A,
        ComplexSparseMatrixType& C
        )
    {
        std::size_t jb, ja;
        const std::size_t dimension = B.size1();

        std::size_t ptr = 0;
        for (std::size_t i = 0; i < dimension; ++i)
        {
            std::size_t b_ptr = B.index1_data()[i];
            std::size_t a_ptr = A.index1_data()[i];
            while (b_ptr < B.index1_data()[i + 1] || a_ptr < A.index1_data()[i + 1])
            {
                jb = (b_ptr < B.index1_data()[i + 1]) ?
                        B.index2_data()[b_ptr] : dimension;
                ja = (a_ptr < A.index1_data()[i + 1]) ?
                        A.index2_data()[a_ptr] : dimension;

                if (jb < ja)
                {
                    C.value_data()[ptr] = z * B(i, jb).ref();
                    b_ptr++;
                }
                else if (jb > ja)
                {
                    C.value_data()[ptr] = -A(i, ja).ref();
                    a_ptr++;
                }
                else
                { // jb == ja
                    C.value_data()[ptr] = z * B(i, jb).ref() - A(i, ja).ref();
                    b_ptr++;
                    a_ptr++;
                }
                ptr++;
            }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

}; // Class FEASTSolver

/// This utility uses the FEAST solver to obtain (estimate) the the condition number of a regular matrix
/**
 * Regular matrix: A*A^H=A^H*A
 */
template<class TSparseSpace = UblasSpace<double, CompressedMatrix, Vector>,
         class TDenseSpace = UblasSpace<double, Matrix, Vector>
         >
class FEASTConditionNumberUtility
{
public:

    ///@name Type Definitions
    ///@{
    
    typedef Matrix                                             MatrixType;

    typedef Vector                                             VectorType;

    typedef std::size_t                                          SizeType;
    
    typedef std::size_t                                         IndexType;

    typedef typename TSparseSpace::MatrixType            SparseMatrixType;

    typedef typename TSparseSpace::VectorType            SparseVectorType;

    typedef typename TDenseSpace::MatrixType              DenseMatrixType;

    typedef typename TDenseSpace::VectorType              DenseVectorType;
    
    typedef std::complex<double>                              ComplexType;
    
    typedef compressed_matrix<ComplexType>        ComplexSparseMatrixType;

    typedef matrix<ComplexType>                    ComplexDenseMatrixType;

    typedef vector<ComplexType>                         ComplexVectorType;

    typedef UblasSpace<ComplexType, ComplexSparseMatrixType, ComplexVectorType> ComplexSparseSpaceType;

    typedef UblasSpace<ComplexType, ComplexDenseMatrixType, ComplexVectorType> ComplexDenseSpaceType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    FEASTConditionNumberUtility()= default;

    /// Destructor.
    virtual ~FEASTConditionNumberUtility()= default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    FEASTConditionNumberUtility& operator=(FEASTConditionNumberUtility const& rOther)
    = default;

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * Computes the condition number using the maximum and minimum eigenvalue of the system (in moduli). It always uses the default linear solver
     * @param InputMatrix: The matrix to obtain the condition number
     * @return condition_number: The condition number obtained
     */
    double GetConditionNumber(const MatrixType& InputMatrix)
    {
        return ConditionNumber(InputMatrix);
    }
    
    /**
     * Computes the condition number using the maximum and minimum eigenvalue of the system (in moduli)
     * @param InputMatrix: The matrix to obtain the condition number
     * @param pLinearSolver: The complex linear solver considered in the FEAST solver
     * @return condition_number: The condition number obtained
     */
    static inline double ConditionNumber(
        const MatrixType& InputMatrix,
        ComplexLinearSolverType::Pointer pLinearSolver = nullptr
        )
    {
        typedef FEASTSolver<TSparseSpace, TDenseSpace> FEASTSolverType;
        
        Parameters this_params(R"(
        {
            "print_feast_output": false,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem": true,
            "lambda_min": 0.0,
            "lambda_max": 1.0,
            "echo_level": 0,
            "number_of_eigenvalues": 0,
            "search_dimension": 10
        })");
        
        const std::size_t size = InputMatrix.size1();
        
        const double normA = TSparseSpace::TwoNorm(InputMatrix);
        this_params["lambda_max"].SetDouble( normA);
        this_params["lambda_min"].SetDouble(-normA);
        double aux_double = 2.0/3.0 * size;
        int aux_int = static_cast<int>(aux_double) - 1;
        this_params["number_of_eigenvalues"].SetInt(aux_int);
        aux_double = 3.0/2.0 * size;
        aux_int = static_cast<int>(aux_double) + 1;
        this_params["search_dimension"].SetInt(aux_int);
        SparseMatrixType copy_matrix = InputMatrix;
        SparseMatrixType identity_matrix = IdentityMatrix(size, size);
        
        // Create the auxilary eigen system
        DenseMatrixType eigen_vectors;
        DenseVectorType eigen_values;
        
        // Create the FEAST solver
        FEASTSolverType FEASTSolver(boost::make_shared<Parameters>(this_params), pLinearSolver);
        
        // Solve the problem
        FEASTSolver.Solve(copy_matrix, identity_matrix, eigen_values, eigen_vectors);
        
        // Size of the eigen values vector
        const int dim_eigen_values = eigen_values.size();
        
        // We get the moduli of the eigen values
        #pragma omp parallel for 
        for (int i = 0; i < dim_eigen_values; i++)
        {
            eigen_values[i] = std::abs(eigen_values[i]);
        }
        
        // Now we sort the vector
        std::sort(eigen_values.begin(), eigen_values.end());
        
        // We compute the eigen value
        double condition_number = 0.0;
        if (dim_eigen_values > 0) condition_number = eigen_values[dim_eigen_values - 1]/eigen_values[0];
        
        return condition_number;
    }
    
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

private:
    
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

}; /* Class FEASTConditionNumberUtility */

///@}

///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(std::istream& rIStream,
        FEASTSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const FEASTSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}// namespace Kratos.

#endif // KRATOS_FEAST_SOLVER  defined
