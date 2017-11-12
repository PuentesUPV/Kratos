//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//

#ifndef APPLY_PERIODIC_CONDITION_PROCESS_H
#define APPLY_PERIODIC_CONDITION_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/binbased_fast_point_locator.h"

namespace Kratos
{

class ApplyPeriodicConditionProcess : public Process
{
  public:
    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPeriodicConditionProcess);

    typedef MpcData::Pointer MpcDataPointerType;
    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;
    typedef std::map<std::string, MpcDataPointerType> MpcDataMapType;
    typedef MpcData::VariableComponentType VariableComponentType;
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef unsigned int IndexType;
    typedef std::vector<MpcDataPointerType> *MpcDataPointerVectorType;
    typedef MpcData::VariableType VariableType;
    typedef boost::shared_ptr<std::vector<MpcDataPointerType>> MpcDataSharedPointerVectorType;
    typedef ModelPart::NodeIterator NodeIterator;

    /// Constructor.
    ApplyPeriodicConditionProcess(ModelPart &model_part,
                                  Parameters rParameters) : Process(Flags()), mrMainModelPart(model_part), mParameters(rParameters)
    {

        Parameters default_parameters(R"(
                                            {
                                            "constraint_set_name":"default",
                                            "master_sub_model_part_name":"default_master",
                                            "slave_sub_model_part_name":"default_slave",                
                                            "variable_names":[""],
                                            "center":[],
                                            "axis_of_rotation":[],
                                            "angle":0.0                                                                                                                             
                                            }  )");

        // Initializing

        // using the given data first rotate the master surface to match with slave surface. obtain a new model part.

        // Use the salve model part and newly created model part to form MPC constraints
        mSlaveSubModelPartName = m_parameters["slave_sub_model_part_name"].GetString();
        mMasterSubModelPartName = m_parameters["master_sub_model_part_name"].GetString();

        mCenterOfRotation.push_back(m_parameters["center"][0].GetDouble());
        mCenterOfRotation.push_back(m_parameters["center"][1].GetDouble());
        mCenterOfRotation.push_back(m_parameters["center"][2].GetDouble());

        mAxisOfRoationVector.push_back(m_parameters["axis_of_rotation"][0].GetDouble());
        mAxisOfRoationVector.push_back(m_parameters["axis_of_rotation"][1].GetDouble());
        mAxisOfRoationVector.push_back(m_parameters["axis_of_rotation"][2].GetDouble());

        // normalizing the axis of roatation
        double norm = 0.0;
        for (unsigned int d = 0; d < 3; ++d)
            norm += mAxisOfRoationVector[d] * mAxisOfRoationVector[d];
        norm = sqrt(norm);
        for (unsigned int d = 0; d < 3; ++d)
            mAxisOfRoationVectorNormalized.push_back(mAxisOfRoationVector[d] / norm);

        mTheta = m_parameters["angle"].GetDouble();
        mrRotatedMasterModelPart = ModelPart::Pointer(new ModelPart())
    }
    ~ApplyPeriodicConditionProcess()
    {
    }

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart &masterModelPart = mr_model_part.GetSubModelPart(m_parameters["master_sub_model_part_name"].GetString());
        ModelPart &slaveModelPart = mr_model_part.GetSubModelPart(m_parameters["slave_sub_model_part_name"].GetString());

        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        int &dim = info->GetValue(DOMAIN_SIZE);

        GetRotatedMaster<dim>(masterModelPart);

        KRATOS_CATCH("");
    }

    // Functions which use two variable components
    template <int TDim>
    void GetRotatedMaster(ModelPart &master_model_part)
    {
        BinBasedFastPointLocator<TDim> *p_point_locator = new BinBasedFastPointLocator<TDim>(master_model_part);
        int numVars = m_parameters["variable_names"].size();
        // iterating over slave nodes to find the corresponding masters
        const int n_slave_nodes = slave_model_part.Nodes().size();
        array_1d<double, TDim + 1> N; // This is only for triangular meshes
        const int max_results = 100;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        for (int i = 0; i < n_slave_nodes; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = slave_model_part.NodesBegin() + i;
            Node<3>::Pointer p_slave_node = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
            Element::Pointer pMasterElement;
        }
    }

  private:
    std::vector<std::vector<std::vector<double>>> vectorOfRotationMatrices;
    Parameters mParameters;
    ModelPart &mrMainModelPart;
    ModelPart::Pointer *mrRotatedMasterModelPart; // This is new modelpart
    std::string mSlaveSubModelPartName;
    std::string mMasterSubModelPartName;
    double mTheta;
    std::vector<double> mCenterOfRotation;
    std::vector<double> mAxisOfRoationVector;
    std::vector<double> mAxisOfRoationVectorNormalized;

    /*
     * Calculates the cross product of two vectors
     *  c = axb
     */
    void CrossProduct(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c)
    {
        assert(a.size() == b.size());
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
    }

    /*
     * Rotates a given point(node_cords) in space around a given mAxisOfRoationVector by an angle thetha
     */
    std::vector<double> RotateNode(std::vector<double> node_cords, double theta)
    {
        std::vector<double> rotatedNode(4);
        std::vector<double> U(3); // normalized axis of rotation
        node_cords.push_back(1.0);

        //std::cout<<"############## HERE"<<std::endl;
        // normalizing the axis of roatation
        double norm = 0.0;
        for (unsigned int d = 0; d < 3; ++d)
            norm += mAxisOfRoationVector[d] * mAxisOfRoationVector[d];
        norm = sqrt(norm);
        for (unsigned int d = 0; d < 3; ++d)
            U[d] = mAxisOfRoationVector[d] / norm;

        // Constructing the transformation matrix
        double A0[4][4];
        double x1 = mCenterOfRotation[0];
        double y1 = mCenterOfRotation[1];
        double z1 = mCenterOfRotation[2];

        double a = U[0];
        double b = U[1];
        double c = U[2];

        double t2 = cos(theta);
        double t3 = sin(theta);
        double t4 = a * a;
        double t5 = b * b;
        double t6 = c * c;
        double t7 = a * b;
        double t8 = t5 + t6;
        double t9 = 1.0 / t8;
        double t10 = a * c;
        double t11 = b * t3;
        double t12 = a * t3 * t5;
        double t13 = a * t3 * t6;
        double t14 = b * c * t2;
        A0[0][0] = t4 + t2 * t8;
        A0[0][1] = t7 - c * t3 - a * b * t2;
        A0[0][2] = t10 + t11 - a * c * t2;
        A0[0][3] = x1 - t4 * x1 - a * b * y1 - a * c * z1 - b * t3 * z1 + c * t3 * y1 - t2 * t5 * x1 - t2 * t6 * x1 + a * b * t2 * y1 + a * c * t2 * z1;
        A0[1][0] = t7 + c * t3 - a * b * t2;
        A0[1][1] = t9 * (t2 * t6 + t5 * t8 + t2 * t4 * t5);
        A0[1][2] = -t9 * (t12 + t13 + t14 - b * c * t8 - b * c * t2 * t4);
        A0[1][3] = -t9 * (-t8 * y1 + t2 * t6 * y1 + t5 * t8 * y1 + a * b * t8 * x1 - b * c * t2 * z1 + b * c * t8 * z1 - a * t3 * t5 * z1 - a * t3 * t6 * z1 + c * t3 * t8 * x1 + t2 * t4 * t5 * y1 - a * b * t2 * t8 * x1 + b * c * t2 * t4 * z1);
        A0[2][0] = t10 - t11 - a * c * t2;
        A0[2][1] = t9 * (t12 + t13 - t14 + b * c * t8 + b * c * t2 * t4);
        A0[2][2] = t9 * (t2 * t5 + t6 * t8 + t2 * t4 * t6);
        A0[2][3] = -t9 * (-t8 * z1 + t2 * t5 * z1 + t6 * t8 * z1 + a * c * t8 * x1 - b * c * t2 * y1 + b * c * t8 * y1 + a * t3 * t5 * y1 + a * t3 * t6 * y1 - b * t3 * t8 * x1 + t2 * t4 * t6 * z1 - a * c * t2 * t8 * x1 + b * c * t2 * t4 * y1);
        A0[3][3] = 1.0;

        // Multiplying the point to get the rotated point
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                rotatedNode[i] += (A0[i][j] * node_cords[j]);
            }
        }

        return rotatedNode;
    }
};
}

#endif