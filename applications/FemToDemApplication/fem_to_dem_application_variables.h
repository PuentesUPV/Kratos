//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Vel�zquez
//

#if !defined(KRATOS_FEM_TO_DEM_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_FEM_TO_DEM_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{
	KRATOS_DEFINE_VARIABLE( double, DAMAGE_EDGE1)
	KRATOS_DEFINE_VARIABLE( double, DAMAGE_EDGE2)
	KRATOS_DEFINE_VARIABLE( double, DAMAGE_EDGE3)
	KRATOS_DEFINE_VARIABLE( double, DAMAGE_ELEMENT)
	KRATOS_DEFINE_VARIABLE(Vector, STRESS_VECTOR);
	KRATOS_DEFINE_VARIABLE(double, YIELD_STRESS_C);
	KRATOS_DEFINE_VARIABLE(double, YIELD_STRESS_T);
	KRATOS_DEFINE_VARIABLE(int, ITER);
	KRATOS_DEFINE_VARIABLE(double, FRAC_ENERGY_T)
	KRATOS_DEFINE_VARIABLE(double, FRAC_ENERGY_C)
	KRATOS_DEFINE_VARIABLE(Vector, STRESS_VECTOR_INTEGRATED);
	KRATOS_DEFINE_VARIABLE(double, THRESHOLD)
	KRATOS_DEFINE_VARIABLE(Vector, SMOOTHED_STRESS_VECTOR);
	KRATOS_DEFINE_VARIABLE(std::string, YIELD_SURFACE);
	KRATOS_DEFINE_VARIABLE(Vector, STRAIN_VECTOR);
	KRATOS_DEFINE_VARIABLE(bool, TANGENT_CONSTITUTIVE_TENSOR);
	KRATOS_DEFINE_VARIABLE(bool, SMOOTHING);
	KRATOS_DEFINE_VARIABLE(double, IS_DAMAGED);
	KRATOS_DEFINE_VARIABLE(double, CHARACTERISTIC_LENGTH);
	KRATOS_DEFINE_VARIABLE(int, MESH_REFINED);
	KRATOS_DEFINE_VARIABLE(int, IS_DYNAMIC);
	KRATOS_DEFINE_VARIABLE(int, STRESS_THRESHOLD);
	KRATOS_DEFINE_VARIABLE(int, INTEGRATION_COEFFICIENT);
	KRATOS_DEFINE_VARIABLE(std::string, MAPPING_PROCEDURE);
}

#endif	/* KRATOS_FEM_TO_DEM_APPLICATION_VARIABLES_H_INCLUDED */
