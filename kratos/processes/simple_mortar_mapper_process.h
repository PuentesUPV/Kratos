//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SIMPLE_MORTAR_MAPPER_PROCESS)
#define KRATOS_SIMPLE_MORTAR_MAPPER_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "geometries/point.h"

/* Custom includes */
#include "includes/mortar_classes.h"

/* Custom utilities */
#include "utilities/exact_mortar_segmentation_utility.h"
#include "utilities/mortar_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
///@}
///@name  Enum's
///@{

#if !defined(HISTORICAL_VALUES)
#define HISTORICAL_VALUES
    enum HistoricalValues {Historical = 0, NonHistorical = 1};
#endif
    
///@}
///@name  Functions
///@{
    
template< int TDim, int TNumNodes, class TVarType, HistoricalValues THist> 
class SimpleMortarMapperProcess
        : public Process
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of SimpleMortarMapperProcess
    KRATOS_CLASS_POINTER_DEFINITION(SimpleMortarMapperProcess);
    
    typedef Point<3>                                     PointType;
    typedef Node<3>                                       NodeType;
    typedef Geometry<NodeType>                        GeometryType;
    typedef Geometry<PointType>                  GeometryPointType;
    typedef ModelPart::NodesContainerType           NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    
    // Type definition for integration methods
    typedef GeometryData::IntegrationMethod      IntegrationMethod;
    
    // Auxiliar geometries
    typedef Line2D2<PointType>                            LineType;
    typedef Triangle3D3<PointType>                    TriangleType;
    
    // Component type
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SimpleMortarMapperProcess( 
        ModelPart& rThisModelPart,
        TVarType& ThisVariable
        ): mrThisModelPart(rThisModelPart),
           mOriginVariable(ThisVariable),
           mDestinyVariable(ThisVariable)
    {
    }
    
    SimpleMortarMapperProcess( 
        ModelPart& rThisModelPart,
        TVarType& OriginVariable,
        TVarType& DestinyVariable
        ): mrThisModelPart(rThisModelPart),
           mOriginVariable(OriginVariable),
           mDestinyVariable(DestinyVariable)
    {
    }

    /// Destructor.
    ~SimpleMortarMapperProcess() override
    = default;

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

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }
    
    ///@}
    ///@name Operations
    ///@{
    
    void Execute() override
    {
        KRATOS_TRY;

        // Create and initialize condition variables:
        MortarKinematicVariables<TNumNodes> rVariables;
    
        // Create the mortar operators
        MortarOperator<TNumNodes> rThisMortarConditionMatrices;
        
        // We call the exact integration utility
        ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(TDim);
        
        // Initialize to zero // TODO: Add OMP
        for(int i=0; i<static_cast<int>(mrThisModelPart.Conditions().size()); i++)
        {
            auto it_cond = mrThisModelPart.ConditionsBegin() + i;
            
            if (it_cond->Is(SLAVE) == true)
            {
                boost::shared_ptr<ConditionMap>& all_conditions_maps = it_cond->GetValue( MAPPING_PAIRS );
                
                for (auto it_pair = all_conditions_maps->begin(); it_pair != all_conditions_maps->end(); ++it_pair )
                {
                    Condition::Pointer p_cond_master = (it_pair->first); // MASTER
                    GeometryType& master_geometry = p_cond_master->GetGeometry();
                    
                    for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                    {
                        MortarUtilities::ResetValue<TVarType, THist>(master_geometry[i_node], mDestinyVariable);
                    }
                }
            }
        }
        
        // We map the values from one side to the other // TODO: Add OMP
        for(int i=0; i<static_cast<int>(mrThisModelPart.Conditions().size()); i++)
        {
            auto it_cond = mrThisModelPart.ConditionsBegin() + i;
            
            if (it_cond->Is(SLAVE) == true)
            {
                const array_1d<double, 3>& slave_normal = it_cond->GetValue(NORMAL);
                GeometryType& slave_geometry = it_cond->GetGeometry();
                
                boost::shared_ptr<ConditionMap>& all_conditions_maps = it_cond->GetValue( MAPPING_PAIRS );
                
                for (auto it_pair = all_conditions_maps->begin(); it_pair != all_conditions_maps->end(); ++it_pair )
                {
                    Condition::Pointer p_cond_master = (it_pair->first); // MASTER
                    const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                    GeometryType& master_geometry = p_cond_master->GetGeometry();
                        
                    IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_2;

                    // Reading integration points
                    std::vector<array_1d<PointType,TDim>> conditions_points_slave;
                    const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, slave_normal, master_geometry, master_normal, conditions_points_slave);
                    
                    if (is_inside == true)
                    {
                        // Initialize general variables for the current master element
                        rVariables.Initialize();
                
                        // Initialize the mortar operators
                        rThisMortarConditionMatrices.Initialize();
                        
                        for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); i_geom++)
                        {
                            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                            for (unsigned int i_node = 0; i_node < TDim; i_node++)
                            {
                                PointType global_point;
                                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                                points_array[i_node] = boost::make_shared<PointType>(global_point);
                            }
                            
                            typename std::conditional<TDim == 2, LineType, TriangleType >::type decomp_geom( points_array );
                            
                            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
                            
                            if (bad_shape == false)
                            {
                                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                                
                                // Integrating the mortar operators
                                for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); point_number++ )
                                {
                                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                                    PointType local_point_parent;
                                    PointType gp_global;
                                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                                    slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
                                    
                                    // Calculate the kinematic variables
                                    const PointType& local_point = integration_points_slave[point_number].Coordinates();

                                    /// SLAVE CONDITION ///
                                    slave_geometry.ShapeFunctionsValues( rVariables.NSlave, local_point.Coordinates() );
                                    rVariables.PhiLagrangeMultipliers = rVariables.NSlave;
                                    rVariables.DetjSlave = slave_geometry.DeterminantOfJacobian( local_point );
                                    
                                    /// MASTER CONDITION ///
                                    PointType projected_gp_global;
                                    const array_1d<double,3> gp_normal = MortarUtilities::GaussPointNormal(rVariables.NSlave, slave_geometry);
                                    
                                    GeometryType::CoordinatesArrayType slave_gp_global;
                                    slave_geometry.GlobalCoordinates( slave_gp_global, local_point );
                                    MortarUtilities::FastProjectDirection( master_geometry, slave_gp_global, projected_gp_global, master_normal, -gp_normal ); // The opposite direction
                                    
                                    GeometryType::CoordinatesArrayType projected_gp_local;
                                    
                                    master_geometry.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;
                                    
                                    // SHAPE FUNCTIONS 
                                    master_geometry.ShapeFunctionsValues( rVariables.NMaster, projected_gp_local );    
                                    
                                    const double integration_weight = integration_points_slave[point_number].Weight();
                                    
                                    rThisMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);   
                                }
                                
                                double aux_det;
                                const bounded_matrix<double, TNumNodes, TNumNodes> inv_d_operator = MathUtils<double>::InvertMatrix<TNumNodes>(rThisMortarConditionMatrices.DOperator, aux_det);
                                const bounded_matrix<double, TNumNodes, TNumNodes> p_operator = prod(inv_d_operator, rThisMortarConditionMatrices.MOperator); 
                                
                                Matrix var_origin_matrix;
                                MortarUtilities::MatrixValue<TVarType, THist>(slave_geometry, mDestinyVariable, var_origin_matrix);
                
                                const Matrix var_destiny_matrix = prod(p_operator, var_origin_matrix);
                                MortarUtilities::AddValue<TVarType, THist>(master_geometry, mDestinyVariable, var_destiny_matrix);
                            }
                        }
                    }
                }
            }
        }
        
        KRATOS_CATCH("");
    }

    /**
     * This method sets both variables (origin and destination) with the same variable
     */
    void SetVariable(TVarType ThisVariable)
    {
        mOriginVariable = ThisVariable;
        mDestinyVariable = ThisVariable;
    }
    
    /**
     * This method sets both variables (origin and destination) in a separated way
     */
    void SetVariables(
        TVarType OriginVariable,
        TVarType DestinyVariable
        )
    {
        mOriginVariable = OriginVariable;
        mDestinyVariable = DestinyVariable;
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SimpleMortarMapperProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SimpleMortarMapperProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{
    
    ModelPart& mrThisModelPart; // The model part to compute
    TVarType mOriginVariable;   // The origin variable to map
    TVarType mDestinyVariable;  // The destiny variable to map

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    
    /// Assignment operator.
    SimpleMortarMapperProcess& operator=(SimpleMortarMapperProcess const& rOther) = delete;

    /// Copy constructor.
    //SimpleMortarMapperProcess(SimpleMortarMapperProcess const& rOther);
    
    ///@}
};// class SimpleMortarMapperProcess

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
#endif /* KRATOS_SIMPLE_MORTAR_MAPPER_PROCESS defined */
