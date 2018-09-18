/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <random>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Mathematics/Filters/createFilter.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"

#include "Tudat/Astrodynamics/SystemModels/controlSystem.h"
#include "Tudat/Astrodynamics/SystemModels/guidanceSystem.h"
#include "Tudat/Astrodynamics/SystemModels/navigationInstrumentsModel.h"
#include "Tudat/Astrodynamics/SystemModels/navigationSystem.h"
#include "Tudat/Astrodynamics/SystemModels/onboardComputerModel.h"

//! Get path for output directory.
static inline std::string getOutputPath( const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) - std::string( "thesisOpenLoop.cpp" ).length( ) );
    std::string outputPath = reducedPath + "SimulationOutputOpenLoop/";
    if ( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if ( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

//! Class to keep attitude constant.
class AerobrakingAerodynamicGuidance: public tudat::aerodynamics::AerodynamicGuidance
{
public:

    //! Constructor
    AerobrakingAerodynamicGuidance( ) { }

    //! Update guidance with aerodynamic angles set to zero.
    void updateGuidance( const double time )
    {
        TUDAT_UNUSED_PARAMETER( time );
        currentAngleOfAttack_ = 0.0;
        currentAngleOfSideslip_ = 0.0;
        currentBankAngle_ = -tudat::mathematical_constants::PI;
    }

};

//! Main code.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::aerodynamics;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::ephemerides;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::filters;
    using namespace tudat::gravitation;
    using namespace tudat::input_output;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::propagators;
    using namespace tudat::simulation_setup;
    using namespace tudat::system_models;
    using namespace tudat::unit_conversions;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings
    const double simulationStartEpoch = 7.0 * physical_constants::JULIAN_YEAR + 30.0 * 6.0 * physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = simulationStartEpoch + 4.5 * physical_constants::JULIAN_DAY;

    // Define body settings for simulation
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );

    // Tabulated atmosphere settings
    std::map< int, std::string > tabulatedAtmosphereFiles;
    tabulatedAtmosphereFiles[ 0 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/density.dat";
    tabulatedAtmosphereFiles[ 1 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/pressure.dat";
    tabulatedAtmosphereFiles[ 2 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/temperature.dat";
    tabulatedAtmosphereFiles[ 3 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/gasConstant.dat";
    tabulatedAtmosphereFiles[ 4 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat";
    std::vector< AtmosphereDependentVariables > atmosphereDependentVariables = {
        density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere,
        gas_constant_dependent_atmosphere, specific_heat_ratio_dependent_atmosphere };
    std::vector< AtmosphereIndependentVariables > atmosphereIndependentVariables = {
        longitude_dependent_atmosphere, latitude_dependent_atmosphere, altitude_dependent_atmosphere };
    std::vector< interpolators::BoundaryInterpolationType > boundaryConditions =
            std::vector< interpolators::BoundaryInterpolationType >( 3, interpolators::use_boundary_value );

    // Create body settings and give default values
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    // Give Mars a more detailed environment
    std::vector< double > vectorOfAtmosphereParameters = { 115.0e3, 2.424e-08, 6533.0, -1.0, 0.0, 0.0 };
    bodySettings[ "Mars" ]->gravityFieldSettings = boost::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
    bodySettings[ "Mars" ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >(
                tabulatedAtmosphereFiles, atmosphereIndependentVariables, atmosphereDependentVariables, boundaryConditions );
//    boost::make_shared< CustomConstantTemperatureAtmosphereSettings >(
//                    three_term_atmosphere_model, 215.0, 197.0, 1.3, vectorOfAtmosphereParameters );

    // Give Earth zero gravity field such that ephemeris is created, but no acceleration
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< CentralGravityFieldSettings >( 0.0 );

    // Create body objects
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Define simulation objects
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Satellite" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Mars" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE INITIAL CONDITIONS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Get and set Mars physical parameters
    const double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    const double marsPropagationAtmosphericInterfaceAltitude = 500.0e3;

    // Set initial Keplerian elements for satellite
    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 25957.742944018e3;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.864414391742801;
    initialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 158.7 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 43.6 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 180.0 );

    // Define initial translational state
    const Eigen::Vector6d translationalInitialState = convertKeplerianToCartesianElements( initialStateInKeplerianElements,
                                                                                           marsGravitationalParameter );

    // Simulation times
    const double simulationConstantStepSize = 0.1; // 10 Hz
    const double simulationConstantStepSizeDuringAtmosphericPhase = 0.002; // 500 Hz

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            DEFINE ONBOARD MODEL          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define spacecraft physical characteristics
    const double spacecraftMass = 1000.0;
    Eigen::Matrix3d spacecraftInertiaTensor = Eigen::Matrix3d::Zero( );
    spacecraftInertiaTensor( 0, 0 ) = 5750.0;
    spacecraftInertiaTensor( 1, 1 ) = 1215.0;
    spacecraftInertiaTensor( 2, 2 ) = 5210.0;
    const double referenceAreaAerodynamic = 37.5;
    const double referenceAreaRadiation = 37.5;
    const double radiationPressureCoefficient = 1.0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object
    bodyMap[ "Satellite" ] = boost::make_shared< Body >( );
    bodyMap[ "Satellite" ]->setConstantBodyMass( spacecraftMass );
    bodyMap[ "Satellite" ]->setBodyInertiaTensor( spacecraftInertiaTensor );

    // Aerodynamic coefficients from file
    std::map< int, std::string > aerodynamicForceCoefficientFiles;
    aerodynamicForceCoefficientFiles[ 0 ] = getTudatRootPath( ) + "External/MRODragCoefficients.txt";
    aerodynamicForceCoefficientFiles[ 2 ] = getTudatRootPath( ) + "External/MROLiftCoefficients.txt";

    // Create aerodynamic coefficient settings
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            readTabulatedAerodynamicCoefficientsFromFiles(
                aerodynamicForceCoefficientFiles, referenceAreaAerodynamic,
                boost::assign::list_of( angle_of_attack_dependent )( altitude_dependent ), true, true );

    // Constant radiation pressure variables
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    boost::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Set aerodynamic coefficient and radiation pressure settings
    bodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface( createAerodynamicCoefficientInterface(
                                                                    aerodynamicCoefficientSettings, "Satellite" ) );
    bodyMap[ "Satellite" ]->setRadiationPressureInterface( "Sun", createRadiationPressureInterface(
                                                               radiationPressureSettings, "Satellite", bodyMap ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define acceleration settings
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    accelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 21, 21 ) );
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        if ( bodiesToCreate.at( i ) != "Mars" )
        {
            accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
    }
    accelerationsOfSatellite[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
    accelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    // Set accelerations settings
    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE DYNAMICS SETTINGS               ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Add aerodynamic guidance
    boost::shared_ptr< AerodynamicGuidance > aerodynamicGuidance = boost::make_shared< AerobrakingAerodynamicGuidance >( );
    setGuidanceAnglesFunctions( aerodynamicGuidance, bodyMap.at( "Satellite" ) );

    // Define termination conditions
    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Satellite", "Mars" );
    terminationSettingsList.push_back( boost::make_shared< PropagationDependentVariableTerminationSettings >(
                                           terminationDependentVariable, 0.95 * marsPropagationAtmosphericInterfaceAltitude, true ) );

    terminationSettingsList.push_back( boost::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );

    // Create integrator settings
//    const unsigned int saveFrequency = static_cast< unsigned int >( onboardComputerRefreshRate );
//    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//                simulationStartEpoch, simulationConstantStepSize, RungeKuttaCoefficients::rungeKuttaFehlberg56,
//                1.0e-5, 5.0, 1.0e-14, 1.0e-14, saveFrequency, false );
    const unsigned int saveFrequency = 1;//10 * static_cast< unsigned int >( onboardComputerRefreshRate );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< IntegratorSettings< > >(
                rungeKutta4, simulationStartEpoch, simulationConstantStepSize, saveFrequency, false );

    // Dependent variables
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                          total_acceleration_dependent_variable, "Satellite", "Mars" ) );
    dependentVariablesList.push_back( boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                          aerodynamic, "Satellite", "Mars", false ) );
    dependentVariablesList.push_back( boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Satellite", reference_frames::angle_of_attack ) );
    dependentVariablesList.push_back( boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Satellite", reference_frames::angle_of_sideslip ) );
    dependentVariablesList.push_back( boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Satellite", reference_frames::bank_angle ) );
    dependentVariablesList.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                          local_density_dependent_variable, "Satellite", "Mars" ) );
    dependentVariablesList.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                          local_dynamic_pressure_dependent_variable, "Satellite", "Mars" ) );
    dependentVariablesList.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                          local_aerodynamic_heat_rate_dependent_variable, "Satellite", "Mars" ) );

    // Create propagation settings for translational dynamics
    boost::shared_ptr< TranslationalStatePropagatorSettings< > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, translationalInitialState,
                terminationSettings, cowell );

    // Set full propagation settings
    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< > > > propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    boost::shared_ptr< PropagatorSettings< > > propagatorSettings =
            boost::make_shared< MultiTypePropagatorSettings< > >(
                propagatorSettingsList, terminationSettings,
                boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create empty dynamics simulator object
    boost::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator;

    // Pre-allocate variables
    std::map< double, Eigen::VectorXd > fullIntegrationResult;
    std::map< double, Eigen::VectorXd > fullDependentVariablesResults;
    std::map< double, Eigen::VectorXd > cartesianTranslationalIntegrationResult;
    std::map< double, Eigen::VectorXd > keplerianTranslationalIntegrationResult;
    std::map< double, Eigen::VectorXd > dependentVariablesResults;

    // Create loop for each orbit and check condition for stopping aerobraking at the end
    std::vector< bool > vectorOfTerminationConditions;
    do
    {
        // Propagate until apoapsis
        dynamicsSimulator = boost::make_shared< SingleArcDynamicsSimulator< > >(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false, false );

        // Retrieve elements from dynamics simulator
        std::map< double, Eigen::VectorXd > currentFullIntegrationResult = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > currentDependentVariablesResults = dynamicsSimulator->getDependentVariableHistory( );
        Eigen::VectorXd finalPropagatedState = currentFullIntegrationResult.rbegin( )->second;

        // Reset time
        integratorSettings->initialTime_ = currentFullIntegrationResult.rbegin( )->first;

        // Reset other conditions based on termination reason
        if ( dynamicsSimulator->getPropagationTerminationReason( )->getPropagationTerminationReason( ) == termination_condition_reached )
        {
            // Check which condition was met
            boost::shared_ptr< PropagationTerminationDetailsFromHybridCondition > hybridTerminationDetails =
                    boost::dynamic_pointer_cast< PropagationTerminationDetailsFromHybridCondition >(
                        dynamicsSimulator->getPropagationTerminationReason( ) );
            vectorOfTerminationConditions = hybridTerminationDetails->getWasConditionMetWhenStopping( );
            if ( vectorOfTerminationConditions.at( 0 ) ) // altitude
            {
                // Extract termination object
                boost::shared_ptr< PropagationDependentVariableTerminationSettings > altitudeTerminationSettings =
                        boost::dynamic_pointer_cast< PropagationDependentVariableTerminationSettings >( terminationSettingsList.at( 0 ) );

                // Reset integration time step
                if ( altitudeTerminationSettings->useAsLowerLimit_ )
                {
                    integratorSettings->initialTimeStep_ = simulationConstantStepSizeDuringAtmosphericPhase;
                    altitudeTerminationSettings->limitValue_ = 1.05 * marsPropagationAtmosphericInterfaceAltitude;
                }
                else
                {
                    integratorSettings->initialTimeStep_ = simulationConstantStepSize;
                    altitudeTerminationSettings->limitValue_ = 0.95 * marsPropagationAtmosphericInterfaceAltitude;
                }

                // Invert flag for next propagation
                altitudeTerminationSettings->useAsLowerLimit_ = !altitudeTerminationSettings->useAsLowerLimit_;

                // Overwrite propagation termination settings
                terminationSettingsList.at( 0 ) = altitudeTerminationSettings;
                terminationSettings = boost::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< > >(
                            propagatorSettings )->resetTerminationSettings( terminationSettings );
            }
        }
        else
        {
            // Inform user on error and break propagation
            std::cerr << "Propagation was interrupted prematurely." << std::endl;

            // Save results
            fullIntegrationResult.insert( currentFullIntegrationResult.begin( ), currentFullIntegrationResult.end( ) );
            fullDependentVariablesResults.insert( currentDependentVariablesResults.begin( ), currentDependentVariablesResults.end( ) );
            break;
        }

        // Reset state
        propagatorSettings->resetInitialStates( finalPropagatedState );

        // Save results
        fullIntegrationResult.insert( currentFullIntegrationResult.begin( ), currentFullIntegrationResult.end( ) );
        fullDependentVariablesResults.insert( currentDependentVariablesResults.begin( ), currentDependentVariablesResults.end( ) );
    }
    while ( !vectorOfTerminationConditions.at( 1 ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             RETRIEVE RESULTS           ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Compute map of Kepler elements
    Eigen::VectorXd currentFullState;
    Eigen::Vector6d currentCartesianState;
    for ( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = fullIntegrationResult.begin( );
          stateIterator != fullIntegrationResult.end( ); stateIterator++ )
    {
        // Get current states
        currentFullState = stateIterator->second;
        currentCartesianState = currentFullState.segment( 0, 6 );

        // Store translational and rotational states
        cartesianTranslationalIntegrationResult[ stateIterator->first ] = currentCartesianState;

        // Compute current Keplerian state
        keplerianTranslationalIntegrationResult[ stateIterator->first ] = convertCartesianToKeplerianElements(
                    currentCartesianState, marsGravitationalParameter );

        // Store dependent variables
        dependentVariablesResults[ stateIterator->first ] = fullDependentVariablesResults[ stateIterator->first ];
    }

    // Write propagation history to files
    writeDataMapToTextFile( cartesianTranslationalIntegrationResult, "cartesianPropagated.dat", getOutputPath( ) );
    writeDataMapToTextFile( keplerianTranslationalIntegrationResult, "keplerianPropagated.dat", getOutputPath( ) );
    writeDataMapToTextFile( dependentVariablesResults, "dependentVariables.dat", getOutputPath( ) );

    // Final statement
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed
    return EXIT_SUCCESS;
}
