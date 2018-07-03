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

#include "Tudat/Astrodynamics/GuidanceNavigationControl/controlSystem.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/guidanceSystem.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigationSystem.h"

#include "Tudat/Astrodynamics/SystemModels/navigationInstrumentsModel.h"
#include "Tudat/Astrodynamics/SystemModels/onboardComputerModel.h"

//! Get path for output directory.
static inline std::string getOutputPath( const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) - std::string( "thesis.cpp" ).length( ) );
    std::string outputPath = reducedPath + "SimulationOutput/";
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
    using namespace tudat::guidance_navigation_control;
    using namespace tudat::input_output;
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
    const double simulationStartEpoch = 7.0 * physical_constants::JULIAN_YEAR +
            30.0 * 6.0 * physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = //600.0 + simulationStartEpoch;
            1.0 * physical_constants::JULIAN_DAY + simulationStartEpoch;

    // Define body settings for simulation
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
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
    std::vector< interpolators::BoundaryInterpolationType > boundaryConditions = {
        interpolators::use_boundary_value, interpolators::use_boundary_value, interpolators::use_default_value };

    // Define default extrapolation values
    std::vector< std::vector< std::pair< double, double > > > extrapolationValues =
            std::vector< std::vector< std::pair< double, double > > >(
                5, std::vector< std::pair< double, double > >( 3, std::make_pair( 0.0, 0.0 ) ) );
    extrapolationValues.at( 0 ).at( 2 ) = { 7.62e-5, 8.0e-18 };
    extrapolationValues.at( 1 ).at( 2 ) = { 2.35, 1.45e-11 };
    extrapolationValues.at( 2 ).at( 2 ) = { 1.61e2, 186.813 };
    extrapolationValues.at( 3 ).at( 2 ) = { 190.7, 8183.0 };
    extrapolationValues.at( 4 ).at( 2 ) = { 1.377, 1.667 };

    // Create body objects
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    bodySettings[ "Mars" ]->gravityFieldSettings = boost::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
    bodySettings[ "Mars" ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >(
                tabulatedAtmosphereFiles, atmosphereIndependentVariables, atmosphereDependentVariables,
                boundaryConditions, extrapolationValues );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object
    bodyMap[ "Satellite" ] = boost::make_shared< Body >( );
    const double spacecraftMass = 1000.0;
    bodyMap[ "Satellite" ]->setConstantBodyMass( spacecraftMass );

    Eigen::Matrix3d spacecraftInertiaTensor = Eigen::Matrix3d::Zero( );
    spacecraftInertiaTensor( 0, 0 ) = 5750.0;
    spacecraftInertiaTensor( 1, 1 ) = 1215.0;
    spacecraftInertiaTensor( 2, 2 ) = 5210.0;
    bodyMap[ "Satellite" ]->setBodyInertiaTensor( spacecraftInertiaTensor );

    // Aerodynamic coefficients from file
    std::map< int, std::string > aerodynamicForceCoefficientFiles;
    std::map< int, std::string > aerodynamicMomentCoefficientFiles;
    aerodynamicForceCoefficientFiles[ 0 ] = getTudatRootPath( ) + "External/MRODragCoefficients.txt";
    aerodynamicForceCoefficientFiles[ 2 ] = getTudatRootPath( ) + "External/MROLiftCoefficients.txt";
    aerodynamicMomentCoefficientFiles[ 1 ] = getTudatRootPath( ) + "External/MROMomentCoefficients.txt";

    // Create aerodynamic coefficient settings
    const double referenceAreaAerodynamic = 37.5;
    const double referenceLengthAerodynamic = 2.5;
    const Eigen::Vector3d momentReferencePoint = Eigen::Vector3d::Zero( );
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            readTabulatedAerodynamicCoefficientsFromFiles(
                aerodynamicForceCoefficientFiles, aerodynamicMomentCoefficientFiles, referenceLengthAerodynamic,
                referenceAreaAerodynamic, referenceLengthAerodynamic, momentReferencePoint,
                boost::assign::list_of( aerodynamics::angle_of_attack_dependent )( aerodynamics::altitude_dependent ),
                true, true );

    // Constant radiation pressure variables
    const double referenceAreaRadiation = 37.5;
    const double radiationPressureCoefficient = 1.0;
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

    // Define simulation objects
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Satellite" );
    centralBodies.push_back( "Mars" );

    // Define acceleration settings
    SelectedAccelerationMap accelerationMap;
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
    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    // Define and set torque settings
    SelectedTorqueMap torqueMap;
    torqueMap[ "Satellite" ][ "Mars" ].push_back( boost::make_shared< TorqueSettings >( aerodynamic_torque ) );
    torqueMap[ "Satellite" ][ "Mars" ].push_back( boost::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
    TorqueModelMap torqueModelMap = createTorqueModelsMap( bodyMap, torqueMap );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE INITIAL CONDITIONS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Get and set Mars physical parameters
    double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    double marsRadius = bodyMap.at( "Mars" )->getShapeModel( )->getAverageRadius( );
    double marsAtmosphericInterfaceAltitude = 500.0e3;

    // Set initial Keplerian elements for satellite
    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 27228500;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.869218;
    initialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 158.7 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 43.6 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 180.0 );

    const Eigen::Vector6d translationalInitialState = convertKeplerianToCartesianElements(
                initialStateInKeplerianElements, marsGravitationalParameter );

    // Define initial rotational state
    Eigen::Quaterniond initialStateInQuaternionElements = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    Eigen::Vector7d rotationalInitialState = Eigen::Vector7d::Zero( );
    rotationalInitialState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( initialStateInQuaternionElements );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE INSTRUMENT MODELS               ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define map of central bodies for onboard instrument system
    std::map< std::string, std::string > centralBodyMap;
    for ( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        centralBodyMap[ bodiesToPropagate.at( i ) ] = centralBodies.at( i );
    }

    // Define error characteristics of instruments
    double simulationConstantStepSize = 1.0;
    double onboardComputerRefreshStepSize = simulationConstantStepSize; // seconds
    double onboardComputerRefreshRate = 1.0 / onboardComputerRefreshStepSize; // Hertz
    Eigen::Vector3d accelerationBias, accelerometerScaleFactor, gyroscopeBias, gyroscopeScaleFactor;
    Eigen::Vector6d accelerometerMisalignment, gyroscopeMisalignment;

    double gyroscopeBiasStandardDeviation = 5.0e-9;
    double gyroscopeScaleFactorStandardDeviation = 1e-4;

    std::default_random_engine generator;
    std::normal_distribution< double > distributionOne( 0.0, gyroscopeBiasStandardDeviation );
    std::normal_distribution< double > distributionTwo( 0.0, gyroscopeScaleFactorStandardDeviation );
    std::normal_distribution< double > distributionThree( 0.0, 1.0e-3 );
    for ( unsigned int i = 0; i < 6; i++ )
    {
        if ( i < 3 )
        {
            accelerationBias[ i ] = distributionTwo( generator );
            accelerometerScaleFactor[ i ] = distributionTwo( generator );
            gyroscopeBias[ i ] = distributionOne( generator );
            gyroscopeScaleFactor[ i ] = distributionTwo( generator );
        }
        accelerometerMisalignment[ i ] = distributionThree( generator );
        gyroscopeMisalignment[ i ] = distributionThree( generator );
    }

    Eigen::Vector3d accelerometerAccuracy = Eigen::Vector3d::Constant( 2.0e-4 / std::sqrt( onboardComputerRefreshRate ) );
    Eigen::Vector3d gyroscopeAccuracy = Eigen::Vector3d::Constant( 3.0e-7 / std::sqrt( onboardComputerRefreshRate ) );
    Eigen::Vector3d starTrackerAccuracy = Eigen::Vector3d::Constant( 20.0 / 3600.0 );

    // Create navigation instruments object
    boost::shared_ptr< NavigationInstrumentsModel > onboardInstruments =
            boost::make_shared< NavigationInstrumentsModel >( bodyMap, accelerationModelMap, "Satellite" );

    // Add inertial measurement unit
    onboardInstruments->addInertialMeasurementUnit( accelerationBias, accelerometerScaleFactor,
                                                    accelerometerMisalignment, accelerometerAccuracy,
                                                    gyroscopeBias, gyroscopeScaleFactor,
                                                    gyroscopeMisalignment, gyroscopeAccuracy );

    // Ad star tracker
    onboardInstruments->addStarTracker( 2, starTrackerAccuracy );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             ONBOARD PARAMETERS                     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Initial conditions
    Eigen::Vector16d initialEstimatedStateVector;
    initialEstimatedStateVector << translationalInitialState, rotationalInitialState.segment( 0, 4 ), // assume perfect initial knowledge
            Eigen::Vector6d::Zero( ); // assume zero errors in gyroscope
    Eigen::Matrix16d initialEstimatedStateCovarianceMatrix = Eigen::Matrix16d::Identity( 16, 16 );

    // System and measurment uncertainties
    double positionStandardDeviation = 10.0;
    double translationalVelocityStandardDeviation = 1.0e-2;
    double attitudeStandardDeviation = 1.0e-6;
    Eigen::Vector16d diagonalOfSystemUncertainty;
    diagonalOfSystemUncertainty << Eigen::Vector3d::Constant( std::pow( positionStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( translationalVelocityStandardDeviation, 2 ) ),
            Eigen::Vector4d::Constant( std::pow( attitudeStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( gyroscopeBiasStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( gyroscopeScaleFactorStandardDeviation, 2 ) );
    Eigen::Matrix16d systemUncertainty = diagonalOfSystemUncertainty.asDiagonal( );

    Eigen::Vector7d diagonalOfMeasurementUncertainty;
    diagonalOfMeasurementUncertainty << gyroscopeAccuracy, Eigen::Vector4d::Constant( 1e-6 ); // starTrackerAccuracy <<<<----
    Eigen::Matrix< double, 7, 7 > measurementUncertainty = diagonalOfMeasurementUncertainty.asDiagonal( );

    // Aerodynamic coefficients
    Eigen::Vector3d onboardAerodynamicCoefficients = Eigen::Vector3d::Zero( );
    onboardAerodynamicCoefficients[ 0 ] = 1.87; // drag coefficient at 0 deg and 125 km
    onboardAerodynamicCoefficients[ 2 ] = 0.013; // lift coefficient at 0 deg and 125 km

    // Controller gains
    Eigen::Vector3d proportionalGain = Eigen::Vector3d::Constant( 0.0 );
    Eigen::Vector3d integralGain = Eigen::Vector3d::Constant( 0.0 );
    Eigen::Vector3d derivativeGain = Eigen::Vector3d::Constant( 0.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            DEFINE ONBOARD MODEL          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create body map for onboard model
    std::map< std::string, boost::shared_ptr< BodySettings > > onboardBodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        onboardBodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        onboardBodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    onboardBodySettings[ "Mars" ]->gravityFieldSettings =
            boost::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
    onboardBodySettings[ "Mars" ]->atmosphereSettings = boost::make_shared< ExponentialAtmosphereSettings >( aerodynamics::mars );
    NamedBodyMap onboardBodyMap = createBodies( onboardBodySettings );
    onboardBodyMap[ "Satellite" ] = boost::make_shared< Body >( );
    onboardBodyMap[ "Satellite" ]->setConstantBodyMass( spacecraftMass );
    onboardBodyMap[ "Satellite" ]->setBodyInertiaTensor( spacecraftInertiaTensor );
    onboardBodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                                                           referenceAreaAerodynamic, onboardAerodynamicCoefficients, true, true ),
                                                       "Satellite" ) );
    setGlobalFrameBodyEphemerides( onboardBodyMap, "SSB", "J2000" );

    // Define acceleration settings for onboard model
    SelectedAccelerationMap onboardAccelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > onboardAccelerationsOfSatellite;
    onboardAccelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 4, 4 ) );
    onboardAccelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    // Set accelerations settings
    onboardAccelerationMap[ "Satellite" ] = onboardAccelerationsOfSatellite;
    AccelerationMap onboardAccelerationModelMap = createAccelerationModelsMap(
                onboardBodyMap, onboardAccelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE GNC MODELS                      ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create control system object
    boost::shared_ptr< ControlSystem > controlSystem = boost::make_shared< ControlSystem >(
                proportionalGain, integralGain, derivativeGain );

    // Create unscented Kalman filter settings object for navigation
    boost::shared_ptr< IntegratorSettings< > > filterIntegratorSettings =
            boost::make_shared< IntegratorSettings< > >( euler, simulationStartEpoch, onboardComputerRefreshStepSize );
    boost::shared_ptr< FilterSettings< > > filteringSettings = boost::make_shared< UnscentedKalmanFilterSettings< > >(
                systemUncertainty, measurementUncertainty, simulationStartEpoch, initialEstimatedStateVector,
                initialEstimatedStateCovarianceMatrix, filterIntegratorSettings );

    // Create giudance system object
    boost::shared_ptr< GuidanceSystem > guidanceSystem = boost::make_shared< GuidanceSystem >( );

    // Create navigation system object
    boost::shared_ptr< NavigationSystem > navigationSystem = boost::make_shared< NavigationSystem >(
                onboardBodyMap, onboardAccelerationModelMap, "Satellite", "Mars",
                filteringSettings, three_term_atmosphere_model, marsAtmosphericInterfaceAltitude );

    // Create onboard computer object
    boost::shared_ptr< OnboardComputerModel > onboardComputer = boost::make_shared< OnboardComputerModel >(
                controlSystem, guidanceSystem, navigationSystem, onboardInstruments );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE DYNAMICS SETTINGS               ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define termination conditions
    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    terminationSettingsList.push_back( boost::make_shared< CustomTerminationSettings >(
                                           boost::bind( &OnboardComputerModel::checkStopCondition, onboardComputer, _1 ) ) );
    terminationSettingsList.push_back( boost::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );

    // Create integrator settings
    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< IntegratorSettings< > >(
                euler, simulationStartEpoch, simulationConstantStepSize );

    // Create propagation settings for translational dynamics
    boost::shared_ptr< TranslationalStatePropagatorSettings< > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, translationalInitialState,
                terminationSettings, cowell );

    // Create propagator settings for rotation
    boost::shared_ptr< RotationalStatePropagatorSettings< > > rotationalPropagatorSettings =
            boost::make_shared< RotationalStatePropagatorSettings< > >(
                torqueModelMap, bodiesToPropagate, rotationalInitialState, terminationSettings, quaternions );

    // Set full propagation settings
    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< > > > propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    propagatorSettingsList.push_back( rotationalPropagatorSettings );
    boost::shared_ptr< PropagatorSettings< > > propagatorSettings =
            boost::make_shared< MultiTypePropagatorSettings< > >( propagatorSettingsList, terminationSettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false, true );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             RETRIEVE RESULTS           ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Get actual states
    std::map< double, Eigen::VectorXd > fullIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > cartesianTranslationalIntegrationResult;
    std::map< double, Eigen::VectorXd > keplerianTranslationalIntegrationResult;
    std::map< double, Eigen::VectorXd > quaternionsRotationalIntegrationResult;

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
        quaternionsRotationalIntegrationResult[ stateIterator->first ] = currentFullState.segment( 6, 7 );

        // Copute current Keplerian state
        keplerianTranslationalIntegrationResult[ stateIterator->first ] = convertCartesianToKeplerianElements( currentCartesianState,
                                                                                                               marsGravitationalParameter );
    }

    // Get estimated states
    std::map< double, std::pair< std::pair< Eigen::Vector6d, Eigen::Vector6d >, Eigen::Vector7d > > fullEstimationResult =
            navigationSystem->getHistoryOfEstimatedStates( );
    std::map< double, Eigen::VectorXd > cartesianTranslationalEstimationResult;
    std::map< double, Eigen::VectorXd > keplerianTranslationalEstimationResult;
    std::map< double, Eigen::VectorXd > quaternionsRotationalEstimationResult;

    // Extract map of translational and rotational results
    for ( std::map< double, std::pair< std::pair< Eigen::Vector6d, Eigen::Vector6d >, Eigen::Vector7d > >::const_iterator
          stateIterator = fullEstimationResult.begin( ); stateIterator != fullEstimationResult.end( ); stateIterator++ )
    {
        // Add results to maps
        cartesianTranslationalEstimationResult[ stateIterator->first ] = stateIterator->second.first.first;
        keplerianTranslationalEstimationResult[ stateIterator->first ] = stateIterator->second.first.second;
        quaternionsRotationalEstimationResult[ stateIterator->first ] = stateIterator->second.second;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Write perturbed satellite propagation history to files
    writeDataMapToTextFile( cartesianTranslationalIntegrationResult, "cartesianPropagated.dat", getOutputPath( ) );
    writeDataMapToTextFile( keplerianTranslationalIntegrationResult, "keplerianPropagated.dat", getOutputPath( ) );
    writeDataMapToTextFile( quaternionsRotationalIntegrationResult, "rotationalPropagated.dat", getOutputPath( ) );

    // Write estimated satellite state hisotry to files
    writeDataMapToTextFile( cartesianTranslationalEstimationResult, "cartesianEstimated.dat", getOutputPath( ) );
    writeDataMapToTextFile( keplerianTranslationalEstimationResult, "keplerianEstimated.dat", getOutputPath( ) );
    writeDataMapToTextFile( quaternionsRotationalEstimationResult, "rotationalEstimated.dat", getOutputPath( ) );

    // Final statement
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed
    return EXIT_SUCCESS;
}
