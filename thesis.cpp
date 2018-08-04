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

//! Function to compute the altimeter error as a function of altitude.
double altimeterErrorAsAFunctionOfAltitude( const double currentPseudoAltitudeMeasurement )
{
    TUDAT_UNUSED_PARAMETER( currentPseudoAltitudeMeasurement );
    return 1.0e-3;//0.001 * currentPseudoAltitudeMeasurement;
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
            1.4 * physical_constants::JULIAN_DAY + simulationStartEpoch;

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
    std::vector< interpolators::BoundaryInterpolationType > boundaryConditions =
            std::vector< interpolators::BoundaryInterpolationType >( 3, interpolators::use_boundary_value );

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
                tabulatedAtmosphereFiles, atmosphereIndependentVariables, atmosphereDependentVariables, boundaryConditions );
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
    const double marsAtmosphericInterfaceAltitude = 250.0e3;
    const double marsReducedAtmosphericInterfaceAltitude = 175.0e3;
    const unsigned int numberOfRequiredAtmosphereSamplesForInitiation = 7;

    // Set initial Keplerian elements for satellite
    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 25956000.0;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.864540;
    initialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 158.7 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 43.6 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 180.0 );

    // Define initial translational state
    const Eigen::Vector6d translationalInitialState = convertKeplerianToCartesianElements( initialStateInKeplerianElements,
                                                                                           marsGravitationalParameter );

    // Find quaternion of initial rotation
    Eigen::Matrix3d directionCosineMatrix = Eigen::Matrix3d::Zero( );

    // Find the trajectory x-axis unit vector
    Eigen::Vector3d initialRadialVector = translationalInitialState.segment( 0, 3 );
    Eigen::Vector3d xUnitVector = ( translationalInitialState.segment( 3, 3 ) -
                                    7.07763225880808e-05 * Eigen::Vector3d::UnitZ( ).cross( initialRadialVector ) ).normalized( );
    directionCosineMatrix.col( 0 ) = xUnitVector;

    // Find trajectory z-axis unit vector
    Eigen::Vector3d zUnitVector = initialRadialVector.normalized( );
    zUnitVector -= zUnitVector.dot( xUnitVector ) * xUnitVector;
    directionCosineMatrix.col( 2 ) = zUnitVector;

    // Find body-fixed y-axis unit vector
    directionCosineMatrix.col( 1 ) = zUnitVector.cross( xUnitVector );

    // Define initial rotational state
    Eigen::Vector7d rotationalInitialState = Eigen::Vector7d::Zero( );
    Eigen::Vector4d initialQuaternion = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( directionCosineMatrix ) );
    rotationalInitialState.segment( 0, 4 ) = initialQuaternion;

    // Simulation times
    double simulationConstantStepSize = 0.5;
    double onboardComputerRefreshStepSize = simulationConstantStepSize; // seconds
    double onboardComputerRefreshRate = 1.0 / onboardComputerRefreshStepSize; // Hertz

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE ONBOARD PARAMETERS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Initial conditions
    Eigen::Vector22d initialEstimatedStateVector = Eigen::Vector22d::Zero( ); // assume zero errors in IMU
    initialEstimatedStateVector.segment( 0, 6 ) = translationalInitialState;
    initialEstimatedStateVector.segment( 6, 4 ) = rotationalInitialState.segment( 0, 4 ); // assume perfect initial knowledge
    Eigen::Matrix22d initialEstimatedStateCovarianceMatrix = Eigen::Matrix22d::Identity( );

    // Altimeter specifications
    Eigen::Vector3d altimeterBodyFixedPointingDirection = -Eigen::Vector3d::UnitZ( );
    std::pair< double, double > altimeterAltitudeRange = std::make_pair( 0.0, 5e6 );

    // Deifine instrument accuracy
    double accelerometerBiasStandardDeviation = 1e-4;
    double accelerometerScaleFactorStandardDeviation = 1e-4;
    double gyroscopeBiasStandardDeviation = 5.0e-9;
    double gyroscopeScaleFactorStandardDeviation = 1e-4;
    Eigen::Vector3d accelerometerAccuracy = Eigen::Vector3d::Constant( 2.0e-4 * std::sqrt( onboardComputerRefreshRate ) );
    Eigen::Vector3d gyroscopeAccuracy = Eigen::Vector3d::Constant( 3.0e-7 * std::sqrt( onboardComputerRefreshRate ) );
    Eigen::Vector3d starTrackerAccuracy = Eigen::Vector3d::Constant( 20.0 / 3600.0 );
    double altimeterAccuracy = 1.0e2;
    boost::function< double( double ) > altimeterAccuracyFunction = boost::bind( &altimeterErrorAsAFunctionOfAltitude, _1 );

    // System and measurment uncertainties
    double positionStandardDeviation = 100.0;
    double translationalVelocityStandardDeviation = 1.0e-1;
    double attitudeStandardDeviation = 1.0e-4;
    Eigen::Vector22d diagonalOfSystemUncertainty;
    diagonalOfSystemUncertainty << Eigen::Vector3d::Constant( std::pow( positionStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( translationalVelocityStandardDeviation, 2 ) ),
            Eigen::Vector4d::Constant( std::pow( attitudeStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( accelerometerBiasStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( accelerometerScaleFactorStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( gyroscopeBiasStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( gyroscopeScaleFactorStandardDeviation, 2 ) );
    Eigen::Matrix22d systemUncertainty = diagonalOfSystemUncertainty.asDiagonal( );

    Eigen::Vector5d diagonalOfMeasurementUncertainty;
    diagonalOfMeasurementUncertainty << altimeterAccuracy, Eigen::Vector4d::Constant( 1e-6 ); // starTrackerAccuracy <<<<----
    Eigen::Matrix5d measurementUncertainty = diagonalOfMeasurementUncertainty.asDiagonal( );

    // Aerodynamic coefficients
    Eigen::Vector3d onboardAerodynamicCoefficients = Eigen::Vector3d::Zero( );
    onboardAerodynamicCoefficients[ 0 ] = 1.87; // drag coefficient at 0 deg and 125 km
    onboardAerodynamicCoefficients[ 2 ] = 0.013; // lift coefficient at 0 deg and 125 km

    // Controller gains
    Eigen::Vector3d proportionalGain = Eigen::Vector3d::Constant( 0.0 );
//    proportionalGain << 1.0e-3, 5.0e-3, 1.0e-3;
    Eigen::Vector3d integralGain = Eigen::Vector3d::Constant( 0.0 );//1.0e-4 );
    Eigen::Vector3d derivativeGain = Eigen::Vector3d::Constant( 0.0 );//1.0e-1 );

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
    const double referenceLengthAerodynamic = 2.5;
    const Eigen::Vector3d momentReferencePoint = Eigen::Vector3d::Zero( );
    const double referenceAreaRadiation = 37.5;
    const double radiationPressureCoefficient = 1.0;

    // Create body map for onboard model
    std::map< std::string, boost::shared_ptr< BodySettings > > onboardBodySettings =
            getDefaultBodySettings( { "Mars" }, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    onboardBodySettings[ "Mars" ]->ephemerisSettings->resetFrameOrientation( "J2000" );
    onboardBodySettings[ "Mars" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    onboardBodySettings[ "Mars" ]->gravityFieldSettings =
            boost::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
    NamedBodyMap onboardBodyMap = createBodies( onboardBodySettings );

    std::vector< double > vectorOfModelSpecificParameters = { 115.0e3, 2.424e-08, 6533.0 };//, -1.0, 0.0, 0.0 };
    onboardBodyMap[ "Mars" ]->setAtmosphereModel(
                boost::make_shared< CustomConstantTemperatureAtmosphere >( exponential_atmosphere_model, 215.0, 197.0, 1.3,
                                                                           vectorOfModelSpecificParameters ) );

    onboardBodyMap[ "Satellite" ] = boost::make_shared< Body >( );
    onboardBodyMap[ "Satellite" ]->setConstantBodyMass( spacecraftMass );
    onboardBodyMap[ "Satellite" ]->setBodyInertiaTensor( spacecraftInertiaTensor );
    onboardBodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                                                           referenceAreaAerodynamic, onboardAerodynamicCoefficients, true, true ),
                                                       "Satellite" ) );

    setGlobalFrameBodyEphemerides( onboardBodyMap, "Mars", "J2000" );

    // Define acceleration settings for onboard model
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > onboardAccelerationsOfSatellite;
    onboardAccelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
    onboardAccelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    // Set accelerations settings
    SelectedAccelerationMap onboardAccelerationMap;
    onboardAccelerationMap[ "Satellite" ] = onboardAccelerationsOfSatellite;
    AccelerationMap onboardAccelerationModelMap = createAccelerationModelsMap(
                onboardBodyMap, onboardAccelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE GNC MODELS                      ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create control system object
    boost::shared_ptr< ControlSystem > controlSystem = boost::make_shared< ControlSystem >(
                proportionalGain, integralGain, derivativeGain );

    // Create giudance system object
    boost::shared_ptr< GuidanceSystem > guidanceSystem = boost::make_shared< GuidanceSystem >( 250e3, 320.0e3, 2800.0, 500.0e3, 0.19 );

    // Create unscented Kalman filter settings object for navigation
    boost::shared_ptr< IntegratorSettings< > > filterIntegratorSettings =
            boost::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, onboardComputerRefreshStepSize );
    boost::shared_ptr< FilterSettings< > > filteringSettings = boost::make_shared< UnscentedKalmanFilterSettings< > >(
                systemUncertainty, measurementUncertainty, simulationStartEpoch, initialEstimatedStateVector,
                initialEstimatedStateCovarianceMatrix, filterIntegratorSettings );

    // Create navigation system object
    boost::shared_ptr< NavigationSystem > navigationSystem = boost::make_shared< NavigationSystem >(
                onboardBodyMap, onboardAccelerationModelMap, "Satellite", "Mars",
                filteringSettings, exponential_atmosphere_model, marsAtmosphericInterfaceAltitude,
                marsReducedAtmosphericInterfaceAltitude, numberOfRequiredAtmosphereSamplesForInitiation,
                altimeterBodyFixedPointingDirection, altimeterAltitudeRange );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object
    bodyMap[ "Satellite" ] = boost::make_shared< Body >( );
    bodyMap[ "Satellite" ]->setConstantBodyMass( spacecraftMass );
    bodyMap[ "Satellite" ]->setBodyInertiaTensor( spacecraftInertiaTensor );

    // Aerodynamic coefficients from file
    std::map< int, std::string > aerodynamicForceCoefficientFiles;
    std::map< int, std::string > aerodynamicMomentCoefficientFiles;
    aerodynamicForceCoefficientFiles[ 0 ] = getTudatRootPath( ) + "External/MRODragCoefficients.txt";
    aerodynamicForceCoefficientFiles[ 2 ] = getTudatRootPath( ) + "External/MROLiftCoefficients.txt";
    aerodynamicMomentCoefficientFiles[ 1 ] = getTudatRootPath( ) + "External/MROMomentCoefficients.txt";

    // Create aerodynamic coefficient settings
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            readTabulatedAerodynamicCoefficientsFromFiles(
                aerodynamicForceCoefficientFiles, aerodynamicMomentCoefficientFiles, referenceLengthAerodynamic,
                referenceAreaAerodynamic, referenceLengthAerodynamic, momentReferencePoint,
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
    bodyMap[ "Satellite" ]->setControlSystem( controlSystem );

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

    // Define and set torque settings
    SelectedTorqueMap torqueMap;
    torqueMap[ "Satellite" ][ "Mars" ].push_back( boost::make_shared< TorqueSettings >( aerodynamic_torque ) );
    torqueMap[ "Satellite" ][ "Mars" ].push_back( boost::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
    torqueMap[ "Satellite" ][ "Satellite" ].push_back( boost::make_shared< TorqueSettings >( control_torque ) );
    TorqueModelMap torqueModelMap = createTorqueModelsMap( bodyMap, torqueMap );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE INSTRUMENT MODELS               ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Instrument errors
    Eigen::Vector3d accelerationBias, accelerometerScaleFactor, gyroscopeBias, gyroscopeScaleFactor;
    Eigen::Vector6d accelerometerMisalignment, gyroscopeMisalignment;

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

    // Create navigation instruments object
    boost::shared_ptr< NavigationInstrumentsModel > onboardInstruments =
            boost::make_shared< NavigationInstrumentsModel >( bodyMap, accelerationModelMap, "Satellite", "Mars" );

    // Add inertial measurement unit
    onboardInstruments->addInertialMeasurementUnit( accelerationBias, accelerometerScaleFactor,
                                                    accelerometerMisalignment, accelerometerAccuracy,
                                                    gyroscopeBias, gyroscopeScaleFactor,
                                                    gyroscopeMisalignment, gyroscopeAccuracy );

    // Add star tracker
    onboardInstruments->addStarTracker( 2, starTrackerAccuracy );

    // Add altimeter
    onboardInstruments->addAltimeter( altimeterBodyFixedPointingDirection, altimeterAltitudeRange, altimeterAccuracyFunction );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE DYNAMICS SETTINGS               ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create onboard computer object
    boost::shared_ptr< OnboardComputerModel > onboardComputer = boost::make_shared< OnboardComputerModel >(
                controlSystem, guidanceSystem, navigationSystem, onboardInstruments );

    // Define termination conditions
    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    terminationSettingsList.push_back( boost::make_shared< PropagationCustomTerminationSettings >(
                                           boost::bind( &OnboardComputerModel::checkStopCondition, onboardComputer, _1 ) ) );
    terminationSettingsList.push_back( boost::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );

    // Create integrator settings
    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< IntegratorSettings< > >(
                rungeKutta4, simulationStartEpoch, simulationConstantStepSize );

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

    // Create propagator settings for rotation
    boost::shared_ptr< RotationalStatePropagatorSettings< > > rotationalPropagatorSettings =
            boost::make_shared< RotationalStatePropagatorSettings< > >(
                torqueModelMap, bodiesToPropagate, rotationalInitialState, terminationSettings, quaternions );

    // Set full propagation settings
    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< > > > propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    propagatorSettingsList.push_back( rotationalPropagatorSettings );
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
    std::map< double, Eigen::VectorXd > cartesianTranslationalIntegrationResult;
    std::map< double, Eigen::VectorXd > keplerianTranslationalIntegrationResult;
    std::map< double, Eigen::VectorXd > quaternionsRotationalIntegrationResult;
    std::map< double, Eigen::VectorXd > dependentVariablesResults;

    // Create loop for each orbit and check condition for stopping aerobraking at the end
    do
    {
        // Propagate until apoapsis
        dynamicsSimulator = boost::make_shared< SingleArcDynamicsSimulator< > >(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false, true );

        // Retrieve elements from dynamics simulator
        std::map< double, Eigen::VectorXd > currentFullIntegrationResult = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > currentDependentVariablesResults = dynamicsSimulator->getDependentVariableHistory( );

        Eigen::VectorXd finalPropagatedState = currentFullIntegrationResult.rbegin( )->second;

        // Add estimated apoapsis maneuver to state and modify state
        integratorSettings->initialTime_ = currentFullIntegrationResult.rbegin( )->first;
        finalPropagatedState.segment( 3, 3 ) += controlSystem->getScheduledApoapsisManeuver( );
        propagatorSettings->resetInitialStates( finalPropagatedState );

        // Save results
        fullIntegrationResult.insert( currentFullIntegrationResult.begin( ), currentFullIntegrationResult.end( ) );
        dependentVariablesResults.insert( currentDependentVariablesResults.begin( ), currentDependentVariablesResults.end( ) );
    }
    while ( !onboardComputer->isAerobrakingComplete( ) );

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
        quaternionsRotationalIntegrationResult[ stateIterator->first ] = currentFullState.segment( 6, 7 );

        // Copute current Keplerian state
        keplerianTranslationalIntegrationResult[ stateIterator->first ] = convertCartesianToKeplerianElements(
                    currentCartesianState, marsGravitationalParameter );
    }

    // Get estimated states
    std::map< double, std::pair< std::pair< Eigen::Vector6d, Eigen::Vector6d >, Eigen::Vector7d > > fullEstimationResult =
            navigationSystem->getHistoryOfEstimatedStates( );
    std::map< double, Eigen::Vector6d > cartesianTranslationalEstimationResult;
    std::map< double, Eigen::Vector6d > keplerianTranslationalEstimationResult;
    std::map< double, Eigen::Vector7d > quaternionsRotationalEstimationResult;

    // Get control torques
    std::vector< Eigen::Vector3d > vectorOfControlTorques = controlSystem->getCurrentOrbitHistoryOfControlVectors( );
    std::map< double, Eigen::Vector3d > controlTorques;

    // Extract map of translational and rotational results
    unsigned int i = 0;
    for ( std::map< double, std::pair< std::pair< Eigen::Vector6d, Eigen::Vector6d >, Eigen::Vector7d > >::const_iterator
          stateIterator = fullEstimationResult.begin( ); stateIterator != fullEstimationResult.end( ); stateIterator++, i++ )
    {
        // Add results to maps
        cartesianTranslationalEstimationResult[ stateIterator->first ] = stateIterator->second.first.first;
        keplerianTranslationalEstimationResult[ stateIterator->first ] = stateIterator->second.first.second;
        quaternionsRotationalEstimationResult[ stateIterator->first ] = stateIterator->second.second;
        if ( i < vectorOfControlTorques.size( ) )
        {
            controlTorques[ stateIterator->first ] = vectorOfControlTorques.at( i );
        }
    }

    // Get instrument measurements and correct them
//    Eigen::Vector6d estimatedAccelerometerErrors = navigationSystem->getEstimatedAccelerometerErrors( );
    Eigen::Vector6d estimatedAccelerometerErrors = navigationSystem->getCurrentEstimatedState( ).segment( 10, 6 );
    Eigen::Vector6d estimatedGyroscopeErrors = navigationSystem->getCurrentEstimatedState( ).segment( 16, 6 );
    std::map< double, Eigen::Vector6d > inertialMeasurementUnitMeasurements =
            onboardInstruments->getCurrentOrbitHistoryOfInertialMeasurmentUnitMeasurements( );
    for ( std::map< double, Eigen::Vector6d >::iterator
          measurementIterator = inertialMeasurementUnitMeasurements.begin( );
          measurementIterator != inertialMeasurementUnitMeasurements.end( ); measurementIterator++ )
    {
        measurementIterator->second.segment( 0, 3 ) =
                removeErrorsFromInertialMeasurementUnitMeasurement( measurementIterator->second.segment( 0, 3 ),
                                                                    estimatedAccelerometerErrors );
        measurementIterator->second.segment( 3, 3 ) =
                removeErrorsFromInertialMeasurementUnitMeasurement( measurementIterator->second.segment( 3, 3 ),
                                                                    estimatedGyroscopeErrors );
    }
    std::map< double, Eigen::Vector3d > onboardExpectedMeasurements =
            navigationSystem->getCurrentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations( );
    std::cout << "Acc. Error: " << estimatedAccelerometerErrors.transpose( ) << std::endl;
    std::cout << "Gyro. Error: " << estimatedGyroscopeErrors.transpose( ) << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Write perturbed satellite propagation history to files
    writeDataMapToTextFile( cartesianTranslationalIntegrationResult, "cartesianPropagated.dat", getOutputPath( ) );
    writeDataMapToTextFile( keplerianTranslationalIntegrationResult, "keplerianPropagated.dat", getOutputPath( ) );
    writeDataMapToTextFile( quaternionsRotationalIntegrationResult, "rotationalPropagated.dat", getOutputPath( ) );
    writeDataMapToTextFile( dependentVariablesResults, "dependentVariables.dat", getOutputPath( ) );

    // Write estimated satellite state hisotry to files
    writeDataMapToTextFile( cartesianTranslationalEstimationResult, "cartesianEstimated.dat", getOutputPath( ) );
    writeDataMapToTextFile( keplerianTranslationalEstimationResult, "keplerianEstimated.dat", getOutputPath( ) );
    writeDataMapToTextFile( quaternionsRotationalEstimationResult, "rotationalEstimated.dat", getOutputPath( ) );

//    // Write estimated states directly from uscented Kalman filter
//    writeDataMapToTextFile( navigationSystem->getHistoryOfEstimatedStatesFromNavigationFilter( ),
//                            "filterStateEstimates.dat", getOutputPath( ) );
//    writeDataMapToTextFile( navigationSystem->getHistoryOfEstimatedCovarianceFromNavigationFilter( ),
//                            "filterCovarianceEstimates.dat", getOutputPath( ) );

    // Write other data to file
    writeDataMapToTextFile( controlTorques, "controlTorques.dat", getOutputPath( ) );
    writeDataMapToTextFile( inertialMeasurementUnitMeasurements, "inertialMeasurements.dat", getOutputPath( ) );
    writeDataMapToTextFile( onboardExpectedMeasurements, "expectedlMeasurements.dat", getOutputPath( ) );

    // Final statement
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed
    return EXIT_SUCCESS;
}
