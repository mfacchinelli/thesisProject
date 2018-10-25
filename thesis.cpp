/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Facchinelli, M. (2018). Aerobraking Navigation, Guidance and Control.
 *          Master Thesis, Delft University of Technology.
 */

#include <random>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Mathematics/Filters/createFilter.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"

#include "Tudat/Astrodynamics/SystemModels/controlSystem.h"
#include "Tudat/Astrodynamics/SystemModels/guidanceSystem.h"
#include "Tudat/Astrodynamics/SystemModels/instrumentsModel.h"
#include "Tudat/Astrodynamics/SystemModels/navigationSystem.h"
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
    using namespace tudat::input_output;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::propagators;
    using namespace tudat::simulation_setup;
    using namespace tudat::system_models;
    using namespace tudat::unit_conversions;

    // Save settings
    bool extractFilterResults = true;
    bool extractMeasurementResults = false;
    bool extractAtmosphericData = true;
    bool extractManeuverInformation = true;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings
    const double simulationStartEpoch = 7.0 * physical_constants::JULIAN_YEAR + 30.0 * 6.0 * physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = simulationStartEpoch + 1.4 * physical_constants::JULIAN_DAY; // 0.1125

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
    std::vector< AtmosphereIndependentVariables > atmosphereIndependentVariables = {
        longitude_dependent_atmosphere, latitude_dependent_atmosphere, altitude_dependent_atmosphere };
    std::vector< AtmosphereDependentVariables > atmosphereDependentVariables = {
        density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere,
        gas_constant_dependent_atmosphere, specific_heat_ratio_dependent_atmosphere };
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
    bodySettings[ "Mars" ]->gravityFieldSettings = boost::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );

    std::vector< double > vectorOfAtmosphereParameters = { 115.0e3, 2.424e-08, 6533.0, -1.0, 0.0, 0.0 };
    bodySettings[ "Mars" ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >(
                tabulatedAtmosphereFiles, atmosphereIndependentVariables, atmosphereDependentVariables, boundaryConditions );

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
    const double marsPropagationAtmosphericInterfaceAltitude = 1000.0e3;
    const double marsAtmosphericInterfaceAltitude = 200.0e3;
    const double marsReducedAtmosphericInterfaceAltitude = 150.0e3;
    const double periapseEstimatorConstant = 0.955;
    const unsigned int frequencyOfDeepSpaceNetworkTracking = -1;
    const unsigned int numberOfRequiredAtmosphereSamplesForInitiation = 7;

    // Set initial Keplerian elements for satellite
    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 26021000.0;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.859882;
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
    Eigen::Vector3d xUnitVector = - ( translationalInitialState.segment( 3, 3 ) -
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
    const bool useUnscentedKalmanFilter = true;
    const double simulationConstantStepSize = 0.2; // 5 Hz
    const double simulationConstantStepSizeDuringAtmosphericPhase = 0.02; // 50 Hz
    const unsigned int ratioOfOnboardOverSimulatedTimes = 1;
    const double onboardComputerRefreshStepSize = simulationConstantStepSize * ratioOfOnboardOverSimulatedTimes; // seconds
    const double onboardComputerRefreshStepSizeDuringAtmosphericPhase =
            simulationConstantStepSizeDuringAtmosphericPhase * ratioOfOnboardOverSimulatedTimes; // seconds
    const double onboardComputerRefreshRate = 1.0 / onboardComputerRefreshStepSize; // Hertz
    const double onboardComputerRefreshRateDuringAtmosphericPhase = 1.0 / onboardComputerRefreshStepSizeDuringAtmosphericPhase; // Hertz

    // Save frequency
    const unsigned int saveFrequency = std::max< unsigned int >( 5 * static_cast< unsigned int >( onboardComputerRefreshRate ), 1 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE ONBOARD PARAMETERS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Initial conditions
    Eigen::Vector16d initialEstimatedStateVector = Eigen::Vector16d::Zero( ); // assume zero errors in IMU
    initialEstimatedStateVector.segment( 0, 6 ) = translationalInitialState; // assume perfect initial knowledge
    initialEstimatedStateVector.segment( 9, 4 ) = rotationalInitialState.segment( 0, 4 ); // assume perfect initial knowledge
    Eigen::Matrix16d initialEstimatedStateCovarianceMatrix = Eigen::Matrix16d::Identity( );

    // Define instrument accuracy
    const double accelerometerBiasStandardDeviation = 1.0e-4;
    const double accelerometerScaleFactorStandardDeviation = 1.0e-4;
    const double gyroscopeBiasStandardDeviation = 5.0e-9;
    const double gyroscopeScaleFactorStandardDeviation = accelerometerScaleFactorStandardDeviation;
    const Eigen::Vector3d accelerometerAccuracy = Eigen::Vector3d::Constant( 2.0e-4 * std::sqrt( onboardComputerRefreshRate ) );
    const Eigen::Vector3d gyroscopeAccuracy = Eigen::Vector3d::Constant( 3.0e-7 * std::sqrt( onboardComputerRefreshRate ) );
    const Eigen::Vector3d starTrackerAccuracy = Eigen::Vector3d::Constant( 20.0 / 3600.0 );

    const Eigen::Vector3d accelerometerAccuracyAtmosphericPhase =
            Eigen::Vector3d::Constant( 2.0e-4 * std::sqrt( onboardComputerRefreshRateDuringAtmosphericPhase ) );
    const Eigen::Vector3d gyroscopeAccuracyAtmosphericPhase =
            Eigen::Vector3d::Constant( 3.0e-7 * std::sqrt( onboardComputerRefreshRateDuringAtmosphericPhase ) );

    // Define Deep Space Network accuracy
    const double deepSpaceNetworkPositionAccuracy = 10.0;
    const double deepSpaceNetworkVelocityAccuracy = 0.01;
    const double deepSpaceNetworkLightTimeAccuracy = 1.0e-3;

    // Define generic ranging system accuracy
    const Eigen::Vector3d positionAccuracy = Eigen::Vector3d::Constant( 5.0e2 );

    // System and measurment uncertainties
    const double positionStandardDeviation = 1.0e2;
    const double translationalVelocityStandardDeviation = 1.0e-1;
    const double attitudeStandardDeviation = 1.0e-4;
    Eigen::Vector16d diagonalOfSystemUncertainty;
    diagonalOfSystemUncertainty << Eigen::Vector3d::Constant( std::pow( positionStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( translationalVelocityStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( accelerometerBiasStandardDeviation, 2 ) ),
            Eigen::Vector4d::Constant( std::pow( attitudeStandardDeviation, 2 ) ),
            Eigen::Vector3d::Constant( std::pow( gyroscopeBiasStandardDeviation, 2 ) );
    Eigen::Matrix16d systemUncertainty = diagonalOfSystemUncertainty.asDiagonal( );

    Eigen::Vector7d diagonalOfMeasurementUncertainty;
    diagonalOfMeasurementUncertainty.segment( 0, 3 ) = 10.0 * ( positionAccuracy.cwiseProduct( positionAccuracy ) );
    diagonalOfMeasurementUncertainty.segment( 3, 4 ) = Eigen::Vector4d::Constant( 1.0e-6 ); // starTrackerAccuracy <<<<----
    Eigen::Matrix7d measurementUncertainty = diagonalOfMeasurementUncertainty.asDiagonal( );

    // Aerodynamic coefficients
    Eigen::Vector3d onboardAerodynamicCoefficients = Eigen::Vector3d::Zero( );
    onboardAerodynamicCoefficients[ 0 ] = 1.9;

    // Controller gains
    Eigen::Vector3d proportionalGain = Eigen::Vector3d::Constant( 1.0e-5 );
    Eigen::Vector3d integralGain = Eigen::Vector3d::Constant( 0.0 );
    Eigen::Vector3d derivativeGain = Eigen::Vector3d::Constant( 0.05 );

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
    onboardBodySettings[ "Mars" ]->gravityFieldSettings = boost::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );

    AvailableConstantTemperatureAtmosphereModels selectedOnboardAtmosphereModel = three_term_atmosphere_model;
    std::vector< double > vectorOfModelSpecificParameters = { 115.0e3, 2.424e-08, 6533.0, -1.0, 0.0, 0.0 };
    onboardBodySettings[ "Mars" ]->atmosphereSettings = boost::make_shared< CustomConstantTemperatureAtmosphereSettings >(
                selectedOnboardAtmosphereModel, 215.0, 197.0, 1.3, vectorOfModelSpecificParameters );

    NamedBodyMap onboardBodyMap = createBodies( onboardBodySettings );

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
    onboardAccelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 4, 4 ) );
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
    boost::shared_ptr< ControlSystem > controlSystem = boost::make_shared< ControlSystem >( proportionalGain, integralGain, derivativeGain );

    // Create guidance system object
    boost::shared_ptr< GuidanceSystem > guidanceSystem = boost::make_shared< GuidanceSystem >( 255.0e3, 320.0e3, 2800.0, 500.0e3, 0.19, 2.0 );

    // Create unscented Kalman filter settings object for navigation
    boost::shared_ptr< FilterSettings< > > filteringSettings;
    if ( useUnscentedKalmanFilter )
    {
        boost::shared_ptr< IntegratorSettings< > > filterIntegratorSettings =
                boost::make_shared< IntegratorSettings< > >( euler, simulationStartEpoch, onboardComputerRefreshStepSize );
        filteringSettings = boost::make_shared< UnscentedKalmanFilterSettings< > >(
                    systemUncertainty, measurementUncertainty, onboardComputerRefreshStepSize, simulationStartEpoch,
                    initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix, filterIntegratorSettings );
    }
    else
    {
        measurementUncertainty *= 1.0e7;
        boost::shared_ptr< IntegratorSettings< > > filterIntegratorSettings =
                boost::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, onboardComputerRefreshStepSize );
        filteringSettings = boost::make_shared< ExtendedKalmanFilterSettings< > >(
                    systemUncertainty, measurementUncertainty, onboardComputerRefreshStepSize, simulationStartEpoch,
                    initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix, filterIntegratorSettings );
    }

    // Create navigation system object
    boost::shared_ptr< NavigationSystem > navigationSystem = boost::make_shared< NavigationSystem >(
                onboardBodyMap, onboardAccelerationModelMap, "Satellite", "Mars",
                filteringSettings, selectedOnboardAtmosphereModel, marsAtmosphericInterfaceAltitude, marsReducedAtmosphericInterfaceAltitude,
                periapseEstimatorConstant, numberOfRequiredAtmosphereSamplesForInitiation, frequencyOfDeepSpaceNetworkTracking );

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
                std::vector< AerodynamicCoefficientsIndependentVariables >{ angle_of_attack_dependent, altitude_dependent } );

    // Constant radiation pressure variables
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    boost::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Set aerodynamic coefficient and radiation pressure settings
    bodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );
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
    Eigen::Vector3d accelerationBias, accelerometerScaleFactor, gyroscopeBias, gyroscopeScaleFactor, positionBias, positionScaleFactor;
    Eigen::Vector6d accelerometerMisalignment, gyroscopeMisalignment, positionMisalignment;

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
    positionBias.setZero( );
    positionScaleFactor.setZero( );
    positionMisalignment.setZero( );

    // Create navigation instruments object
    boost::shared_ptr< InstrumentsModel > onboardInstruments = boost::make_shared< InstrumentsModel >(
                bodyMap, accelerationModelMap, "Satellite", "Mars" );

    // Add inertial measurement unit
    onboardInstruments->addInertialMeasurementUnit( accelerationBias, accelerometerScaleFactor,
                                                    accelerometerMisalignment, accelerometerAccuracy,
                                                    gyroscopeBias, gyroscopeScaleFactor,
                                                    gyroscopeMisalignment, gyroscopeAccuracy );

    // Add star tracker
    onboardInstruments->addStarTracker( 2, starTrackerAccuracy );

    // Add Deep Space Network tracking
    onboardInstruments->addDeepSpaceNetwork( deepSpaceNetworkPositionAccuracy, deepSpaceNetworkVelocityAccuracy,
                                             deepSpaceNetworkLightTimeAccuracy );

    // Add generic ranging system
    onboardInstruments->addGenericRangingSystem( positionBias, positionScaleFactor, positionMisalignment, positionAccuracy );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE DYNAMICS SETTINGS               ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create onboard computer object
    boost::shared_ptr< OnboardComputerModel > onboardComputer = boost::make_shared< OnboardComputerModel >(
                controlSystem, guidanceSystem, navigationSystem, onboardInstruments, saveFrequency );

    // Define termination conditions
    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    terminationSettingsList.push_back( boost::make_shared< PropagationCustomTerminationSettings >(
                                           boost::bind( &OnboardComputerModel::checkStopCondition, onboardComputer, _1 ) ) );

    boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Satellite", "Mars" );
    terminationSettingsList.push_back( boost::make_shared< PropagationDependentVariableTerminationSettings >(
                                           terminationDependentVariable, 0.95 * marsPropagationAtmosphericInterfaceAltitude, true ) );

    terminationSettingsList.push_back( boost::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );
    terminationSettingsList.push_back( boost::make_shared< PropagationCPUTimeTerminationSettings >( 10800.0 ) );
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );

    // Create integrator settings
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

    // Create propagator settings for rotation
    boost::shared_ptr< RotationalStatePropagatorSettings< > > rotationalPropagatorSettings =
            boost::make_shared< RotationalStatePropagatorSettings< > >(
                torqueModelMap, bodiesToPropagate, rotationalInitialState, terminationSettings, quaternions );

    // Set full propagation settings
    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< > > > fullPropagatorSettings;
    fullPropagatorSettings.push_back( translationalPropagatorSettings );
    fullPropagatorSettings.push_back( rotationalPropagatorSettings );
    boost::shared_ptr< PropagatorSettings< > > propagatorSettings = boost::make_shared< MultiTypePropagatorSettings< > >(
                fullPropagatorSettings, terminationSettings,
                boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create empty dynamics simulator object
    boost::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator =
            boost::make_shared< SingleArcDynamicsSimulator< > >( bodyMap, integratorSettings, propagatorSettings, false );

    // Update onboard models
    onboardInstruments->updateInstruments( simulationStartEpoch );

    // Pre-allocate variables
    std::map< double, Eigen::VectorXd > fullIntegrationResult;
    std::map< double, Eigen::VectorXd > fullDependentVariablesResults;
    std::map< double, Eigen::VectorXd > cartesianTranslationalIntegrationResult;
    std::map< double, Eigen::VectorXd > keplerianTranslationalIntegrationResult;
    std::map< double, Eigen::VectorXd > rotationalIntegrationResult;
    std::map< double, Eigen::VectorXd > dependentVariablesResults;

    // Add initial time and state
    fullIntegrationResult[ simulationStartEpoch ] =
            ( Eigen::VectorXd( 13 ) << translationalInitialState, rotationalInitialState ).finished( );

    // Create loop for each orbit and check condition for stopping aerobraking at the end
    std::vector< bool > vectorOfTerminationConditions;
    do
    {
        // Propagate until a termination condition is reached
        dynamicsSimulator->integrateEquationsOfMotion( propagatorSettings->getInitialStates( ) );

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
            if ( vectorOfTerminationConditions.at( 0 ) ) // onboard computer
            {
                if ( !guidanceSystem->getIsAerobrakingPhaseActive( GuidanceSystem::termination_phase ) )
                {
                    // Add estimated apoapsis maneuver to state
                    finalPropagatedState.segment( 3, 3 ) += controlSystem->getScheduledApsisManeuver( );
                    currentFullIntegrationResult.rbegin( )->second.segment( 3, 3 ) += controlSystem->getScheduledApsisManeuver( );
                }
                else
                {
                    // Add estimated periapsis maneuver to state
                    finalPropagatedState.segment( 3, 3 ) += controlSystem->getScheduledApsisManeuver( );
                    currentFullIntegrationResult.rbegin( )->second.segment( 3, 3 ) += controlSystem->getScheduledApsisManeuver( );
                }
            }
            else if ( vectorOfTerminationConditions.at( 1 ) ) // altitude
            {
                // Extract termination object
                boost::shared_ptr< PropagationDependentVariableTerminationSettings > altitudeTerminationSettings =
                        boost::dynamic_pointer_cast< PropagationDependentVariableTerminationSettings >( terminationSettingsList.at( 1 ) );

                // Reset simulation variables
                if ( altitudeTerminationSettings->useAsLowerLimit_ )
                {
                    // Change step size
                    integratorSettings->initialTimeStep_ = simulationConstantStepSizeDuringAtmosphericPhase;
                    navigationSystem->resetNavigationRefreshStepSize( onboardComputerRefreshStepSizeDuringAtmosphericPhase );
                    onboardInstruments->resetInertialMeasurementUnitRandomNoiseDistribution( accelerometerAccuracyAtmosphericPhase,
                                                                                             gyroscopeAccuracyAtmosphericPhase );

                    // Change altitude termination settings
                    altitudeTerminationSettings->limitValue_ = 1.05 * marsPropagationAtmosphericInterfaceAltitude;
                }
                else
                {
                    // Change step size
                    integratorSettings->initialTimeStep_ = simulationConstantStepSize;
                    navigationSystem->resetNavigationRefreshStepSize( onboardComputerRefreshStepSize );
                    onboardInstruments->resetInertialMeasurementUnitRandomNoiseDistribution( accelerometerAccuracy,
                                                                                             gyroscopeAccuracy );

                    // Change altitude termination settings
                    altitudeTerminationSettings->limitValue_ = 0.95 * marsPropagationAtmosphericInterfaceAltitude;
                }

                // Invert flag for next propagation
                altitudeTerminationSettings->useAsLowerLimit_ = !altitudeTerminationSettings->useAsLowerLimit_;

                // Overwrite propagation termination settings
                terminationSettingsList.at( 1 ) = altitudeTerminationSettings;
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

        // Overwrite dynamics simulator
        dynamicsSimulator = boost::make_shared< SingleArcDynamicsSimulator< > >( bodyMap, integratorSettings, propagatorSettings, false );

        // Save results
        if ( !fullIntegrationResult.empty( ) )
        {
            fullIntegrationResult.erase( integratorSettings->initialTime_ );
            fullDependentVariablesResults.erase( integratorSettings->initialTime_ );
            currentFullIntegrationResult.erase( currentFullIntegrationResult.begin( )->first );
            currentDependentVariablesResults.erase( currentDependentVariablesResults.begin( )->first );
        }
        if ( !fullIntegrationResult.empty( ) )
        {
            fullIntegrationResult.insert( currentFullIntegrationResult.begin( ), currentFullIntegrationResult.end( ) );
            fullDependentVariablesResults.insert( currentDependentVariablesResults.begin( ), currentDependentVariablesResults.end( ) );
        }
    }
    while ( !( onboardComputer->isAerobrakingComplete( ) || vectorOfTerminationConditions.at( 2 ) || vectorOfTerminationConditions.at( 3 ) ) );

    // Propagate for one final orbit, to check orbit configuration
    if ( ( ( simulationEndEpoch - simulationStartEpoch ) / physical_constants::JULIAN_DAY ) > 50.0 )
    {
        terminationSettings = boost::make_shared< PropagationTimeTerminationSettings >(
                    integratorSettings->initialTime_ + physical_constants::JULIAN_DAY );
        boost::dynamic_pointer_cast< SingleArcPropagatorSettings< > >( propagatorSettings )->resetTerminationSettings( terminationSettings );
        dynamicsSimulator = boost::make_shared< SingleArcDynamicsSimulator< > >( bodyMap, integratorSettings, propagatorSettings, false );
        dynamicsSimulator->integrateEquationsOfMotion( propagatorSettings->getInitialStates( ) );

        // Retrieve elements from dynamics simulator
        std::map< double, Eigen::VectorXd > currentFullIntegrationResult = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > currentDependentVariablesResults = dynamicsSimulator->getDependentVariableHistory( );

        // Save results
        if ( !fullIntegrationResult.empty( ) )
        {
            fullIntegrationResult.erase( integratorSettings->initialTime_ );
            fullDependentVariablesResults.erase( integratorSettings->initialTime_ );
            currentFullIntegrationResult.erase( currentFullIntegrationResult.begin( )->first );
            currentDependentVariablesResults.erase( currentDependentVariablesResults.begin( )->first );
        }
        if ( !fullIntegrationResult.empty( ) )
        {
            fullIntegrationResult.insert( currentFullIntegrationResult.begin( ), currentFullIntegrationResult.end( ) );
            fullDependentVariablesResults.insert( currentDependentVariablesResults.begin( ), currentDependentVariablesResults.end( ) );
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             RETRIEVE AND SAVE RESULTS           ///////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Output path
    std::string outputPath = getOutputPath( );

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
        rotationalIntegrationResult[ stateIterator->first ] = currentFullState.segment( 6, 7 );

        // Compute current Keplerian state
        keplerianTranslationalIntegrationResult[ stateIterator->first ] = convertCartesianToKeplerianElements(
                    currentCartesianState, marsGravitationalParameter );

        // Store dependent variables
        dependentVariablesResults[ stateIterator->first ] = fullDependentVariablesResults[ stateIterator->first ];
    }

    // Write propagation history to files
    writeDataMapToTextFile( cartesianTranslationalIntegrationResult, "cartesianPropagated.dat", outputPath );
    writeDataMapToTextFile( keplerianTranslationalIntegrationResult, "keplerianPropagated.dat", outputPath );
    writeDataMapToTextFile( rotationalIntegrationResult, "rotationalPropagated.dat", getOutputPath( ) );
    writeDataMapToTextFile( dependentVariablesResults, "dependentVariables.dat", outputPath );

    // Get estimated states
    std::map< double, std::pair< std::pair< Eigen::Vector6d, Eigen::Vector6d >, Eigen::Vector7d > > fullEstimationResult =
            navigationSystem->getHistoryOfEstimatedStates( );
    std::map< double, Eigen::Vector6d > cartesianTranslationalEstimationResult;
    std::map< double, Eigen::Vector6d > keplerianTranslationalEstimationResult;
    std::map< double, Eigen::Vector7d > rotationalEstimationResult;

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
        rotationalEstimationResult[ stateIterator->first ] = stateIterator->second.second;
        if ( i < vectorOfControlTorques.size( ) )
        {
            controlTorques[ stateIterator->first ] = vectorOfControlTorques.at( i );
        }
    }

    // Check that very last estimate is present
    if ( ( simulationEndEpoch - simulationStartEpoch ) < 50.0 )
    {
        double finalTime = cartesianTranslationalIntegrationResult.rbegin( )->first;
        if ( cartesianTranslationalEstimationResult.count( finalTime ) == 0 )
        {
            std::pair< Eigen::Vector6d, Eigen::Vector6d > finalEstimatedTranslationalState =
                    navigationSystem->getCurrentEstimatedTranslationalState( );
            cartesianTranslationalEstimationResult[ finalTime ] = finalEstimatedTranslationalState.first;
            keplerianTranslationalEstimationResult[ finalTime ] = finalEstimatedTranslationalState.second;
        }
    }

    // Write estimated satellite state hisotry to files
    writeDataMapToTextFile( cartesianTranslationalEstimationResult, "cartesianEstimated.dat", outputPath );
    writeDataMapToTextFile( keplerianTranslationalEstimationResult, "keplerianEstimated.dat", outputPath );
    writeDataMapToTextFile( rotationalEstimationResult, "rotationalEstimated.dat", outputPath );

    // Extract commanded quaternion states
    std::map< unsigned int, Eigen::Vector4d > historyOfCommandedQuaternionStates = controlSystem->getHistoryOfCommandedQuaternions( );
    writeDataMapToTextFile( historyOfCommandedQuaternionStates, "commandedQuaternions.dat", outputPath );

    // Filter results
    if ( extractFilterResults )
    {
        // Get estimated states from filter
        std::map< double, Eigen::VectorXd > fullFilterStateEstimates =
                navigationSystem->getHistoryOfEstimatedStatesFromNavigationFilter( );
        std::map< double, Eigen::MatrixXd > fullFilterCovarianceEstimates =
                navigationSystem->getHistoryOfEstimatedCovarianceFromNavigationFilter( );

        // Extract states and covariances from filter
        unsigned int i = 0;
        std::map< double, Eigen::VectorXd > filterStateEstimates;
        std::map< double, Eigen::VectorXd > filterCovarianceEstimates;
        for ( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = fullFilterStateEstimates.begin( );
              stateIterator != fullFilterStateEstimates.end( ); stateIterator++, i++ )
        {
            if ( ( i % saveFrequency ) == 0 )
            {
                filterStateEstimates[ stateIterator->first ] = stateIterator->second;
                filterCovarianceEstimates[ stateIterator->first ] =
                        Eigen::VectorXd( fullFilterCovarianceEstimates[ stateIterator->first ].diagonal( ) );
            }
        }

        // Write estimated states directly from navigation filter
        writeDataMapToTextFile( filterStateEstimates, "filterStateEstimates.dat", outputPath );
        writeDataMapToTextFile( filterCovarianceEstimates, "filterCovarianceEstimates.dat", outputPath );
    }

    // Get estimated accelerometer errors
    Eigen::Vector3d estimatedAccelerometerErrors = navigationSystem->getEstimatedAccelerometerErrors( );
    std::cout << "Acc. Error: " << estimatedAccelerometerErrors.transpose( ) << std::endl;
    std::cout << "Gyr. Error: " << navigationSystem->getEstimatedGyroscopeErrors( ).transpose( ) << std::endl;

    // Measurements
    if ( extractMeasurementResults )
    {
        // Get instrument measurements and correct them
        unsigned int i = 0;
        std::map< double, Eigen::Vector6d > fullInertialMeasurementUnitMeasurements =
                onboardInstruments->getCurrentOrbitHistoryOfInertialMeasurmentUnitMeasurements( );
        std::map< double, Eigen::Vector3d > fullOnboardExpectedMeasurements =
                navigationSystem->getCurrentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations( );
        std::map< double, Eigen::Vector3d > accelerometerMeasurements;
        std::map< double, Eigen::Vector3d > onboardExpectedMeasurements;
        for ( std::map< double, Eigen::Vector6d >::const_iterator
              measurementIterator = fullInertialMeasurementUnitMeasurements.begin( );
              measurementIterator != fullInertialMeasurementUnitMeasurements.end( ); measurementIterator++, i++ )
        {
            if ( ( i % saveFrequency ) == 0 )
            {
                accelerometerMeasurements[ measurementIterator->first ] =
                        removeErrorsFromInertialMeasurementUnitMeasurement( measurementIterator->second.segment( 0, 3 ),
                                                                            estimatedAccelerometerErrors );
                onboardExpectedMeasurements[ measurementIterator->first ] =
                        fullOnboardExpectedMeasurements[ measurementIterator->first ];
            }
        }

        // Write other data to file
        writeDataMapToTextFile( accelerometerMeasurements, "accelerometerMeasurements.dat", outputPath );
        writeDataMapToTextFile( onboardExpectedMeasurements, "expectedlMeasurements.dat", outputPath );
    }

    // Extract atmosphere data
    if ( extractAtmosphericData )
    {
        std::map< unsigned int, Eigen::VectorXd > atmosphericParameters = navigationSystem->getHistoryOfEstimatedAtmosphereParameters( );
        writeDataMapToTextFile( atmosphericParameters, "atmosphericParameters.dat", outputPath );
    }

    // Extract maneuver information
    if ( extractManeuverInformation )
    {
        std::map< unsigned int, std::pair< double, double > > pairOfPeriapsisCorridorBoundaries =
                guidanceSystem->getHistoryOfEstimatedPeriapsisCorridorBoundaries( );
        std::map< unsigned int, Eigen::Vector2d > periapsisCorridorBoundaries;
        for ( std::map< unsigned int, std::pair< double, double > >::const_iterator mapIterator = pairOfPeriapsisCorridorBoundaries.begin( );
              mapIterator != pairOfPeriapsisCorridorBoundaries.end( ); mapIterator++ )
        {
            periapsisCorridorBoundaries[ mapIterator->first ] =
                    ( Eigen::VectorXd( 2 ) << mapIterator->second.first, mapIterator->second.second ).finished( );
        }
        std::map< unsigned int, double > apoaspsisManeuverMagnitudes = guidanceSystem->getHistoryOfApsisManeuverMagnitudes( ).second;
        writeDataMapToTextFile( periapsisCorridorBoundaries, "periapsisCorridors.dat", outputPath );
        writeDataMapToTextFile( apoaspsisManeuverMagnitudes, "apoapsisManeuver.dat", outputPath );
    }

    // Final statement
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed
    return EXIT_SUCCESS;
}
