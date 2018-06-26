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
#include "Tudat/Mathematics/Filters/unscentedKalmanFilter.h"

#include "Tudat/Astrodynamics/GuidanceNavigationControl/control.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/guidance.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/instrument.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigation.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/onboardComputer.h"

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

//! Function to compute the state derivative for both the translational and rotational motion of the spacecraft.
typedef Eigen::Matrix< double, 13, 1 > Eigen::Vector13d;
Eigen::Vector13d stateDerivativeFromOnboardModels( const double currentTime, const Eigen::VectorXd& currentStateVector,
                                                   const Eigen::VectorXd& currentControlVector, const double marsGravitationalParameter,
                                                   const double marsRadius, const double secondOrderSphericalHarmonicCoefficient,
                                                   const double currentAirDensity, const double ballisticCoefficient,
                                                   const double dragCoefficient, const double liftCoefficient,
                                                   const Eigen::Matrix3d& inertiaMatrix )
{
    // Declare state derivative vector
    Eigen::Vector13d stateDerivative;

    // Add translational motion
    stateDerivative[ 0 ] = currentStateVector[ 3 ];
    stateDerivative[ 1 ] = currentStateVector[ 4 ];
    stateDerivative[ 2 ] = currentStateVector[ 5 ];

    double currentPositionInZDirectionSquared = 5.0 * currentStateVector[ 2 ] * currentStateVector[ 2 ];
    double currentRadialDistance = currentStateVector.segment( 0, 3 ).norm( );
    double currentRadialDistanceSquared = currentRadialDistance * currentRadialDistance;
    double centralGravitationalTerm = - marsGravitationalParameter / currentRadialDistanceSquared;
    double gravityRecurringTerm = 3.0 / 2.0 * secondOrderSphericalHarmonicCoefficient * std::pow( marsRadius /
                                                                                           currentRadialDistanceSquared, 2 );
    stateDerivative[ 3 ] = centralGravitationalTerm * ( 1.0 + gravityRecurringTerm * ( currentRadialDistanceSquared -
                                                                                       currentPositionInZDirectionSquared ) );
    stateDerivative[ 4 ] = stateDerivative[ 3 ];
    stateDerivative[ 5 ] = centralGravitationalTerm * ( 1.0 + gravityRecurringTerm * ( 3.0 * currentRadialDistanceSquared -
                                                                                       currentPositionInZDirectionSquared ) );
    stateDerivative.segment( 3, 3 ) *= currentStateVector.segment( 0, 3 ).normalized( );

    double aerodynamicRecurringTerm = - 0.5 * currentAirDensity / ballisticCoefficient *
            currentStateVector.segment( 3, 3 ).squaredNorm( );

    // Add rotational motion

    // Give back result
    return stateDerivative;
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
    using namespace tudat::unit_conversions;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings
    const double simulationStartEpoch = 7.0 * tudat::physical_constants::JULIAN_YEAR +
            30.0 * 6.0 * tudat::physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = 50.0 * tudat::physical_constants::JULIAN_DAY + simulationStartEpoch;

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
    std::vector< double > extrapolationValues = { 0.0, 0.0, 186.813, 8183.0, 1.667 };

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
    const double vehicleMass = 1000.0;
    bodyMap[ "Satellite" ]->setConstantBodyMass( vehicleMass );

    Eigen::Matrix3d inertiaTensor = Eigen::Matrix3d::Zero( );
    inertiaTensor( 0, 0 ) = 5750.0;
    inertiaTensor( 1, 1 ) = 1200.0;
    inertiaTensor( 2, 2 ) = 5200.0;
    bodyMap[ "Satellite" ]->setBodyInertiaTensor( inertiaTensor );

    // Aerodynamic coefficients from file
    std::map< int, std::string > aerodynamicForceCoefficientFiles;
    std::map< int, std::string > aerodynamicMomentCoefficientFiles;
    aerodynamicForceCoefficientFiles[ 0 ] = "/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/"
                                            "University/Master Thesis/Code/MATLAB/data/MRODragCoefficients.txt";
    aerodynamicForceCoefficientFiles[ 2 ] = "/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/"
                                            "University/Master Thesis/Code/MATLAB/data/MROLiftCoefficients.txt";
    aerodynamicMomentCoefficientFiles[ 1 ] = "/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/"
                                             "University/Master Thesis/Code/MATLAB/data/MROMomentCoefficients.txt";

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

    // Set initial Keplerian elements for satellite
    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 27228500;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.869218;
    initialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 158.7 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 23.4 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 180.0 );

    double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    double marsRadius = bodyMap.at( "Mars" )->getShapeModel( )->getAverageRadius( );
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
    double onboardComputerRefreshRate = 0.1;
    Eigen::Vector3d accelerationBias, accelerometerScaleFactor, accelerometerAccuracy,
            gyroscopeBias, gyroscopeScaleFactor, gyroscopeAccuracy, starTrackerAccuracy;
    Eigen::Vector6d accelerometerMisalignment, gyroscopeMisalignment;

    std::default_random_engine generator;
    std::normal_distribution< double > distributionOne( 0.0, 5.0e-9 );
    std::normal_distribution< double > distributionTwo( 0.0, 3.0e-7 / std::sqrt( onboardComputerRefreshRate ) );
    std::normal_distribution< double > distributionThree( 0.0, 1.0e-4 );
    std::normal_distribution< double > distributionFour( 0.0, 2.0e-4 / std::sqrt( onboardComputerRefreshRate ) );
    std::normal_distribution< double > distributionFive( 0.0, 1.0e-3 );
    for ( unsigned int i = 0; i < 6; i++ )
    {
        if ( i < 3 )
        {
            accelerationBias[ i ] = distributionThree( generator );
            accelerometerScaleFactor[ i ] = distributionThree( generator );
            accelerometerAccuracy[ i ] = distributionFour( generator );
            gyroscopeBias[ i ] = distributionOne( generator );
            gyroscopeScaleFactor[ i ] = distributionThree( generator );
            gyroscopeAccuracy[ i ] = distributionTwo( generator );
            starTrackerAccuracy[ i ] = distributionThree( generator );
        }
        accelerometerMisalignment[ i ] = distributionFive( generator );
        gyroscopeMisalignment[ i ] = distributionFive( generator );
    }

    // Create intertial measurement unit object
    boost::shared_ptr< InstrumentSystem > onboardInstruments = boost::make_shared< InstrumentSystem >(
                bodyMap, accelerationMap, centralBodyMap, "Satellite" );
    onboardInstruments->addInertialMeasurementUnit( accelerationBias, accelerometerScaleFactor,
                                                    accelerometerMisalignment, accelerometerAccuracy,
                                                    gyroscopeBias, gyroscopeScaleFactor,
                                                    gyroscopeMisalignment, gyroscopeAccuracy );
    onboardInstruments->addStarTracker( 2, starTrackerAccuracy );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE GNC MODELS                      ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create control system object
    boost::shared_ptr< ControlSystem > controlSystem = boost::make_shared< ControlSystem >( );

    // Create unscented Kalman filter object for navigation
    Eigen::VectorXd initialEstimatedStateVector;
    initialEstimatedStateVector << translationalInitialState, rotationalInitialState;
    Eigen::MatrixXd initialEstimatedStateCovarianceMatrix = Eigen::MatrixXd::Identity(
                initialEstimatedStateVector.rows( ), initialEstimatedStateVector.rows( ) );

    Eigen::MatrixXd systemUncertainty;
    double positionStandardDeviation = 1.0;
    double translationalVelocityStandardDeviation = 1.0e-3;
    double attitudeStandardDeviation = 1.0e-6;
    double rotationalVelocityStandardDeviation = 1.0e-9;
    Eigen::VectorXd diagonalOfSystemUncertainty;
    diagonalOfSystemUncertainty << Eigen::VectorXd::Constant( 3, std::pow( positionStandardDeviation, 2 ) ),
            Eigen::VectorXd::Constant( 3, std::pow( translationalVelocityStandardDeviation, 2 ) ),
            Eigen::VectorXd::Constant( 4, std::pow( attitudeStandardDeviation, 2 ) ),
            Eigen::VectorXd::Constant( 3, std::pow( rotationalVelocityStandardDeviation, 2 ) );
    systemUncertainty = diagonalOfSystemUncertainty.asDiagonal( );

    Eigen::MatrixXd measurementUncertainty;

    boost::shared_ptr< IntegratorSettings< > > kalmanFilterIntegratorSettings =
            boost::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, onboardComputerRefreshRate );

    boost::shared_ptr< UnscentedKalmanFilter< > > unscentedKalmanFilter = boost::make_shared< UnscentedKalmanFilter< > >(
                boost::bind( &stateDerivativeFromOnboardModels, _1, _2, controlSystem->getCurrentControlVector( ),
                             marsGravitationalParameter, marsRadius ),
                boost::lambda::constant( &onboardInstruments->getCurrentInertialMeasurementUnitMeasurement ),
                systemUncertainty, measurementUncertainty,
                simulationStartEpoch, initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix,
                kalmanFilterIntegratorSettings );

    // Create giudance system object
    boost::shared_ptr< GuidanceSystem > guidanceSystem = boost::make_shared< GuidanceSystem >( );

    // Create navigation system object
    boost::shared_ptr< NavigationSystem > navigationSystem = boost::make_shared< NavigationSystem >(
                unscentedKalmanFilter, &stateTransitionMatrix, marsGravitationalParameter, marsRadius );

    // Create onboard computer object
    boost::shared_ptr< OnboardComputer > onboardComputer = boost::make_shared< OnboardComputer >(
                controlSystem, guidanceSystem, navigationSystem, onboardInstruments );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE DYNAMICS SETTINGS               ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define termination conditions
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< CustomTerminationSettings >( boost::bind( &onboardComputer->checkStoppingCondition, _1 ) );

    // Create integrator settings
    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                rungeKuttaVariableStepSize, simulationStartEpoch, 100.0,
                RungeKuttaCoefficients::rungeKuttaFehlberg56, 1e-5, 1e5, 1e-13, 1e-13 );

    // Create propagation settings for translational dynamics
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, translationalInitialState,
                terminationSettings, unified_state_model_exponential_map );

    // Create propagator settings for rotation
    boost::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            boost::make_shared< RotationalStatePropagatorSettings< double > >(
                torqueModelMap, bodiesToPropagate, rotationalInitialState, terminationSettings, exponential_map );

    // Set full propagation settings
    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    propagatorSettingsList.push_back( rotationalPropagatorSettings );
    boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            boost::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsList, terminationSettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false, true );
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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Write perturbed satellite propagation history to file
    writeDataMapToTextFile( cartesianTranslationalIntegrationResult, "CartesianTranslational.dat", getOutputPath( ) );
    writeDataMapToTextFile( keplerianTranslationalIntegrationResult, "KeplerianTranslational.dat", getOutputPath( ) );
    writeDataMapToTextFile( quaternionsRotationalIntegrationResult, "rotational.dat", getOutputPath( ) );

    // Final statement
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed
    return EXIT_SUCCESS;
}
