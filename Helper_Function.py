import numpy as np
from tudatpy.io import save2txt
from tudatpy.kernel import constants
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import estimation_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel.astro import two_body_dynamics
from tudatpy.kernel.astro import element_conversion

# Define departure/arrival epoch - in seconds since J2000
departure_epoch = 3813.075185 * constants.JULIAN_DAY
time_of_flight = 185.4590589 * constants.JULIAN_DAY
arrival_epoch = departure_epoch + time_of_flight
target_body = "Venus"
global_frame_orientation = 'ECLIPJ2000'
fixed_step_size = 3600.0*3

# Departure conditons
perigee_radius_at_earth = 500*1000
apogee_radius_at_earth = 500*1000

# Arrival Conditions
perigee_radius_at_venus = 250*1000
apogee_radius_at_venus = 66000*1000

# To be changed
def write_propagation_results_to_file(
        dynamics_simulator: numerical_simulation.SingleArcSimulator,
        lambert_arc_ephemeris: environment.Ephemeris,
        file_output_identifier: str,
        output_directory: str):

    # Save numerical states
    simulation_result = dynamics_simulator.state_history
    save2txt(solution=simulation_result,
             filename=output_directory + file_output_identifier + "_numerical_states.dat",
             directory="./")

    # Save dependent variables
    dependent_variables = dynamics_simulator.dependent_variable_history
    if len(dependent_variables.keys()) > 0:
        save2txt(solution=dependent_variables,
                 filename=output_directory + file_output_identifier + "_dependent_variables.dat",
                 directory="./")

    # Save Lambert arc states
    lambert_arc_states = get_lambert_arc_history(lambert_arc_ephemeris, simulation_result)

    save2txt(solution=lambert_arc_states,
             filename=output_directory + file_output_identifier + "_lambert_states.dat",
             directory="./")

    return

def get_lambert_problem_result(
        bodies: environment.SystemOfBodies,
        target_body: str,
        departure_epoch: float,
        arrival_epoch: float ) -> environment.Ephemeris:

    # Gravitational parameter of the Sun
    central_body_gravitational_parameter = bodies.get_body("Sun").gravitational_parameter

    # Set initial and final positions for Lambert targeter
    initial_state = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name="Earth",
        observer_body_name="Sun",
        reference_frame_name=global_frame_orientation,
        aberration_corrections="NONE",
        ephemeris_time=departure_epoch)

    final_state = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name=target_body,
        observer_body_name="Sun",
        reference_frame_name=global_frame_orientation,
        aberration_corrections="NONE",
        ephemeris_time=arrival_epoch)

    # Create Lambert targeter
    lambertTargeter = two_body_dynamics.LambertTargeterIzzo(
        initial_state[:3], final_state[:3], arrival_epoch - departure_epoch, central_body_gravitational_parameter);

    # Compute initial Cartesian state of Lambert arc
    lambert_arc_initial_state = initial_state
    lambert_arc_initial_state[3:] = lambertTargeter.get_departure_velocity()

    # Compute Keplerian state of Lambert arc
    lambert_arc_keplerian_elements = element_conversion.cartesian_to_keplerian(lambert_arc_initial_state,
                                                                       central_body_gravitational_parameter)

    # Setup Keplerian ephemeris model that describes the Lambert arc
    kepler_ephemeris = environment_setup.create_body_ephemeris(
        environment_setup.ephemeris.keplerian(lambert_arc_keplerian_elements, departure_epoch,
                                              central_body_gravitational_parameter), "")

    return kepler_ephemeris

def get_lambert_arc_history(
        lambert_arc_ephemeris: environment.Ephemeris,
        simulation_result: dict ) -> dict:

    lambert_arc_states = dict()
    for state in simulation_result:
        lambert_arc_states[state] = lambert_arc_ephemeris.cartesian_state(state)

    return lambert_arc_states

def create_simulation_bodies( ) -> environment.SystemOfBodies:

    bodies_to_create = ['Sun', 'Earth', 'Venus']
    global_frame_origin = 'Sun'
    global_frame_orientation = 'ECLIPJ2000'
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation)

    bodies = environment_setup.create_system_of_bodies(body_settings)

    bodies.create_empty_body('Spacecraft')

    bodies.get_body('Spacecraft').mass = 1000.0

    return bodies

def propagate_trajectory(
        initial_time: float,
        final_time: float,
        bodies: environment.SystemOfBodies,
        lambert_arc_ephemeris: environment.Ephemeris) -> numerical_simulation.SingleArcSimulator:

    lambert_arc_initial_state = lambert_arc_ephemeris.cartesian_state(initial_time)

    propagator_settings = get_unperturbed_propagator_settings(
            bodies, lambert_arc_initial_state, final_time)


    integrator_settings = propagation_setup.integrator.runge_kutta_4(initial_time, fixed_step_size)
    dynamics_simulator = numerical_simulation.SingleArcSimulator(bodies, integrator_settings, propagator_settings)

    return dynamics_simulator

def get_unperturbed_propagator_settings(
        bodies: environment.SystemOfBodies,
        initial_state: np.array,
        termination_time: float ) -> propagation_setup.propagator.SingleArcPropagatorSettings:

    bodies_to_propagate = ['Spacecraft']
    central_bodies = ['Sun']
    acceleration_settings_on_spacecraft = dict(
        Sun=[
            propagation_setup.acceleration.point_mass_gravity()
        ]
    )
    acceleration_settings = {'Spacecraft': acceleration_settings_on_spacecraft}
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, bodies_to_propagate, central_bodies)

    dependent_variables_to_save = [
        propagation_setup.dependent_variable.relative_position('Earth', 'Sun'),
        propagation_setup.dependent_variable.relative_position('Venus', 'Sun'),
        propagation_setup.dependent_variable.relative_distance('Spacecraft', 'Venus')
    ]

    # Create propagation settings.
    termination_settings = propagation_setup.propagator.time_termination(termination_time)
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        termination_settings,
        output_variables = dependent_variables_to_save
    )

    return propagator_settings

def departure_delta_v(
        lambert_states: dict):

    mu_earth = spice_interface.get_body_gravitational_parameter('Earth')
    velocity_vector = np.vstack(list(lambert_states.values()))[0, 3:6]
    initial_epoch = [t for t in lambert_states.keys()][0]

    radius_at_perigee = spice_interface.get_average_radius('Earth') + perigee_radius_at_earth
    radius_at_apogee = spice_interface.get_average_radius('Earth') + apogee_radius_at_earth
    eccentricity = (radius_at_apogee - radius_at_perigee) / (radius_at_apogee + radius_at_apogee)

    state_earth = spice_interface.get_body_cartesian_state_at_epoch('Earth', 'SSB', 'ECLIPJ2000', 'NONE',
                                                                    initial_epoch)
    velocity_vector_earth = state_earth[3:6]
    excess_velocity = velocity_vector - velocity_vector_earth
    excess_velocity_magnitude = np.sqrt(excess_velocity[0] ** 2 + excess_velocity[1] ** 2 + excess_velocity[2] ** 2)

    delta_v_earth_departure = np.sqrt(((2 * mu_earth) / radius_at_perigee) + excess_velocity_magnitude ** 2) - (
            1 + eccentricity) * np.sqrt(mu_earth / radius_at_perigee)

    return delta_v_earth_departure

def arrival_delta_v(
        lambert_states: dict):

    # Delta V required for orbit insertion around venus
    mu_venus = spice_interface.get_body_gravitational_parameter('Venus')
    velocity_vector = np.vstack(list(lambert_states.values()))[-1, 3:6]
    final_epoch = [t for t in lambert_states.keys()][-1]

    radius_at_perigee = spice_interface.get_average_radius('Venus') + perigee_radius_at_venus
    radius_at_apogee = spice_interface.get_average_radius('Venus') + apogee_radius_at_venus
    eccentricity = (radius_at_apogee - radius_at_perigee) / (radius_at_apogee + radius_at_apogee)
    
    state_venus = spice_interface.get_body_cartesian_state_at_epoch('Venus', 'SSB', 'ECLIPJ2000', 'NONE',
                                                                   final_epoch)
    velocity_vector_venus = state_venus[3:6]
    excess_velocity = velocity_vector - velocity_vector_venus
    excess_velocity_magnitude = np.sqrt(excess_velocity[0] ** 2 + excess_velocity[1] ** 2 + excess_velocity[2] ** 2)

    delta_v_venus_arrival = np.sqrt(((2 * mu_venus) / radius_at_perigee) + excess_velocity_magnitude ** 2) - (
            1 + eccentricity) * np.sqrt(
        mu_venus / radius_at_perigee)

    return abs(delta_v_venus_arrival)