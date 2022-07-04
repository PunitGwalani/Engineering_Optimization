import numpy as np
import matplotlib.pyplot as plt
import time
from Helper_Function import *
from tudatpy.kernel.astro import time_conversion

spice_interface.load_standard_kernels( )
output_directory = "./Results101/"

def compute_total_delta_v(epoch, tof):
    # start_time = time.time()

    departure_epoch = epoch * constants.JULIAN_DAY
    time_of_flight = tof * constants.JULIAN_DAY
    arrival_epoch = departure_epoch + time_of_flight

    target_body = 'Venus'

    bodies = create_simulation_bodies()
    lambert_ephemeris = get_lambert_problem_result(bodies, target_body, departure_epoch, arrival_epoch)
    dynamics_simulator = propagate_trajectory(departure_epoch, arrival_epoch, bodies, lambert_ephemeris)

    state_history = dynamics_simulator.state_history

    lambert_history = get_lambert_arc_history(lambert_ephemeris, state_history)

    write_propagation_results_to_file(
        dynamics_simulator, lambert_ephemeris, "Trial",output_directory)

    # Delta V calculations
    earth_escape_deltav = departure_delta_v(lambert_history)[0]
    v_infinity_earth = departure_delta_v(lambert_history)[1]
    venus_arrival_deltav = arrival_delta_v(lambert_history)[0]
    v_infinity = arrival_delta_v(lambert_history)[1]
    total_delta_v = earth_escape_deltav + venus_arrival_deltav

    lambert_states = np.vstack(list(lambert_history.values()))
    position_magnitudes = np.sqrt(lambert_states[:, 0]**2 + lambert_states[:, 1]**2 + lambert_states[:, 2]**2)
    maximum_position_magnitude = np.amax(position_magnitudes)
    minimum_position_magnitude = np.amin(position_magnitudes)

    # end_time = time.time()

    # print('Earth Escape Delta V = ', earth_escape_deltav, ' m/s')
    # print('Venus Capture Delta V = ', venus_arrival_deltav, ' m/s')
    # print('Total Delta V = ', total_delta_v, ' m/s')
    # print('Total CPU Time = ', (end_time - start_time), ' s')

    print('Epoch: ', epoch)
    print('ToF: ', tof)

    time = lambert_history.keys()
    time_days = [(t / constants.JULIAN_DAY) - (departure_epoch / constants.JULIAN_DAY) for t in time]

    lambert_cartesian_state = np.vstack(list(lambert_history.values()))

    # fig = plt.figure(figsize=(20, 17))
    # ax = plt.axes(projection='3d')
    #
    # dependent_variables = dynamics_simulator.dependent_variable_history
    # dependent_variable_history = np.vstack(list(dependent_variables.values()))
    #
    # ax.plot(lambert_cartesian_state[:, 0], lambert_cartesian_state[:, 1], lambert_cartesian_state[:, 2],
    #         color='g', label='Lambert Trajectory', linewidth=1.5)
    #
    # ax.plot(dependent_variable_history[:, 0], dependent_variable_history[:, 1], dependent_variable_history[:, 2],
    #         color='b', label='Earth\'s Orbit')
    # ax.plot(dependent_variable_history[:, 3], dependent_variable_history[:, 4], dependent_variable_history[:, 5],
    #         color='tab:orange', label='Venus\'s Orbit')
    #
    # ax.scatter3D(lambert_cartesian_state[0, 0], lambert_cartesian_state[0, 1], lambert_cartesian_state[0, 2],
    #              s=100, color='c', label='Earth')
    # ax.scatter3D(lambert_cartesian_state[-1, 0], lambert_cartesian_state[-1, 1], lambert_cartesian_state[-1, 2],
    #              s=90, color='r', label='Venus')
    # ax.scatter3D(0, 0, 0, s=90, color='y', label='Sun')
    # ax.set_xlabel('X [m]')
    # ax.set_ylabel('Y [m]')
    # ax.set_zlabel('Z [m]')
    # ax.legend()
    # ax.set_title("Optimal Lambert arc from Earth to Venus", fontweight='bold', fontsize=15)
    # plt.show()

    print(time_conversion.julian_day_to_calendar_date(epoch + 2451545.0))

    return [total_delta_v, earth_escape_deltav, venus_arrival_deltav, v_infinity, v_infinity_earth,
            maximum_position_magnitude, minimum_position_magnitude]

# Trajectory Plot for visualization
# compute_total_delta_v(9108.466675980939, 159.62500120650586)

