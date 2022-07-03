import numpy as np
from scipy import optimize
from Execution import *
from tudatpy.kernel.astro import time_conversion
import datetime
import math

def delta_v(X):
    epoch = X[0]
    tof = X[1]
    total_wet_mass = 1245  # kg
    mass_fuel = 570  # kg
    Isp = 317  # s
    g0 = 9.81  # m/s
    au = 1.495978 * 10 ** (11)  # m
    max_trajectory_limit = 1.25  # AU
    min_trajectory_limit = 0.56  # AU
    penalties = 0

    compute_deltav_func = compute_total_delta_v(epoch, tof)
    delta_v_venus_capture = compute_deltav_func[2]

    constraint_delta_v = Isp * g0 * np.log(total_wet_mass / mass_fuel)

    if delta_v_venus_capture > constraint_delta_v:
        mass_penalty = 1e30
        penalties += mass_penalty

    max_position = compute_deltav_func[5] / au  # in AU
    min_position = compute_deltav_func[6] / au  # in AU

    if max_position > max_trajectory_limit:
        maximum_limit_penalty = 10000 * abs(max_position - max_trajectory_limit)
        penalties += maximum_limit_penalty
        print('MAX DISTANCE LIMIT WAS BREACHED BY: ', abs(max_position - max_trajectory_limit), 'AU')
    if min_position < min_trajectory_limit:
        minimum_limit_penalty = 10000 * abs(min_position - min_trajectory_limit)
        penalties += minimum_limit_penalty
        print('MIN DISTANCE LIMIT WAS BREACHED BY: ', abs(min_position - min_trajectory_limit), ' AU')

    fitness = delta_v_venus_capture + penalties

    return fitness

def banana_function(X):
    x = X[0]
    y = X[1]
    a = 1
    b = 100
    fitness = (a - x) ** 2 + b * (y - x ** 2) ** 2
    return fitness

def schaffer_func(X):
    x = X[0]
    y = X[1]
    fitness = 1 + ((x ** 2 + y ** 2) ** (0.25)) * ((math.sin(50 * (x ** 2 + y ** 2) ** 0.1)) ** 2 + 1)
    return fitness

calendar_date_start = datetime.datetime(2024, 1, 1, 00, 00, 00)
julian_date_start = time_conversion.calendar_date_to_julian_day_since_epoch(calendar_date_start)

calendar_date_end = datetime.datetime(2025, 1, 1, 00, 00, 00)
julian_date_end = time_conversion.calendar_date_to_julian_day_since_epoch(calendar_date_end)

bounds = [[julian_date_start, julian_date_end], [120, 500]]

bounds_test_func = [[-100, 100], [-100, 100]]

result = optimize.differential_evolution(delta_v, bounds, popsize = 50, tol = 1e-2, disp=True)

print(result)