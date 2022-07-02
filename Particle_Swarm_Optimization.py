import copy
import sys
import datetime

import numpy as np
import random
import math
import matplotlib.animation as animation

from Execution import *
from tudatpy.kernel.astro import time_conversion

# Fitness function

test_function = False

if test_function:
    # rastrigin function
    def fitness_rastrigin(position1, position2):
        fitnessVal = 0.0
        for i in range(2):
            if i == 0:
                xi = position1
            else:
                xi = position2
            fitnessVal += (xi * xi) - (10 * math.cos(2 * math.pi * xi)) + 10
        return fitnessVal
else:
    def delta_v(epoch, tof):

        total_wet_mass = 1245 #kg
        mass_fuel = 570 #kg
        Isp = 317 #s
        g0 = 9.81 #m/s
        au = 1.495978 * 10**(11) #m
        max_trajectory_limit = 1.25 #AU
        min_trajectory_limit = 0.56 #AU
        penalties = 0

        compute_deltav_func = compute_total_delta_v(epoch, tof)
        delta_v_venus_capture = compute_deltav_func[2]

        constraint_delta_v = Isp*g0*np.log(total_wet_mass/mass_fuel)

        if delta_v_venus_capture > constraint_delta_v:
            mass_penalty = 1e30
            penalties += mass_penalty

        max_position = compute_deltav_func[5]/au    # in AU
        min_position = compute_deltav_func[6]/au    # in AU

        if max_position > max_trajectory_limit:
            maximum_limit_penalty = 10000*abs(max_position - max_trajectory_limit)
            penalties += maximum_limit_penalty
            print('MAX DISTANCE LIMIT WAS BREACHED BY: ', abs(max_position - max_trajectory_limit), 'AU')
        if min_position < min_trajectory_limit:
            minimum_limit_penalty = 10000*abs(min_position - min_trajectory_limit)
            penalties += minimum_limit_penalty
            print('MIN DISTANCE LIMIT WAS BREACHED BY: ', abs(min_position - min_trajectory_limit), ' AU')


        fitness = delta_v_venus_capture + penalties

        return fitness

# Defining the class for PSO

class Particles:
    def __init__(self, fitness_func, lower_bound, upper_bound):
        self.position_indv = []  # indv --> individual
        self.velocity_indv = []
        self.best_position_indv = []
        self.best_indv_fitness = sys.float_info.max

        for i in range(len(lower_bound)):
            random_pos_value = (upper_bound[i] - lower_bound[i])*random.random() + lower_bound[i]
            self.position_indv.append(random_pos_value)
            random_vel_value = 0 # (upper_bound[i] - lower_bound[i])*random.random() + lower_bound[i]
            self.velocity_indv.append(random_vel_value)

        self.best_position_indv = copy.deepcopy(self.position_indv)
        self.fitness = fitness_func(self.position_indv[0], self.position_indv[1])
        self.best_indv_fitness = copy.deepcopy(self.fitness)

# Defining the fuction for PSO, w = 0.729, c1 = 1.49445, c2 = 1.49445

def pso(fitness_func, iterations, particles, lower_bound, upper_bound, tol, p, w = 0.729, c1 = 1.49445, c2 = 1.49445):

    best_position_global = []
    best_fitness_global = sys.float_info.max
    x_positions_for_plot = np.zeros((iterations, particles))
    y_positions_for_plot = np.zeros((iterations, particles))
    x_velocity_for_plot = np.zeros((iterations, particles))
    y_velocity_for_plot = np.zeros((iterations, particles))

    best_fitness_for_plot = np.zeros((iterations,))

    swarm = [Particles(fitness_func, lower_bound, upper_bound) for i in range(particles)]

    for i in range(particles):
        if swarm[i].fitness < best_fitness_global:
            best_position_global = copy.deepcopy(swarm[i].position_indv)
            best_fitness_global = swarm[i].fitness
        x_positions_for_plot[0, i] = swarm[i].position_indv[0]
        y_positions_for_plot[0, i] = swarm[i].position_indv[1]

        x_velocity_for_plot[0, i] = swarm[i].velocity_indv[0]
        y_velocity_for_plot[0, i] = swarm[i].velocity_indv[1]

    best_fitness_for_plot[0] = best_fitness_global

    k = 0

    while k < iterations:
        # Print statement to be ADDED

        fitness_each_iteration = []

        print('Iteration = ', k)
        print('Best Position = ', best_position_global)
        print('Best Fitness = ', best_fitness_global)
        print('---------------')

        for i in range(particles):
           for j in range(len(lower_bound)):
               r1 = random.random()
               r2 = random.random()

               inertial_term = (w*swarm[i].velocity_indv[j])
               congnitive_term = (c1*r1*(swarm[i].best_position_indv[j] - swarm[i].position_indv[j]))
               social_term = (c2*r2*(best_position_global[j] - swarm[i].position_indv[j]))

               swarm[i].velocity_indv[j] = inertial_term + congnitive_term + social_term


               # if swarm[i].velocity_indv[j] < lower_bound[j]:
               #     swarm[i].velocity_indv[j] = lower_bound[j]
               if swarm[i].velocity_indv[j] > upper_bound[j]:
                   swarm[i].velocity_indv[j] = upper_bound[j]

           for j in range(len(lower_bound)):
               swarm[i].position_indv[j] += swarm[i].velocity_indv[j]

               if swarm[i].position_indv[j] < lower_bound[j]:
                   swarm[i].position_indv[j] = lower_bound[j]
               elif swarm[i].position_indv[j] > upper_bound[j]:
                   swarm[i].position_indv[j] = upper_bound[j]


           swarm[i].fitness = fitness_func(swarm[i].position_indv[0], swarm[i].position_indv[1])

           if swarm[i].fitness < swarm[i].best_indv_fitness:
               swarm[i].best_position_indv = swarm[i].fitness
               swarm[i].best_position_indv = copy.deepcopy(swarm[i].position_indv)

           if swarm[i].fitness < best_fitness_global:
               best_fitness_global = swarm[i].fitness
               best_position_global = copy.deepcopy(swarm[i].position_indv)

           x_positions_for_plot[k, i] = swarm[i].position_indv[0]
           y_positions_for_plot[k, i] = swarm[i].position_indv[1]
           x_velocity_for_plot[k, i] = swarm[i].velocity_indv[0]
           y_velocity_for_plot[k, i] = swarm[i].velocity_indv[1]

           fitness_each_iteration.append(swarm[i].fitness)

        best_fitness_for_plot[k] = best_fitness_global
        # print(fitness_each_iteration)

        # implementing a stopping criteria
        print('Minimum Fitness each iteration:', min(fitness_each_iteration))
        print('Maximum Fitness each iteration:', max(fitness_each_iteration))
        print('Difference between fitnesses: ', abs(min(fitness_each_iteration) - max(fitness_each_iteration)))
        print('***********')
        if abs(min(fitness_each_iteration) - max(fitness_each_iteration)) < tol:
            counter = 0
            for flag in range(int(p*particles)):
                value = random.randint(0, particles-1)
                if abs(min(fitness_each_iteration) - fitness_each_iteration[value]) < tol:
                    counter = counter + 1
            if int(counter) == int(p*particles):
                break

        k = k + 1

    return [best_position_global, best_fitness_global, x_positions_for_plot, y_positions_for_plot,
            best_fitness_for_plot, x_velocity_for_plot, y_velocity_for_plot, k]

# Test Function

calendar_date_start = datetime.datetime(2024, 1, 1, 00, 00, 00)
julian_date_start = time_conversion.calendar_date_to_julian_day_since_epoch(calendar_date_start)

calendar_date_end = datetime.datetime(2025, 1, 1, 00, 00, 00)
julian_date_end = time_conversion.calendar_date_to_julian_day_since_epoch(calendar_date_end)

if test_function:
    lowerbound = [-10, -10]
    upperbound = [10, 10]

    optimals = pso(fitness_rastrigin, iterations=100, particles=50, lower_bound=lowerbound,
                   upper_bound=upperbound, tol=1e-5, p=0.9)
else:

    lowerbound = [julian_date_start, 120]
    upperbound = [julian_date_end, 500]

    optimals = pso(delta_v,iterations=100, particles=70, lower_bound=lowerbound,
                   upper_bound=upperbound, tol=1e-1, p=0.9)



print('Optimal Position: ', optimals[0])
print('Optimal fitness: ', optimals[1])
print('Iterations taken: ', optimals[7])

np.savetxt('./Results101/All_files/x_position.dat', optimals[2])
np.savetxt('./Results101/All_files/y_position.dat', optimals[3])
np.savetxt('./Results101/All_files/Fitness_best.dat', optimals[4])
np.savetxt('./Results101/All_files/x_velocity.dat', optimals[5])
np.savetxt('./Results101/All_files/y_velocity.dat', optimals[6])


