import random
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import math

def schaffer_func(x, y):
    fitness = 1 + ((x**2 + y**2)**(0.25))*((math.sin(50*(x**2 + y**2)**0.1))**2 + 1)
    return fitness

def banana_function(x, y):
    a = 1
    b = 100
    fitness = (a - x)**2 + b*(y - x**2)**2
    return fitness

def fitness_rastrigin(position1, position2):
    fitnessVal = 0.0
    for i in range(2):
        if i == 0:
            xi = position1
        else:
            xi = position2
        fitnessVal += (xi * xi) - (10 * math.cos(2 * math.pi * xi)) + 10
    return fitnessVal

x = np.linspace(-2, 2, 50)
y = np.linspace(-1, 3, 50)

X, Y = np.meshgrid(x, y)
f2 = np.vectorize(banana_function)
Z = f2(X,Y)

three_D = False

if three_D:
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    # ax.contour3D(X, Y, Z, 100, cmap='jet')
    ax.plot_surface(X, Y, Z, cmap="jet", lw=0.5, rstride=1, cstride=1)
    ax.set_zlabel('z', fontweight='bold')
else:
    fig, ax = plt.subplots()
    ax.contour(X, Y, Z, 30, cmap='jet')
# ax.plot_surface(X, Y, Z, cmap="jet", lw=0.5, rstride=1, cstride=1)
# ax.contour(X, Y, Z, 10, lw=3, cmap="jet", linestyles="solid", offset=-1)
# ax.contour(X, Y, Z, 10, lw=3, colors="k", linestyles="solid")

ax.set_xlabel('x', fontweight ='bold')
ax.set_ylabel('y', fontweight ='bold')
ax.set_title(r'Rosenbrock’s Banana Test Function', fontweight ='bold')
plt.show()

# fig, ax = plt.subplots()
# pop_sizes = [10, 20, 35, 50, 70, 100]
# fitnesses = [0.02682974755878579, 1.059683087379614e-07,  8.648911190989914e-10, 6.187334410657202e-10,
#              1.6286571568886448e-10,
#              4.0897589853516964e-11]
# ax.plot(pop_sizes, fitnesses, '-o')
# ax.set_yscale('log')
# ax.set_xlabel('Population Sizes', fontweight ='bold')
# ax.set_ylabel('Average Optimal Values', fontweight ='bold')
# ax.set_title('PSO on Rosenbrock’s Banana Test Function', fontweight ='bold')
# plt.show()

