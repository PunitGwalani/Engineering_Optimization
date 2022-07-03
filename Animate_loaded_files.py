import matplotlib.animation as animation
import matplotlib.pyplot as plt
from Plotting_LoadedFiles import *
import numpy as np

from tudatpy.kernel.astro import time_conversion
import datetime

x_positions = np.loadtxt('./Results101/All_files/x_position_final.dat')
y_positions = np.loadtxt('./Results101/All_files/y_position_final.dat')
x_velocity = np.loadtxt('./Results101/All_files/x_velocity_final.dat')
y_velocity = np.loadtxt('./Results101/All_files/y_velocity_final.dat')
fitness = np.loadtxt('./Results101/All_files/Fitness_best_final.dat')

x_pos_updated = np.add(x_positions, x_velocity)
y_pos_updated = np.add(y_positions, y_velocity)

fitness = np.loadtxt('./Results101/All_files/Fitness_best_final.dat')

calendar_date_start = datetime.datetime(2024, 1, 1, 00, 00, 00)
julian_date_start = time_conversion.calendar_date_to_julian_day_since_epoch(calendar_date_start)

calendar_date_end = datetime.datetime(2025, 1, 1, 00, 00, 00)
julian_date_end = time_conversion.calendar_date_to_julian_day_since_epoch(calendar_date_end)

lower_bound=[julian_date_start, 120]
upper_bound=[julian_date_end, 500]

fig, ax = plt.subplots()

def animate(i):
    ax.cla()
    ax.scatter(x_positions[i, :], y_positions[i, :], marker='x', facecolor='red', zorder=5)
    plot_contour(fig, ax)
    ax.set_xlim(lower_bound[0], upper_bound[0])
    ax.set_ylim(lower_bound[1], upper_bound[1])
    ax.set_title(r'Particle Swarm Optimization [$\Delta$ V (km/s)] - Iterations: ' + str(i), fontweight = 'bold')
    # ax.quiver(x_positions[i,:], y_positions[i, :], x_pos_updated[i, :],
    #           y_pos_updated[i, :], color = 'b', scale = 1)

# plt.grid(b=None)
ani = animation.FuncAnimation(fig, animate, interval=250, frames=range(len(fitness)), repeat=False, blit=False)

ani.save('PSO_animation.gif', writer='pillow', dpi=720)

# Plotting the fitness plot

fig2, ax2 = plt.subplots()
ax2.plot(range(len(fitness)), fitness/1000, linewidth = 1.5)
ax2.set_xlim(0, 100)
# ax2.set_ylim(0.02, np.amax(fitness/1000))
# ax2.set_yscale('log')
ax2.set_xlabel('Iterations [-]', fontweight = 'bold')
ax2.set_ylabel('Best Global Fitness [km/s]', fontweight = 'bold')
ax2.set_title('Fitness Curve - Particle Swarm Optimization', fontweight = 'bold')
fig2.savefig('Fitness_Curve.png')
plt.show()