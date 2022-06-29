import datetime
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from Execution import *
from tudatpy.kernel.astro import time_conversion

# epoch = 3813.075185
# tof = 185.4590589

calendar_date_start = datetime.datetime(2024, 1, 1, 00, 00, 00)
julian_date_start = time_conversion.calendar_date_to_julian_day_since_epoch(calendar_date_start)
print(julian_date_start)

calendar_date_end = datetime.datetime(2025, 1, 1, 00, 00, 00)
julian_date_end = time_conversion.calendar_date_to_julian_day_since_epoch(calendar_date_end)
print(julian_date_end)

epoch_size = 150
tof_size = 300
rows = epoch_size*tof_size

epoch_array = np.linspace(julian_date_start, julian_date_end, epoch_size)
tof_array = np.linspace(120, 500, tof_size)

# data_list = []
#
#
# for epoch in epoch_array:
#     for tof in tof_array:
#         final_delta_v = compute_total_delta_v(epoch, tof)
#         data_list.append([epoch, tof, final_delta_v])
#
# data_array = np.array(data_list)

X, Y = np.meshgrid(epoch_array, tof_array)
final_delta_v = np.zeros((tof_size, epoch_size))
earth_escape = np.zeros((tof_size, epoch_size))
venus_arrival = np.zeros((tof_size, epoch_size))
excess_velocity = np.zeros((tof_size, epoch_size))
excess_velocity_earth = np.zeros((tof_size, epoch_size))

for i in range(epoch_size):
    for j in range(tof_size):
        start_time = time.time()
        output = compute_total_delta_v(X[0, i], Y[j, 0])
        final_delta_v[j, i] = (output[0])/1000
        earth_escape[j, i] = (output[1])/1000
        venus_arrival[j, i] = (output[2])/1000
        excess_velocity[j, i] = (output[3])/1000
        excess_velocity_earth[j, i] = (output[4]) / 1000
        end_time = time.time()
        print('Time:', end_time-start_time, 's')

# Saving the data generated
np.savetxt('./Results101/All_files/Final_DepartureEpochs.dat', X)
np.savetxt('./Results101/All_files/Final_TimeOfFlight.dat', Y)
np.savetxt('./Results101/All_files/Final_DeltaV.dat', final_delta_v)
np.savetxt('./Results101/All_files/Earth_Escape_DeltaV.dat', earth_escape)
np.savetxt('./Results101/All_files/Venus_Capture_DeltaV.dat', venus_arrival)
np.savetxt('./Results101/All_files/Excess_Velocity_Venus.dat', excess_velocity)
np.savetxt('./Results101/All_files/Excess_Velocity_Earth.dat', excess_velocity_earth)

# Plotting the Porkchop plot
fig1, ax = plt.subplots()
CS = ax.contour(X, Y, final_delta_v, cmap='jet')
ax.clabel(CS, inline=True, fontsize=8)
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel(r'$\Delta V$ [km/s]', fontweight = 'bold')
ax.set_xlabel('Departure Epoch [Julian Days from J2000]', fontweight = 'bold')
ax.set_ylabel('Time of Flight [Days]', fontweight = 'bold')
plt.show()