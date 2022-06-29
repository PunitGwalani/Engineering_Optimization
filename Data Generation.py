import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from Execution import *
from tudatpy.kernel.astro import time_conversion

# epoch = 3813.075185
# tof = 185.4590589

calendar_date_start = datetime.datetime(2023, 1, 1, 00, 00, 00)
julian_date_start = time_conversion.calendar_date_to_julian_day_since_epoch(calendar_date_start)
print(julian_date_start)

calendar_date_end = datetime.datetime(2025, 1, 1, 00, 00, 00)
julian_date_end = time_conversion.calendar_date_to_julian_day_since_epoch(calendar_date_end)
print(julian_date_end)

epoch_size = 200
tof_size = 400
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

for i in range(epoch_size):
    for j in range(tof_size):
        final_delta_v[j, i] = compute_total_delta_v(X[0, i], Y[j, 0])/1000

# Saving the data generated
np.savetxt('./Results101/Final_DepartureEpochs_2years.dat', X)
np.savetxt('./Results101/Final_TimeOfFlight_2years.dat', Y)
np.savetxt('./Results101/Final_DeltaV_2years.dat', final_delta_v)

# Plotting the Porkchop plot
fig1, ax = plt.subplots()
CS = ax.contour(X, Y, final_delta_v, cmap='jet')
ax.clabel(CS, inline=True, fontsize=8)
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel(r'$\Delta V$ [km/s]', fontweight = 'bold')
ax.set_xlabel('Departure Epoch [Julian Days from J2000]', fontweight = 'bold')
ax.set_ylabel('Time of Flight [Days]', fontweight = 'bold')
plt.show()