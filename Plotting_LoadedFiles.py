import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from tudatpy.kernel.astro import time_conversion

X = np.loadtxt('./Results101/All_files/Final_DepartureEpochs_final.dat')
Y = np.loadtxt('./Results101/All_files/Final_TimeOfFlight_final.dat')
Z = np.loadtxt('./Results101/All_files/Final_DeltaV_final.dat')

fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z, levels = 20, cmap='jet')
ax.set_xticks(np.linspace(min(X[0,:]), max(X[0,:]), 5))
print(np.linspace(min(X[0,:]), max(X[0,:]),5))
X_calendardates = np.empty(np.shape(np.linspace(min(X[0,:]), max(X[0,:]), 5)), dtype = object)
X_juliandays = np.linspace(min(X[0,:]), max(X[0,:]), 5) + 2451545.0

for j in range(len(X_calendardates)):
    print('[',j,']')
    X_calendardates[j] = time_conversion.julian_day_to_calendar_date(X_juliandays[j]).strftime("%d %b, %Y")
    print(X_calendardates[j])

ax.set_xticklabels(X_calendardates)
ax.clabel(CS, inline=True, fontsize=5)
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel(r'$\Delta V$ [km/s]', fontweight = 'bold')
ax.set_xlabel('Departure Dates', fontweight = 'bold')
ax.set_ylabel('Time of Flight [Days]', fontweight = 'bold')
plt.show()