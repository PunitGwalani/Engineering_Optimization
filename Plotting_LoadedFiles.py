import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from tudatpy.kernel.astro import time_conversion

# X = np.load('./Results101/Final_DepartureEpochs_og.dat', allow_pickle=True)
# Y = np.load('./Results101/Final_TimeOfFlight_og.dat', allow_pickle=True)
# Z = np.load('./Results101/Final_DeltaV_og.dat', allow_pickle=False)

X = np.loadtxt('./Results101/Final_DepartureEpochs_og.dat')
Y = np.loadtxt('./Results101/Final_TimeOfFlight_og.dat')
Z = np.loadtxt('./Results101/Final_DeltaV_og.dat')

# X = np.genfromtxt('./Results101/Final_DepartureEpochs_og.dat',
#                      skip_header=1,
#                      skip_footer=1,
#                      names=True,
#                      dtype=None,
#                      delimiter=' ')
#
# Y = np.genfromtxt('./Results101/Final_TimeOfFlight_og.dat',
#                      skip_header=1,
#                      skip_footer=1,
#                      names=True,
#                      dtype=None,
#                      delimiter=' ')
#
# Z = np.genfromtxt('./Results101/Final_DeltaV_og.dat',
#                      skip_header=1,
#                      skip_footer=1,
#                      names=True,
#                      dtype=None,
#                      delimiter=' ')

Z1 = np.transpose(Z)

# X_juliandays = X + 2451545.0
# print(time_conversion.julian_day_to_calendar_date(X_juliandays[0,0]))
# X_calendardates = np.empty(np.shape(X_juliandays), dtype = object)
#
# for i in range(len(X_calendardates[:,0])):
#     for j in range(len(X_calendardates[0, :])):
#         print('[', i,j, ']')
#         X_calendardates[i, j] = time_conversion.julian_day_to_calendar_date(X_juliandays[i,j]).strftime("%d %b, %Y")
#         print(X_calendardates[i,j])

fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z1, levels = 20, cmap='jet')
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