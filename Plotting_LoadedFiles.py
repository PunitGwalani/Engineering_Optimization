import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib import cm
from tudatpy.kernel.astro import time_conversion

def plot_contour(fig, ax):
    X = np.loadtxt('./Results101/All_files/Final_DepartureEpochs_final.dat')
    Y = np.loadtxt('./Results101/All_files/Final_TimeOfFlight_final.dat')
    Z2 = np.loadtxt('./Results101/All_files/Venus_Capture_DeltaV_final.dat')
    print(np.amin(Z2))

    # fig, ax = plt.subplots()
    loc = matplotlib.ticker.MaxNLocator(30)
    lvls = loc.tick_values(Z2.min(), Z2.max())

    CS = ax.contour(X, Y, Z2, colors='black', levels = [0.05, 0.1, 0.3, 0.5, 1, 2, 4, 6, 8, 10, 30], zorder = 1)

    for line, lvl in zip(CS.collections, CS.levels):
        if lvl < 1:
            line.set_linestyle(':')
        elif lvl <= 5 and lvl > 1:
            line.set_linestyle('--')
        else:
            line.set_linestyle('-')

    # ax.set_facecolor('#d0d0d0')
    ax.set_xticks(np.linspace(min(X[0,:]), max(X[0,:]), 5))
    print(np.linspace(min(X[0,:]), max(X[0,:]),5))
    X_calendardates = np.empty(np.shape(np.linspace(min(X[0,:]), max(X[0,:]), 5)), dtype = object)
    X_juliandays = np.linspace(min(X[0,:]), max(X[0,:]), 5) + 2451545.0

    for j in range(len(X_calendardates)):
        print('[',j,']')
        X_calendardates[j] = time_conversion.julian_day_to_calendar_date(X_juliandays[j]).strftime("%d %b, %Y")
        print(X_calendardates[j])

    ax.set_xticklabels(X_calendardates)
    ax.clabel(CS, inline=True, fontsize=9)
    # cbar = fig.colorbar(CS)
    # cbar.ax.set_ylabel(r'$\Delta V$ [km/s]', fontweight = 'bold')
    ax.set_xlabel('Departure Dates', fontweight = 'bold')
    ax.set_ylabel('Time of Flight [Days]', fontweight = 'bold')
# plt.show()