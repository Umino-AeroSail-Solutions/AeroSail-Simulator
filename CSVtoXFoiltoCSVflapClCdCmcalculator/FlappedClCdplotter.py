from MakeE473flapAirfoil import makeflap
from pyxfoil import Xfoil, set_workdir, set_xfoilexe
import numpy as np
import matplotlib.pyplot as plt

set_workdir('C:/Xfoil699src')
set_xfoilexe('C:/Xfoil699src/xfoil.exe')
chordfraction = 0.4
mach = 0.1
re = 1000000
testDeltas = [0, 3, 11, 15, 20, 35]
for delta in testDeltas:
    foil = makeflap('E473coordinates.txt', chordfraction, np.radians(delta))
    x = foil[:, 0].tolist()
    y = foil[:, 1].tolist()
    xfoil = Xfoil('Flapped E473')
    # xfoil.points_from_dat('airfoil.dat')
    xfoil.set_points(x, y)
    xfoil.set_ppar(160) #Number of panels
    # Plots xfoil airfoil profile
    # ax1 = xfoil.plot_profile(ls='-')
    rescase = xfoil.run_result(5, mach=mach, re=re)
    # Run the xfoil in polar mode from alpha min to alpha max with interval
    almin = -5
    almax = 25
    alint = 0.5
    polar1 = xfoil.run_polar(almin, almax, alint, mach=mach, re=re)
    ax2 = None
    axp1 = None
    axp2 = None
    axp3 = None
    for resname in xfoil.results:
        # ax2 = xfoil.results[resname].plot_result(yaxis='cp', ax=ax2, ls='-x')
        # _ = ax2.legend()
        # axp1 = polar1.plot_polar(ax=axp1, xaxis='cl', yaxis='cd', ls='-o')
        # _ = axp1.legend()
        # axp2 = polar1.plot_polar(ax=axp2, xaxis='alpha', yaxis='clocd', ls='-o')
        # _ = axp2.legend()
        axp3 = polar1.plot_polar(ax=axp3, xaxis='alpha', yaxis='cl', ls='-o')
        _ = axp3.legend()
        plt.show()
# Shows plots for cases in xfoil cases


plt.show()