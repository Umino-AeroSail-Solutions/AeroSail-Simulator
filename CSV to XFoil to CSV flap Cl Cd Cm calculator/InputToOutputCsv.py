import numpy as np
import matplotlib.pyplot as plt

from pyxfoil import Xfoil, set_workdir, set_xfoilexe

from MakeE473flapAirfoil import makeflap
from DATtoArray import DATtoArray
from GetExcelAlphaDelta import getExcelAlphaDelta

plainAirfoil = DATtoArray("E473coordinates.txt")
AlphasDeltas = getExcelAlphaDelta("Input.txt")

ClCdCm = np.array((AlphasDeltas.len(), 3))

mach = 0
re = 1000000

chordratio = 0.4

for AlphaDelta in AlphasDeltas:
    Airfoil = makeflap(plainAirfoil, chordratio , AlphaDelta[1])
    x = Airfoil[:, 0].tolist()
    y = Airfoil[:, 1].tolist()
    xfoil = Xfoil('Flapped E473')
    xfoil.set_points(x, y)
    xfoil.set_ppar(160)  # Number of panels
    rescase = xfoil.run_result(AlphaDelta[0], mach=mach, re=re)
    almin = -5
    almax = 25
    alint = 0.5
    polar1 = xfoil.run_polar(almin, almax, alint, mach=mach, re=re)