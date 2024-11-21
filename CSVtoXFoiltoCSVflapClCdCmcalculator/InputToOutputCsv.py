import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from pyxfoil import Xfoil, set_workdir, set_xfoilexe

from MakeE473flapAirfoil import makeflap
from DATtoArray import DATtoArray
from GetExcelAlphaDelta import getExcelAlphaDelta

set_workdir('C:/Xfoil699src')
set_xfoilexe('C:/Xfoil699src/xfoil.exe')

plainAirfoil = "E473coordinates.txt"
Alphas, Deltas = getExcelAlphaDelta("Input.txt")

ClCdCm = []

mach = 0
re = 1000000

chordratio = 0.4


for index in range(np.size(Alphas)):
    print(Alphas[index], Deltas[index])
    Airfoil = makeflap(plainAirfoil, chordratio , Deltas[index])
    x = Airfoil[:, 0].tolist()
    y = Airfoil[:, 1].tolist()
    xfoil = Xfoil('Flapped E473')
    xfoil.set_points(x, y)
    xfoil.set_ppar(160)  # Number of panels
    rescase = xfoil.run_result(Alphas[index], mach=mach, re=re)
    almin = Alphas[index]
    almax = Alphas[index]
    alint = 0
    polar = xfoil.run_polar(almin, almax, alint, mach=mach, re=re)
    ClCdCm.append( [polar.cl[0], polar.cd[0], polar.cm[0]] )
Data = pd.DataFrame(np.array(ClCdCm))
Data.to_csv('OutputCsv.csv', header=False,index=False)
print(ClCdCm)

