import numpy as np
import matplotlib.pyplot as plt
import XFLR5_Scraper as XFScrp
import os
from pathlib import Path

# Filenames have to end in -flapdeflection as (_XX) followed by .txt: randomassfilename-04.txt
# Alphas shall be equal arrays for all files
dir = 'Data/XFLR5_5_30_0,5_10m_s_INTERPOLATION/'


def crt_XFLR5_interpolation(dir):
    files = [name for name in os.listdir(dir) if os.path.isfile(os.path.join(dir, name))]
    InterpFlaps = np.zeros(len(files), dtype=float)
    results_list = []
    InterpAlphas = None

    for i, file in enumerate(files):
        flapdeflection = np.radians(float((Path(file).stem).split('-')[-1]))
        alphas, CL, CDi, CDv, CD, CY, Cl, Cm, Cn, Cni, QInf, XCP = XFScrp.read_polar_data(dir + '/' + file)

        if InterpAlphas is None:
            InterpAlphas = alphas

        results_list.append((CL, CDi, CDv, CD, CY, Cl, Cm, Cn, Cni, QInf, XCP))
        InterpFlaps[i] = flapdeflection

    Alphas, Flaps = np.meshgrid(InterpAlphas, InterpFlaps)
    Cl = np.zeros((len(files), len(InterpAlphas)), dtype=float)
    Cd = np.zeros((len(files), len(InterpAlphas)), dtype=float)
    CloCd = np.zeros((len(files), len(InterpAlphas)), dtype=float)

    for i, result in enumerate(results_list):
        Cl[i, :] = result[0]
        Cd[i, :] = result[3]
        CloCd[i, :] = Cl[i,:]

    CloCd = np.divide(Cl, Cd, out=np.zeros_like(Cl), where=Cd != 0)
    return InterpAlphas, InterpFlaps, Cl, Cd, CloCd


# Example usage:
InterpAlphas, InterpFlaps, Cl, Cd, CloCd = crt_XFLR5_interpolation(dir)

# Print some example data
# print('InterpAlphas:', InterpAlphas)
# print('InterpFlaps:', InterpFlaps)
