import numpy as np
import matplotlib.pyplot as plt
from DATtoArray import DATtoArray

def makeflap(datfile, chordratio, deltaradians):

    # WARNING: THIS ONLY WORKS UP TO A DEFLECTION AROUND 30 DEGREES

    plainfoil = DATtoArray(datfile)
    flapfoil = plainfoil
    npoints = len(plainfoil)
    for i in range(npoints):
        rotatedpointx = (1 - chordratio) + ((np.cos(-deltaradians) * (plainfoil[i, 0] - (1 - chordratio))) - (np.sin(-deltaradians) * (plainfoil[i, 1])))
        rotatedpointy = ((np.sin(-deltaradians) * (plainfoil[i, 0] - (1 - chordratio))) + (np.cos(-deltaradians) * (plainfoil[i, 1])))
        if flapfoil[i, 0] > (1-chordratio):
            if deltaradians > 0:
                if rotatedpointy < plainfoil[i, 1]:
                    flapfoil[i, 0] = rotatedpointx
                    flapfoil[i, 1] = rotatedpointy
                elif deltaradians > np.radians(36):
                    flapfoil[i, 1] = -plainfoil[i, 1]

            elif rotatedpointy > plainfoil[i, 1]:
                flapfoil[i, 0] = rotatedpointx
                flapfoil[i, 1] = rotatedpointy
            elif deltaradians < np.radians(-36):
                flapfoil[i, 1] = -plainfoil[i, 1]

    np.savetxt('airfoil.dat', flapfoil, fmt=['%.3f', '%.3f'])
    return flapfoil

# Generating new foil coordinates
# newfoil = makeflap("E473coordinates.txt", 0.5, np.radians(20))
# np.savetxt('output.dat', newfoil, fmt=['%.3f','%.3f'])
# # Scatter plot
# plt.scatter(newfoil[:, 0], newfoil[:, 1])
# plt.title('New Flap Airfoil Scatter Plot')
# plt.xlabel('x-coordinate')
# plt.ylabel('y-coordinate')
# plt.grid(True)
# plt.axis('equal')  # Ensure the aspect ratio is equal
# plt.show()
