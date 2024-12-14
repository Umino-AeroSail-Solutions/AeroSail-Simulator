import numpy as np
import matplotlib.pyplot as plt
from Sail import Sail_Class
import Profile as P
import Container_Load_Computations as ContLoad

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

viridis = cm.get_cmap('viridis', 12)

# Initialize Xfoil
P.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')

# Create a Sail instance and load the interpolation
sail_instance = Sail_Class('Data/E473coordinates.txt', 5, 0.4, height=30, panels=20)  # Example height of 20
interpolation = 'Data/interpolationCR4sail_XFLR5.npz'
sail_instance.load_interpolation(interpolation)

# Set up the parameters
height = 30  # Example height
chord = 5
windspeed = 30 / 1.944  # Max windspeed in m/s

Stackheight = 4
SF = 1.5
full_container_weight = 24390.4
container_load_ratio = 0.8
real_container_weight = full_container_weight * container_load_ratio

def drawenvelope(sail_instance, height, chord, windspeedknots, ax, real_container_weight=real_container_weight, SF=1.5, Stackheight=Stackheight):
    failure = False
    windspeed = windspeedknots / 1.944
    maxcf = 0
    maxcf_step = 0.001
    while not failure:
        sail_instance.set_p('height', height)
        sail_instance.set_p('chord', chord)
        for direction in range(360):
            print("Testing height: ", height, "     Testing direction: ", direction)
            forcemag = maxcf * 0.5 * 1.225 * (windspeed ** 2) * chord * height
            force = np.array([forcemag * np.sin(direction * np.pi / 180), forcemag * np.cos(direction * np.pi / 180)])
            ok = ContLoad.CheckContainer(force, height / 2, Stackheight, SF=SF, Containerweight=real_container_weight)
            if not ok:
                failure = True
                break
        maxcf += maxcf_step

    contour = sail_instance.plot_cf_level_curve(maxcf, ax=ax, title='Envelopes', filling=True)
    for c in contour.collections:
        c.set_label(f"{windspeedknots} knots")
    # Add text labels on the curve for i, c in enumerate(contour.collections):
    for i, c in enumerate(contour.collections):
        for path in c.get_paths():
            vertices = path.vertices
            midpoint = vertices[len(vertices) // 2]
            ax.text(midpoint[0], midpoint[1], f"{windspeedknots} knots", color=color, fontsize=8, ha='center', va='center', backgroundcolor='white')
        return contour


# testing_windspeeds = np.arange(20, 50, 10)
testing_windspeeds = [20, 30, 50]
# Plot the allowed envelopes for different wind speeds on the same figure
fig, ax = plt.subplots()
colors = viridis(np.linspace(0, 1, len(testing_windspeeds)))

# Dictionary to hold contour plots for the legend
contour_handles = {}

for windspeed, color in zip(testing_windspeeds, colors):
    contour = drawenvelope(sail_instance, height, chord, windspeed, ax)
    for c in contour.collections:
        c.set_label(f"{windspeed} knots")

# Add a legend to the plot
ax.legend()
plt.show()
