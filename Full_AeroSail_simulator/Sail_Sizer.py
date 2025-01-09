import numpy as np
import matplotlib.pyplot as plt
from Sail import Sail_Class
import Profile as P
import Container_Load_Computations as ContLoad

P.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')
sail_instance = Sail_Class('Data/E473coordinates.txt', 5, 0.4, height=None, panels=20)
interpolation = 'Data/interpolationCR4sail_XFLR5.npz'
sail_instance.load_interpolation(interpolation)

Stackheight = 4
SF=2
max_windspeedknots = 30
cruise_aws = 20
max_windspeed= max_windspeedknots/1.944
full_container_weight = 24390.4

container_load_ratio = 0.8

real_container_weight = full_container_weight*container_load_ratio


panamax_thrust = 2599653.764

initial_height = 20
height_step = 0.01

aspectratio = 60/5

height = initial_height

chord = height*2/aspectratio

failure = False
max_cf = sail_instance.get_cf(plot=True) # Cf with full flaps and alpha
# max_cf = 1.55 # Override Cf
sail_instance.plot_cf_level_curve(max_cf) # Plots allowed area
maxload = 0
while not failure:
    chord = height * 2 / aspectratio
    sail_instance.set_p('height', height)
    sail_instance.set_p('chord', chord)
    for direction in range(360):
        print("Testing height: ", height, "     Testing direction: ", direction)
        forcemag = max_cf * 0.5 * 1.225 * (max_windspeed**2) * chord * height
        if forcemag > maxload:
            maxload = forcemag
        force = np.array(([ forcemag*np.sin(direction*np.pi/180), forcemag*np.cos(direction*np.pi/180) ]))
        ok = ContLoad.CheckContainer(force, height/2, Stackheight, SF=SF, Containerweight=real_container_weight)
        if not ok:
            failure = True
            break
    height += height_step

thrust = 0.5 * 1.225 * (max_windspeed**2) * chord * height * sail_instance.get_max_ct()
thrust_at_cruise = 0.5 * 1.225 * ((cruise_aws / 1.944) ** 2) * chord * height * sail_instance.get_max_ct()
print()
print("Max aparent windspeed: ", max_windspeedknots, " kt")
print("Stack height: ", Stackheight, " Containers")
print("Container load ratio: ", container_load_ratio*100, " %")
print()
print("Max thrust: ", thrust, " N; ", ((thrust/panamax_thrust)*100), "% of a panamax ship")
print("Max thrust at cruise of ", cruise_aws, " kt:", thrust_at_cruise, " N; ", ((thrust_at_cruise / panamax_thrust) * 100), "% of a panamax ship")
print()
print("Height: ", height, " m")
print("Chord: ", chord, " m")
print("Surface area: ", sail_instance.get_p('area'), " m2")
print()
print("Failure: ", failure)
print("Max load: ", maxload, " N")
print()