import numpy as np
import matplotlib.pyplot as plt
from Sail import Sail_Class
import Profile as P
import Container_Load_Computations as ContLoad

P.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')
sail_instance = Sail_Class('Data/E473coordinates.txt', 5, 0.4, height=None, panels=20)
interpolation = 'Data/interpolationCR4sail_XFLR5.npz'
sail_instance.load_interpolation(interpolation)

Stackheight = 2
SF=1.3
max_windspeedknots = 30
max_windspeed= max_windspeedknots/1.944

panamax_thrust =

initial_height = 5
height_step = 0.01

aspectratio = 60/5

height = initial_height

chord = height*2/aspectratio

failure = False

maxload = 0
while not failure:
    chord = height * 2 / aspectratio
    sail_instance.set_p('height', height)
    sail_instance.set_p('chord', chord)
    max_cf = sail_instance.get_cf()
    for direction in range(360):
        print("Testing height: ", height, "     Testing direction: ", direction)
        forcemag = max_cf * 0.5 * 1.225 * (max_windspeed**2) * chord * height
        if forcemag > maxload:
            maxload = forcemag
        force = np.array(([ forcemag*np.sin(direction*np.pi/180), forcemag*np.cos(direction*np.pi/180) ]))
        ok = ContLoad.CheckContainer(force, height/2, Stackheight, SF=SF)
        if not ok:
            failure = True
            break
    height += height_step
print("Max aparent windspeed: ", max_windspeedknots, " kt")
print("Height: ", height, " m")
print("Chord: ", chord, " m")
print("Surface area: ", sail_instance.get_p('area'), " m2")
print("Failure: ", failure)
print("Max load: ", maxload, " N")
print("Stack height: ", Stackheight, " Containers")
print()
print()