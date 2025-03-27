import numpy as np
import matplotlib.pyplot as plt
from Multi_segment_mast_sizer import Segment

cf = 0.92
windspeed_knots = 100
windspeed = windspeed_knots/1.944
rho = 1.225

exposed_height = 10-2.6

chord = 5

area = exposed_height*chord

force_mag = 0.5*rho*cf*(windspeed**2)
total_force = np.array(([force_mag, 0]))

mast = Segment(114*4, 10, 2.6,0,0,0,total_force,exposed_height)
mast.compute_internal_loads(plot=True)
force_top, force_bottom = mast.compute_reaction_loads()

print("Top ring support max operational load is: ", np.linalg.norm(force_top), " N")
print("Bottom ring support max operational load is: ", np.linalg.norm(force_bottom), " N")