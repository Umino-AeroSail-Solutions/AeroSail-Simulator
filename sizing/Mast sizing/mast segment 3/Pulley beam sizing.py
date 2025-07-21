import numpy as np
import beam_sizer_copy as bs

alu_6063_T66 = bs.Material("alu_6063_T66",Density=2710, YieldStrength=200000000.0, max_shear=152000000.0)

Pulley_beam = bs.Square_beam(0.274, 0.06, 0.06, material=alu_6063_T66)

aeeoplatform_weight = 114
Mast_segment_4_weight = 195 + aeeoplatform_weight

total_pulley_load = 2*Mast_segment_4_weight*9.81/2
print("Total pulley load: ", total_pulley_load, "N")

shear = [total_pulley_load,0]
moment = [0, total_pulley_load*Pulley_beam.length/2]

Pulley_beam.size_profile(shear[0], shear[1], 0, moment[0], moment[1], 1.5, 1.5/1000,plot_convergence=True)