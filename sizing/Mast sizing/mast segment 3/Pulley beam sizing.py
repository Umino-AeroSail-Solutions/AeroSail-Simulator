import numpy as np
import beam_sizer_copy as bs

alu_6063_T66 = bs.Material("alu_6063_T66",Density=2710, YieldStrength=200000000.0, max_shear=152000000.0)

Pulley4_beam = bs.Square_beam(0.274, 0.06, 0.06, material=alu_6063_T66)

aeeoplatform_weight = 114
Mast_segment_4_weight = 195 + aeeoplatform_weight

total_pulley_load = 2*Mast_segment_4_weight*9.81/2
print("Total pulley 4 load: ", total_pulley_load, "N")

shear = [total_pulley_load,0]
moment = [0, total_pulley_load*Pulley4_beam.length/2]

Pulley4_beam.size_profile(shear[0], shear[1], 0, moment[0], moment[1], 1.5, 1.5/1000,plot_convergence=True)

Mast_segment_3_weight = 540 + aeeoplatform_weight
total_pulley_load = 2*Mast_segment_3_weight*9.81/2 + 4*Mast_segment_4_weight*9.81/2
print("Total pulley 3 load: ", total_pulley_load, "N")

Pulley3_beam = bs.Square_beam(0.438, 0.06, 0.06, material=alu_6063_T66)

shear = [total_pulley_load/2,0]
moment = [0, total_pulley_load*Pulley3_beam.length/4]

Pulley3_beam.size_profile(shear[0], shear[1], 0, moment[0], moment[1], 1.5, 3/1000,plot_convergence=True)

Mast_segment_2_weight = 788.4 + aeeoplatform_weight
total_pulley_load = 2*Mast_segment_2_weight*9.81/2 + 4*Mast_segment_3_weight*9.81/2 - 2*Mast_segment_3_weight*9.81/2
print("Total pulley 2 load: ", total_pulley_load, "N")

Pulley2_beam = bs.Square_beam(0.642, 0.06, 0.06, material=alu_6063_T66)

shear = [total_pulley_load/3,0]
moment = [0, total_pulley_load*Pulley3_beam.length/6]

Pulley2_beam.size_profile(shear[0], shear[1], 0, moment[0], moment[1], 1.5, 3/1000,plot_convergence=True)

bottom_support_beam = bs.Square_beam(1.804, 0.04, 0.04, material=alu_6063_T66)


print("Mast support beam \n")

shear = [0.577*1000, 0]
moment= [0, 1.227*1000]
bottom_support_beam.size_profile(shear[0], shear[1], 0, moment[0], moment[1], 1.5, 4/1000,plot_convergence=True)