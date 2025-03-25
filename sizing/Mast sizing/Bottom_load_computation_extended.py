import numpy as np
import matplotlib.pyplot as plt
from Multi_segment_mast_sizer import Segment

total_force_vector = np.array(([3200.99, 32009.88]))

segment_4_added_weight = 114
w_4, h_4, d_4 = 1.24, 0.5,0.02
OL_34 = 1.03
L_4 = 10.0

segment_3_added_weight = 114
w_3, h_3, d_3 = 1.28, 0.54, 0.04
OL_23 = 1.4
L_3 = 9.5

segment_2_added_weight = 114
w_2, h_2, d_2 = 1.36, 0.62, 0.08
OL_12 = 1.97
L_2 = 9.0

segment_1_added_weight = 114
w_1, h_1, d_1 = 1.44, 0.7, 0.08
OL_01 = 2.6
L_1 = 8.5

segment_4 = Segment(segment_4_added_weight, L_4,OL_34, 0, 0,0, total_force_vector, 30)
segment_4.compute_internal_loads(plot=True)
segment_3_top_force_top, segment_3_top_force_bottom = segment_4.compute_reaction_loads()

segment_3 = Segment(segment_3_added_weight, L_3,OL_23, OL_34, segment_3_top_force_top,segment_3_top_force_bottom , total_force_vector, 30)
segment_3.compute_internal_loads(plot=True)
segment_2_top_force_top, segment_2_top_force_bottom = segment_3.compute_reaction_loads()

segment_2 = Segment(segment_2_added_weight, L_2,OL_12, OL_23, segment_2_top_force_top,segment_2_top_force_bottom , total_force_vector, 30)
segment_2.compute_internal_loads(plot=True)
segment_1_top_force_top, segment_1_top_force_bottom = segment_2.compute_reaction_loads()

segment_1 = Segment(segment_1_added_weight, L_1,OL_01, OL_12, segment_1_top_force_top,segment_1_top_force_bottom , total_force_vector, 30)
segment_1.compute_internal_loads(plot=True)
force_top, force_bottom = segment_2.compute_reaction_loads()

print("Top ring support max operational load is: ", np.linalg.norm(force_top), " N")
print("Bottom ring support max operational load is: ", np.linalg.norm(force_bottom), " N")