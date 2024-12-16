import Cross_section_analysis as cs
import numpy as np

# Finds suitable design with the least posible area

# Vslot options [[area, d, Ixx/Iyy], [...,...]]
# t-slot 2020, t-slot 4040, t-slot 8080 (quad 4040 version)
v_slots = np.array([[.00016210, 0.02, .0000000066881] , [0.000504, 0.04, 0.000000073272], [.0016122, 0.08, .0000010225207]])

w, h = 0.6, 0.8

for v_slot in v_slots:
    a = v_slot[0]
    d = v_slot[1]
    cornerIxx = v_slot[3]
