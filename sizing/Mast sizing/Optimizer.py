import numpy as np
import Cross_section_analysis as cs
# Finds suitable design with the least possible area

# V-slot options [[area, d, Ixx/Iyy], [...,...]]
# t-slot 2020, t-slot 4040, t-slot 8080 (quad 4040 version)
v_slots = np.array([[.00016210, 0.02, .0000000066881], [0.000504, 0.04, 0.000000073272], [.0016122, 0.08, .0000010225207]])

w, h = 0.6, 0.8

SF = 1.5

max_tension = 200000000.0
max_shear = 283000000.0

Vx, Vy = 500000, 500000
Mx, My = -500000, 500000

Vx, Vy, Mx, My = Vx * SF, Vy * SF, Mx * SF, My * SF

skin_step_thickness = 0.0005
min_skin_thickness = 0.0005
max_thickness = 0.01

possible_designs = []

subdivisions = 100

for v_slot in v_slots:
    a = v_slot[0]
    d = v_slot[1]
    cornerIxx = v_slot[2]
    for top_thickness in np.arange(min_skin_thickness, max_thickness, skin_step_thickness):
        for side_thickness in np.arange(min_skin_thickness, max_thickness, skin_step_thickness):
            areas, thicknesses = cs.create_areas_and_thicknesses(w, h, d, a, subdivisions, top_thickness, side_thickness)
            Ixx, Iyy, Ixy = cs.compute_Ixx_Iyy_Ixy(w, h, top_thickness, side_thickness, d, a, cornerIxx=cornerIxx)
            bending_ok, bending_sf = cs.check_bending_ok(w, h, top_thickness, side_thickness, d, a, Mx, My, max_tension, cornerIxx=cornerIxx)
            cs.update_areas_bending(areas, thicknesses, Mx, My, Ixx, Iyy, Ixy)
            qb, qbmax = cs.compute_shear_flows(areas, Vx, Vy, subdivisions, thicknesses, h, d, w)
            top_ok, sides_ok, top_SF, sides_SF = cs.check_shear_ok(qb, top_thickness, side_thickness, max_shear, thicknesses)
            if top_ok and sides_ok and bending_ok:
                bending_sf=bending_sf*SF
                top_SF=top_SF*SF
                sides_SF=sides_SF*SF
                total_area = 2 * (top_thickness * w + side_thickness * h) + 4 * a
                possible_designs.append([d, top_thickness, side_thickness, sides_SF, top_SF, bending_sf,total_area])
                break
possible_designs = np.array(possible_designs)
# print(possible_designs)
if np.size(possible_designs) == 0:
    print("There is no possible designs available :(")
    print("Tip: Increase the thickness limit or add bigger vslot channels")
    print("\nTip: Changing weapons is faster than reloading.......")
else:
    # Sort the possible designs by total area (last column)
    sorted_designs = possible_designs[possible_designs[:, -1].argsort()]

    # Output the design with the lowest total area
    lowest_area_design = sorted_designs[0]
    print("Design with the Lowest Total Area:")
    print(f"d: {lowest_area_design[0]}")
    print(f"Top Thickness: {lowest_area_design[1]}")
    print(f"Side Thickness: {lowest_area_design[2]}")
    print(f"Sides SF: {lowest_area_design[3]}")
    print(f"Top SF: {lowest_area_design[4]}")
    print(f"Bending SF: {lowest_area_design[5]}")
    print(f"Total Area: {lowest_area_design[6]}")

