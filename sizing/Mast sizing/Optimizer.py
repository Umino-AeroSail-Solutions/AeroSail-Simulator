import numpy as np
import Cross_section_analysis as cs
# Finds suitable design with the least possible area

# V-slot options [[area, d, Ixx/Iyy], [...,...]]
# t-slot 2020, t-slot 4040, t-slot 8080 (quad 4040 version)
v_slots = np.array([[.00016210, 0.02, .0000000066881], [0.000504, 0.04, 0.000000073272], [.0016122, 0.08, .0000010225207]])
# v_slots = np.array([[.0016122, 0.08, .0000010225207]])
# w, h = 0.288, 0.14 # I'M FUCKING RETARDED IM SO FUCKING RETARDED
w, h = 1.44, 0.7 # Proper values

SF = 3

max_tension = 200000000.0 # https://struxure.com/wp-content/uploads/2020/05/01-_-6063-T6-Aluminum-MSDS.pdf
max_shear = 152000000.0

windspeed_kt = 30
windspeed = windspeed_kt/1.944
density = 1.22522568
cf = 0.49
height = 30
chord = 5
surface_area = height*chord

F_applied_magnitude = 0.5*density*(windspeed**2) * surface_area * cf
print(F_applied_magnitude)
# print(F_applied_magnitude)
# Shear computations --> see NVM diagrams

G = height/2
A_top = 2
A_bot = 1
L = 9

F_1 = F_applied_magnitude * (((G+L-A_top)/A_bot)-1)
F_2 = F_applied_magnitude * ((G+L-A_top)/A_bot)
F_3 = F_applied_magnitude * (G/A_top)

shear_1 = -1*F_1
shear_2 = shear_1 + F_2
shear_3 = shear_2 + F_3

moment_1 = shear_1 * A_bot
moment_2 = shear_2 * A_top

shears = [shear_1, shear_2, shear_3]
moments = [moment_1, moment_2]
lod_assumed = 10
angle = np.arctan(lod_assumed)

added_weight = 1000*9.81 # Asumes 10kg of extra weight above the first mast segment

def optimize(Vx, Vy, Mx, My, material_density, L, segmentheight=L):
    total_area = 0
    Nz = 0
    prev_Nz = 2
    possible_designs = []
    while abs(prev_Nz-Nz) > .1 :
        prev_Nz  = Nz
        Nz = total_area * material_density *L + added_weight
        Vx, Vy, Mx, My, Nz = Vx * SF, Vy * SF, Mx * SF, My * SF, Nz * SF

        skin_step_thickness = 0.0005
        min_skin_thickness = 0.0005
        max_thickness = 0.01


        subdivisions = 100

        for v_slot in v_slots:
            a = v_slot[0]
            d = v_slot[1]
            cornerIxx = v_slot[2]
            for top_thickness in np.arange(min_skin_thickness, max_thickness, skin_step_thickness):
                for side_thickness in np.arange(min_skin_thickness, max_thickness, skin_step_thickness):
                    areas, thicknesses = cs.create_areas_and_thicknesses(w, h, d, a, subdivisions, top_thickness, side_thickness)
                    Ixx, Iyy, Ixy = cs.compute_Ixx_Iyy_Ixy(w, h, top_thickness, side_thickness, d, a, cornerIxx=cornerIxx)
                    total_area = 2 * (top_thickness * w + side_thickness * h) + 4 * a
                    bending_ok, bending_sf = cs.check_bending_ok(w, h, top_thickness, side_thickness, d, a, Mx, My, max_tension, cornerIxx=cornerIxx, Nz=Nz, area=total_area)
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
        print(f"Total mass: {(lowest_area_design[6]*L*material_density)}")

aludenisy= 2710
# Case 0

print("\nCase 1:\n")
Mx, My = 0.001, 0.001 # Added to avoid division by zero in safety
Vx, Vy = shears[0]*np.cos(angle), shears[0]*np.sin(angle)
optimize(Vx, Vy, Mx, My, aludenisy, (L))

# Case 1
print("\nCase 2:\n")
Vx, Vy = shears[0]*np.cos(angle), shears[0]*np.sin(angle)
Mx, My = moments[0]*np.sin(angle), -1*moments[0]*np.cos(angle)

optimize(Vx, Vy, Mx, My, aludenisy, (L-A_bot))

# Case 2
print("\nCase 3:\n")

Vx, Vy = shears[1]*np.cos(angle), shears[1]*np.sin(angle)
Mx, My = moments[1]*np.sin(angle), -1*moments[1]*np.cos(angle)
optimize(Vx, Vy, Mx, My, aludenisy, (A_top))

# Case 3
print("\nCase 4:\n")
Vx, Vy = shears[2]*np.cos(angle), shears[2]*np.sin(angle)
Mx, My = moments[1]*np.sin(angle), -1*moments[1]*np.cos(angle)
optimize(Vx, Vy, Mx, My, aludenisy, (A_top))