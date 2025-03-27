import numpy as np
import Cross_section_analysis as cs
import matplotlib.pyplot as plt
from tqdm import tqdm

class Segment():
    def __init__(self, added_weight, length, bottom_overlap, top_overlap, top_overlap_top_force,
                 top_overlap_bottom_force, total_force_vector, total_height):
        self.added_weight = added_weight
        self.length = length
        self.bottom_overlap = bottom_overlap
        self.top_overlap = top_overlap
        self.top_overlap_top_force = top_overlap_top_force
        self.top_overlap_bottom_force = top_overlap_bottom_force
        self.applied_force_vector = ((length-bottom_overlap)/total_height)*total_force_vector
        self.width = None
        self.height = None
        self.positions = np.array(
            [0, self.bottom_overlap, self.bottom_overlap, self.length - ((self.length - self.bottom_overlap) * 0.5),
             self.length - ((self.length - self.bottom_overlap) * 0.5), self.length - self.top_overlap,
             self.length - self.top_overlap, self.length])
        self.shears = None
        self.moments = None
    def reset_parameters(self, added_weight=None, length=None, bottom_overlap=None, top_overlap=None,
                        top_overlap_top_force=None, top_overlap_bottom_force=None, total_force_vector=None,
                        total_height=None):
        if added_weight is not None:
            self.added_weight = added_weight
        if length is not None:
            self.length = length
        if bottom_overlap is not None:
            self.bottom_overlap = bottom_overlap
        if top_overlap is not None:
            self.top_overlap = top_overlap
        if top_overlap_top_force is not None:
            self.top_overlap_top_force = top_overlap_top_force
        if top_overlap_bottom_force is not None:
            self.top_overlap_bottom_force = top_overlap_bottom_force
        if total_force_vector is not None:
            self.applied_force_vector = ((self.length - self.bottom_overlap) / (
                total_height if total_height is not None else height)) * total_force_vector
    def compute_reaction_loads(self):
        self.bottom_overlap_top_force = ((self.top_overlap_top_force * self.length) + self.applied_force_vector * (self.length - ((self.length - self.bottom_overlap) * 0.5)) - self.top_overlap_bottom_force * (self.length - self.top_overlap))/self.bottom_overlap
        self.bottom_overlap_bottom_force = self.bottom_overlap_top_force + self.top_overlap_bottom_force - self.applied_force_vector - self.top_overlap_top_force
        return self.bottom_overlap_top_force, self.bottom_overlap_bottom_force
    def compute_internal_loads(self, plot=False):
        self.compute_reaction_loads()
        # Let's first do the shears, they go from 0 to 7 at both sides of the interest points
        shears = np.zeros((8,2))
        shears[7] = 0-self.top_overlap_top_force
        shears[6] = shears[7]
        shears[5] = shears[6]+self.top_overlap_bottom_force
        shears[4] = shears[5]
        shears[3] = shears[4]-self.applied_force_vector
        shears[2] = shears[3]
        shears[1] = shears[2]+self.bottom_overlap_top_force
        shears[0] = shears[1]

        self.shears = shears

        if np.linalg.norm(shears[0]-self.bottom_overlap_bottom_force) > 0.1:
            error = "There is no shear equilibrium in the segment :("
            print(error)

        # Now we do the internal Moments
        # Ik there is a more efficient way to write this but this is as fast when running
        moments = np.zeros((8, 2))
        moments[7] = 0
        moments[6] = self.top_overlap_top_force * (self.length - self.positions[6])
        moments[5] = moments[6]
        moments[4] = self.top_overlap_top_force * (self.length - self.positions[4]) - self.top_overlap_bottom_force * (self.positions[5]-self.positions[4])
        moments[3] = moments[4]
        moments[2] = self.top_overlap_top_force * (self.length - self.positions[2]) - self.top_overlap_bottom_force * (self.positions[5]-self.positions[2]) + self.applied_force_vector * (self.positions[4]-self.positions[2])
        moments[1] = moments[2]
        moments[0] = self.top_overlap_top_force * self.length - self.top_overlap_bottom_force * self.positions[5]  + self.applied_force_vector * self.positions[3] - self.bottom_overlap_top_force * self.positions[1]
        self.moments = moments
        if np.linalg.norm(moments[0]) > 0.1:
            error = "There is no moment equilibrium in the segment :("
            print(error)
        # And now we plot

        if plot:
            fig, ax1 = plt.subplots()

            color = 'tab:blue'
            ax1.set_xlabel('Position')
            ax1.set_ylabel('Shear Force', color=color)
            # Plots x coordinates
            ax1.plot(self.positions, shears[:,0], marker='o', color=color)
            ax1.fill_between(self.positions, shears[:,0], color=color, alpha=0.3)
            ax1.tick_params(axis='y', labelcolor=color)
            for idx, (pos, shear) in enumerate(zip(self.positions, shears[:,0])):
                ax1.annotate(f'{idx}', xy=(pos, shear), textcoords='offset points', xytext=(0, 5), ha='center')

            ax2 = ax1.twinx()
            color = 'tab:red'
            ax2.set_ylabel('Moment', color=color)
            # Plots x coordinates
            ax2.plot(self.positions, moments[:,0], marker='x', linestyle='--', color=color)
            ax2.fill_between(self.positions, moments[:,0], color=color, alpha=0.1)
            ax2.tick_params(axis='y', labelcolor=color)
            for idx, (pos, moment) in enumerate(zip(self.positions,moments[:,0])):
                ax2.annotate(f'{idx}', xy=(pos, moment), textcoords='offset points', xytext=(0, 5), ha='center')

            # Aligning the y=0 axis for both shears and moments
            shear_min, shear_max = ax1.get_ylim()
            moment_min, moment_max = ax2.get_ylim()
            combined_min = min(shear_min, moment_min)
            combined_max = max(shear_max, moment_max)
            ax1.set_ylim(combined_min, combined_max)
            ax2.set_ylim(combined_min, combined_max)

            fig.tight_layout()
            plt.title('Shear Force and Moment Distribution')
            plt.grid(True)
            plt.show()
    def set_width_height(self, width, height):
        self.width = width
        self.height = height

    def optimize_cross_section(self, Vx, Vy, Mx, My, material_density, L, added_weight, v_slots,subdivisions=15,
                               skin_step_thickness=0.0005, min_skin_thickness=0.0005, max_thickness=0.01, Print=False):
        #To avoid dividing by zero
        if Mx == 0:
            Mx += 0.01
        if My == 0:
            My += 0.01
        if Vx == 0:
            Vx += 0.01
        if Vy == 0:
            Vy += 0.01

        w = self.width
        h = self.height
        total_area = 0
        Nz = 0
        prev_Nz = 2
        possible_designs = []
        while abs(prev_Nz - Nz) > .1:
            prev_Nz = Nz
            Nz = total_area * material_density * L + added_weight
            Vx, Vy, Mx, My, Nz = Vx * SF, Vy * SF, Mx * SF, My * SF, Nz * SF

            for v_slot in v_slots:
                a = v_slot[0]
                d = v_slot[1]
                cornerIxx = v_slot[2]
                for top_thickness in np.arange(min_skin_thickness, max_thickness, skin_step_thickness):
                    for side_thickness in np.arange(min_skin_thickness, max_thickness, skin_step_thickness):
                        areas, thicknesses = cs.create_areas_and_thicknesses(w, h, d, a, subdivisions, top_thickness,
                                                                             side_thickness)
                        Ixx, Iyy, Ixy = cs.compute_Ixx_Iyy_Ixy(w, h, top_thickness, side_thickness, d, a,
                                                               cornerIxx=cornerIxx)
                        total_area = 2 * (top_thickness * w + side_thickness * h) + 4 * a
                        bending_ok, bending_sf = cs.check_bending_ok(w, h, top_thickness, side_thickness, d, a, Mx, My,
                                                                     max_tension, cornerIxx=cornerIxx, Nz=Nz,
                                                                     area=total_area)
                        cs.update_areas_bending(areas, thicknesses, Mx, My, Ixx, Iyy, Ixy)
                        qb, qbmax = cs.compute_shear_flows(areas, Vx, Vy, subdivisions, thicknesses, h, d, w)
                        top_ok, sides_ok, top_SF, sides_SF = cs.check_shear_ok(qb, top_thickness, side_thickness,
                                                                               max_shear, thicknesses)
                        if top_ok and sides_ok and bending_ok:
                            bending_sf = bending_sf * SF
                            top_SF = top_SF * SF
                            sides_SF = sides_SF * SF
                            total_area = 2 * (top_thickness * w + side_thickness * h) + 4 * a
                            possible_designs.append(
                                [d, top_thickness, side_thickness, sides_SF, top_SF, bending_sf, total_area])
                            break
        possible_designs = np.array(possible_designs)
        # print(possible_designs)
        if np.size(possible_designs) == 0:
            print("There is no possible designs available :(")
            print("Tip: Increase the thickness limit or add bigger vslot channels")
            print("Tip 2: Changing weapons is faster than reloading.......")
            return None
        else:
            # Sort the possible designs by total area (last column)
            sorted_designs = possible_designs[possible_designs[:, -1].argsort()]
            # Output the design with the lowest total area
            lowest_area_design = sorted_designs[0]
            if Print:
                print("Design with the Lowest Total Area:")
                print(f"d: {lowest_area_design[0]}")
                print(f"Top Thickness: {lowest_area_design[1]}")
                print(f"Side Thickness: {lowest_area_design[2]}")
                print(f"Sides SF: {lowest_area_design[3]}")
                print(f"Top SF: {lowest_area_design[4]}")
                print(f"Bending SF: {lowest_area_design[5]}")
                print(f"Total Area: {lowest_area_design[6]}")
                print(f"Total mass: {(lowest_area_design[6] * self.length * material_density)}")
            return lowest_area_design
    def size_it(self, material_density,max_tension, max_shear, v_slots, skin_step_thickness=0.0005,
                min_skin_thickness=0.0005, max_thickness=0.01, Print=False):
        self.compute_internal_loads()
        cross_sections = []
        for idx in range(self.shears.shape[0]):
            # print("Analyizing position: ", idx)
            cross_sections.append(self.optimize_cross_section(self.shears[idx,0], self.shears[idx,1], self.moments[idx,0], self.moments[idx,1], material_density, self.length-self.positions[idx], self.added_weight, v_slots, skin_step_thickness=skin_step_thickness, min_skin_thickness=min_skin_thickness, max_thickness=max_thickness, Print=Print))
        # Now we choose the final cross-section based on the maximum of every parameter
        cross_sections = np.array(cross_sections)

        # [d, top_thickness, side_thickness, sides_SF, top_SF, bending_sf,total_area]
        final_cross_section = [max(cross_sections[:,0]), max(cross_sections[:,1]), max(cross_sections[:,2]), None, None, None, None]
        a = 0
        # V-slot options [[area, d, Ixx/Iyy], [...,...]]
        used_v_slot = np.zeros_like(v_slots[0])
        for v_slot in v_slots:
            if v_slot[1] == final_cross_section[0]:
                a = v_slot[0]
                used_v_slot = v_slot
                # print("Chosen v-slot: ")
                # print(used_v_slot)

        final_cross_section[-1] = 2 * (final_cross_section[1] * self.width + final_cross_section[2] * self.height) + 4 * a

        # and now we should find the actual safety factors with the updated crossection applied but i'm finished for today <-- lazy ass bastard

        w = self.width
        h = self.height
        d = used_v_slot[1]
        subdivisions = 100
        top_thickness = final_cross_section[1]
        side_thickness = final_cross_section[2]
        cornerIxx = used_v_slot[2]

        safety_factors = []

        for idx in range(self.shears.shape[0]):
            Nz = self.added_weight + (self.length - self.positions[idx]) * final_cross_section[-1] * material_density * 9.81
            Mx = self.moments[idx,0]
            My = self.moments[idx,1]
            Vx = self.shears[idx,0]
            Vy = self.shears[idx,1]
            # To avoid dividing by zero
            if abs(Mx) < 0.1:
                Mx += 0.1
            if abs(My) < 0.1:
                My += 0.1
            if abs(Vx) < 0.1:
                Vx += 0.1
            if abs(Vy) < 0.1:
                Vy += 0.1

            areas, thicknesses = cs.create_areas_and_thicknesses(w, h, d, a, subdivisions, top_thickness,
                                                                 side_thickness)
            Ixx, Iyy, Ixy = cs.compute_Ixx_Iyy_Ixy(w, h, top_thickness, side_thickness, d, a, cornerIxx=cornerIxx)
            total_area = 2 * (top_thickness * w + side_thickness * h) + 4 * a
            bending_ok, bending_sf = cs.check_bending_ok(w, h, top_thickness, side_thickness, d, a, Mx, My, max_tension,
                                                         cornerIxx=cornerIxx, Nz=Nz, area=total_area)
            cs.update_areas_bending(areas, thicknesses, Mx, My, Ixx, Iyy, Ixy)
            qb, qbmax = cs.compute_shear_flows(areas, Vx, Vy, subdivisions, thicknesses, h, d, w)
            top_ok, sides_ok, top_SF, sides_SF = cs.check_shear_ok(qb, top_thickness, side_thickness, max_shear,
                                                                   thicknesses)

            safety_factors.append([bending_sf, top_SF, sides_SF])

        safety_factors = np.array(safety_factors)
        min_safety_factors = np.array([np.min(safety_factors[:,0]), np.min(safety_factors[:,1]), np.min(safety_factors[:,2])])

        final_cross_section[3] = min_safety_factors[2]
        final_cross_section[4] = min_safety_factors[1]
        final_cross_section[5] = min_safety_factors[0]
        if Print:
            print('\n Final design:')
        labels = [
            "V-slot side size",
            "Top Thickness",
            "Side Thickness",
            "Side Safety Factor",
            "Top Safety Factor",
            "Bending Safety Factor",
            "Total Area"
        ]
        if Print:
            for label, value in zip(labels, final_cross_section):
                print(f"{label}: {value:.6f}")
        total_mass = (self.length * final_cross_section[-1] * material_density)
        if Print:
            print("Total mass: ", (self.length * final_cross_section[-1] * material_density))
        return total_mass, final_cross_section
    # def optimize_bottom_overlap(self): # Does not work
    #     self.compute_internal_loads()
    #     top_moment = np.linalg.norm(self.moments[5])
    #     bot_moment = np.linalg.norm(self.moments[2])
    #     ratio = bot_moment / top_moment
    #     while ratio < 0.9 or ratio > 1.1:
    #         old_bottom_overlap = self.bottom_overlap
    #         self.bottom_overlap = self.bottom_overlap * ratio
    #         self.top_overlap = self.top_overlap - (old_bottom_overlap - self.bottom_overlap)
    #         self.compute_internal_loads()
    #         top_moment = np.linalg.norm(self.moments[5])
    #         bot_moment = np.linalg.norm(self.moments[2])
    #         ratio = bot_moment / top_moment
    #     print("Bottom overlap: ", self.bottom_overlap)
    #     print("Ratio: ", ratio)

# Example usage -----------------------------

# V-slot options [[area, d, Ixx/Iyy], [...,...]]-------
# t-slot 2020, t-slot 4040, t-slot 8080 (quad 4040 version)---------
# v_slot_options = np.array([[.00016210, 0.02, .0000000066881], [0.000504, 0.04, 0.000000073272], [.0016122, 0.08, .0000010225207]])
# segment = Segment(100, 10, 2, 0, np.array([0,0]), np.array([0,0]))
# segment.compute_internal_loads(plot=True)
# segment.set_width_height(bottom_w, bottom_h)
# aludenisy= 2710
# segment.size_it(aludenisy,  max_tension, max_shear, v_slot_options)

# Actual Multi-segment mast sizer -------------------------------------------------------------------------------------

# V-slot options [[area, d, Ixx/Iyy], [...,...]]-------
# t-slot 2020, t-slot 4040, t-slot 8080 (quad 4040 version)---------

SF = 3

max_tension = 200000000.0 # https://struxure.com/wp-content/uploads/2020/05/01-_-6063-T6-Aluminum-MSDS.pdf
max_shear = 152000000.0

bottom_w, bottom_h = 1.44, 0.7

windspeed_kt = 30
windspeed = windspeed_kt/1.944
density = 1.22522568
cf = 0.49
height = 30
chord = 5
surface_area = height*chord

F_applied_magnitude = 0.5*density*(windspeed**2) * surface_area * cf * SF
print("Force applied: ", str(F_applied_magnitude), "N")

lod_assumed = 10
angle = np.arctan(lod_assumed)
print("angle: ", np.degrees(angle), "degree")
total_force_vector = np.array([F_applied_magnitude*np.cos(angle), F_applied_magnitude*np.sin(angle)])
print("total force vector: ", total_force_vector)

v_slot_options = np.array([[.00016210, 0.02, .0000000066881], [0.000504, 0.04, 0.000000073272], [.0016122, 0.08, .0000010225207]])

aludenisy = 2710

max_segment_length = 10
min_segment_length_difference = 0.5 # Space to allow for ribs and storage

# First guess

segment_4_length = 10
segment_3_length = segment_4_length - min_segment_length_difference
segment_2_length = segment_3_length - min_segment_length_difference
segment_1_length = segment_2_length - min_segment_length_difference

overlap_01 = 2.6
total_overlap_possible = (segment_4_length +segment_3_length +segment_2_length +segment_1_length -overlap_01) - height
print("total overlap possible: ", total_overlap_possible)

# Equal overlaps:

# overlap_34 = total_overlap_possible/3
# overlap_23 = total_overlap_possible/3
# overlap_12 = total_overlap_possible/3

# Different overlaps overwrite:

overlap_34 = 1.03
overlap_23 = 1.40
overlap_12 = total_overlap_possible - overlap_34 - overlap_23
print()
print("overlap_34: ", overlap_34)
print("overlap_23: ", overlap_23)
print("overlap_12: ", overlap_12)
print("overlap_01: ", overlap_01)
print()
print()

w1, h1, d1 = bottom_w, bottom_h, 0.08
w2, h2, d2 = bottom_w-2*d1, bottom_h-2*d1, 0.08
w3, h3, d3 = w2-2*d2, h2-2*d2, 0.04
w4, h4, d4 = w3-2*d3, h3-2*d3, 0.02

def get_possible_overlap():
    return (segment_4_length +segment_3_length +segment_2_length +segment_1_length -overlap_01 -overlap_12 -overlap_23 -overlap_34) - height
# print("possible overlap: ", get_possible_overlap())

sail_weight_4 = (456/4)*9.81 # Mass numbers for the aero platform subsystem
sail_weight_3 = (456/4)*9.81
sail_weight_2 = (456/4)*9.81
sail_weight_1 = (456/4)*9.81

def get_bool(prompt):
    while True:
        try:
            return {"y":True,"n":False}[input(prompt).lower()]
        except KeyError:
            print("Invalid input please enter y or n!")
            return False

compute_cross_section = get_bool("Compute cross section? [y/N]")
def optimize_mast(Print=False,plot=False):
    segment_4 = Segment(sail_weight_4, segment_4_length, overlap_34, 0, np.array([0,0]), np.array([0,0]), total_force_vector, height)
    segment_4.compute_internal_loads(plot)
    segment_3_top_force_top, segment_3_top_force_bottom = segment_4.compute_reaction_loads()

    segment_4_mass = 0

    if compute_cross_section:
        if Print:
            print()
            print("Segment 4 design: \n")

        segment_4.set_width_height(w4, h4)

        segment_4_mass, segment_4_design = segment_4.size_it(aludenisy, max_tension, max_shear, v_slot_options, Print=Print)
        # pbar.update(1) # Pbar option for once every section
    added_weight = sail_weight_4*2 + sail_weight_3 + 9.81*segment_4_mass*2 # Weight from segment 4 is doubled due to the pulley effects
    segment_3 = Segment(added_weight, segment_3_length, overlap_23, overlap_34, segment_3_top_force_top, segment_3_top_force_bottom,total_force_vector, height)
    # segment_3.optimize_bottom_overlap()
    segment_3.compute_internal_loads(plot)
    segment_2_top_force_top, segment_2_top_force_bottom = segment_3.compute_reaction_loads()

    segment_3_mass = 0

    if compute_cross_section:
        if Print:
            print()
            print("Segment 3 design: \n")

        segment_3.set_width_height(w3, h3)

        segment_3_mass, segment_3_design = segment_3.size_it(aludenisy, max_tension, max_shear, v_slot_options, Print=Print)
        # pbar.update(1) # Pbar option for once every section
    added_weight = sail_weight_4 + sail_weight_3*2 + 9.81*segment_4_mass + sail_weight_2 + 9.81*segment_3_mass*2 # Check pulley effect
    segment_2 = Segment(added_weight, segment_2_length, overlap_12, overlap_23, segment_2_top_force_top, segment_2_top_force_bottom,total_force_vector, height)
    # segment_2.optimize_bottom_overlap()
    segment_2.compute_internal_loads(plot)
    segment_1_top_force_top, segment_1_top_force_bottom = segment_2.compute_reaction_loads()

    segment_2_mass = 0

    if compute_cross_section:
        if Print:
            print()
            print("Segment 2 design: \n")

        segment_2.set_width_height(w2, h2)

        segment_2_mass, segment_2_design = segment_2.size_it(aludenisy, max_tension, max_shear, v_slot_options, Print=Print)
        # pbar.update(1) # Pbar option for once every section
    added_weight = sail_weight_4 + sail_weight_3 + 9.81*segment_4_mass + sail_weight_2*2 + 9.81*segment_3_mass + sail_weight_1 + 9.91*segment_2_mass*2 # Check pulley effect
    segment_1 = Segment(added_weight, segment_1_length, overlap_01, overlap_12, segment_1_top_force_top, segment_1_top_force_bottom,total_force_vector, height)
    # segment_1.optimize_bottom_overlap()
    segment_1.compute_internal_loads(plot)
    reaction_force_top, reaction_force_bottom = segment_1.compute_reaction_loads()

    segment_1_mass = 0

    if compute_cross_section:
        if Print:
            print()
            print("Segment 1 design: \n")

        segment_1.set_width_height(w1, h1)

        segment_1_mass, segment_1_design = segment_1.size_it(aludenisy, max_tension, max_shear, v_slot_options, Print=Print)
    if compute_cross_section and Print:
        print("\nTotal mast mass: ", (segment_1_mass + segment_2_mass + segment_3_mass + segment_4_mass))
        print("Total sail planform mass: ", ((sail_weight_1 + sail_weight_2 + sail_weight_3 + sail_weight_4))/9.81)

        print("Total sail mass: ", ((segment_1_mass + segment_2_mass + segment_3_mass + segment_4_mass)+((sail_weight_1 + sail_weight_2 + sail_weight_3 + sail_weight_4))/9.81))
    return ((segment_1_mass + segment_2_mass + segment_3_mass + segment_4_mass)+((sail_weight_1 + sail_weight_2 + sail_weight_3 + sail_weight_4))/9.81)

# Single mast size

optimize_mast(Print=True, plot=True)

# # Seeding code
#
# overlap_34_range = [0.9, 1.5]
# overlap_23_range = [1.2, 1.6] # Refined from multiple unfinished runs
#
# ol_34_iterations = 4
# ol_23_iterations = 4
#
# seed_deepness = 3
#
# inseed_zoom = 0.75  # how much around the point does it go
#
# optimum_overlaps = []
# optimum_mass = 4000
#
# # Collect seeding locations for plotting
# seeding_locations = []
#
# for i in range(seed_deepness):
#     overlaps_34 = np.linspace(overlap_34_range[0], overlap_34_range[1], ol_34_iterations)
#     overlaps_23 = np.linspace(overlap_23_range[0], overlap_23_range[1], ol_23_iterations)
#     seed_spacing_34 = (overlap_34_range[1] - overlap_34_range[0])/ol_34_iterations
#     seed_spacing_23 = (overlap_23_range[1] - overlap_23_range[0])/ol_23_iterations
#     # total_iterations = len(overlaps_34) * len(overlaps_23) # Once for every mast
#     total_iterations = len(overlaps_34) * len(overlaps_23) * 4 # Once for every segment
#     # total_iterations = len(overlaps_34) * len(overlaps_23) * 4 * 8 # Once for every cross section
#     with tqdm(total=total_iterations, desc=f'Seed Iteration {i + 1}/{seed_deepness}') as pbar:
#         for overlap_34 in overlaps_34:
#             for overlap_23 in overlaps_23:
#                 overlap_12 = total_overlap_possible - overlap_34 - overlap_23
#                 mass = optimize_mast(Print=True)
#                 # pbar.update(1) # Pbar option for once every mast
#                 print("Run completed. Mass:", mass)
#                 if mass < optimum_mass:
#                     optimum_mass = mass
#                     optimum_overlaps = [overlap_23, overlap_34]
#                 seeding_locations.append((overlap_23, overlap_34, mass))  # Store the seeding location
#                 # Plotting the seeding locations
#                 seeding_locations_array = np.array(seeding_locations)
#
#                 fig, ax = plt.subplots(figsize=(10, 6))
#
#                 # Create a scatter plot with color grading based on the third dimension (e.g., seeding_locations_array[:, 2])
#                 scatter = ax.scatter(seeding_locations_array[:, 0], seeding_locations_array[:, 1],
#                                      c=seeding_locations_array[:, 2], cmap='viridis', marker='o',
#                                      label='Seeding Locations')
#
#                 # Add a colorbar to show the color grading
#                 colorbar = plt.colorbar(scatter, ax=ax)
#                 colorbar.set_label('Mass')
#
#                 # Plot the optimal point in red
#                 ax.scatter(optimum_overlaps[0], optimum_overlaps[1], c='red', marker='x', s=50, label='Optimal Point')
#
#                 for point in seeding_locations:
#                     ax.text(point[0], point[1], f'{str(int(point[2]))}', fontsize=8, ha='right')
#
#                 ax.set_xlabel('Overlap 23')
#                 ax.set_ylabel('Overlap 34')
#                 ax.set_title('Seeding Locations for Overlap Optimization')
#                 ax.legend()
#                 ax.grid(True)
#
#                 plt.show()
#
#     overlap_34_range = [optimum_overlaps[1]-(seed_spacing_34*inseed_zoom),
#                         optimum_overlaps[1]+(seed_spacing_34*inseed_zoom)]
#     overlap_23_range = [optimum_overlaps[0] - (seed_spacing_23 * inseed_zoom),
#                         optimum_overlaps[0] + (seed_spacing_23 * inseed_zoom)]
# print("Optimum overlap is: ", optimum_overlaps)
#
# # Optimum overlap is:  [np.float64(1.3989583333333333), np.float64(1.0296875)]