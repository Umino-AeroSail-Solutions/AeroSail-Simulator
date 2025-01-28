import numpy as np
import Cross_section_analysis as cs
import matplotlib.pyplot as plt

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

class Segment():
    def __init__(self, added_weight, length, bottom_overlap, top_overlap, top_overlap_top_force,
                 top_overlap_bottom_force, total_force_vector=total_force_vector, total_height=height):
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

    def optimize_cross_section(self, Vx, Vy, Mx, My, material_density, L, added_weight, v_slots,subdivisions=100,
                               skin_step_thickness=0.0005, min_skin_thickness=0.0005, max_thickness=0.01):
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
            print("\nTip: Changing weapons is faster than reloading.......")
            return None
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
            print(f"Total mass: {(lowest_area_design[6] * self.length * material_density)}")
            return lowest_area_design
    def size_it(self, material_density,max_tension, max_shear, v_slots, skin_step_thickness=0.0005,
                min_skin_thickness=0.0005, max_thickness=0.01):
        self.compute_internal_loads()
        cross_sections = []
        for idx in range(self.shears.shape[0]):
            print("\n Analyizing position: ", idx)
            cross_sections.append(self.optimize_cross_section(self.shears[idx,0], self.shears[idx,1], self.moments[idx,0], self.moments[idx,1], material_density, self.length-self.positions[idx], self.added_weight, v_slots, skin_step_thickness=skin_step_thickness, min_skin_thickness=min_skin_thickness, max_thickness=max_thickness))
        # Now we choose the final cross-section based on the maximum of every parameter
        cross_sections = np.array(cross_sections)

        # [d, top_thickness, side_thickness, sides_SF, top_SF, bending_sf,total_area]
        final_cross_section = [max(cross_sections[0,:]), max(cross_sections[1,:]), max(cross_sections[2,:]), None, None, None, None]
        a = 0

        for v_slot in v_slots:
            if v_slot[1] == final_cross_section[0]:
                a = v_slot[0]

        final_cross_section[-1] = 2 * (final_cross_section[1] * self.width + final_cross_section[2] * self.height) + 4 * a

        # and now we should find the actual safety factors with the updated crossection applied but i'm finished for today

# Example usage

# V-slot options [[area, d, Ixx/Iyy], [...,...]]
# t-slot 2020, t-slot 4040, t-slot 8080 (quad 4040 version)
v_slot_options = np.array([[.00016210, 0.02, .0000000066881], [0.000504, 0.04, 0.000000073272], [.0016122, 0.08, .0000010225207]])
segment = Segment(100, 10, 2, 0, np.array([0,0]), np.array([0,0]))
segment.compute_internal_loads(plot=True)
segment.set_width_height(bottom_w, bottom_h)
aludenisy= 2710
segment.size_it(aludenisy,  max_tension, max_shear, v_slot_options)