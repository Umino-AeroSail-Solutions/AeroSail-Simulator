import numpy as np
import matplotlib.pyplot as plt
import math
import os


# Example of some added weight list: [[mass1, position], [mass2,position]]

class Segment():
    def __init__(self, length, angle, hold_1_coord, hold_2_coord, mass, addedweightlist, subdivisions=1000, g = 9.81):
        self.length = length
        self.angle = angle
        self.hold_1_coord = hold_1_coord
        self.hold_2_coord = hold_2_coord
        self.weight = mass * g
        self.addedweights = np.array(addedweightlist)
        self.subdivisions = subdivisions
        self.g = g
        vertical_loads = []
        for weight in self.addedweights:
            vertical_loads.append(weight)

        self_weight_per_subdivision = self.weight / self.subdivisions
        x = 0
        deltax = self.length / self.subdivisions
        for i in range(subdivisions):
            x = x + deltax
            vertical_loads.append([self_weight_per_subdivision, x])

        self.vertical_loads = np.array(vertical_loads)
        self.vertical_loads = self.vertical_loads[self.vertical_loads[:, 1].argsort()]
        self.R1 = 0
        self.R2 = 0

    def compute_reaction_loads(self):
        """
        Compute vertical reactions R1 and R2 for a simply-supported beam
        between hold_1_coord and hold_2_coord under vertical loads,
        including overhangs beyond supports.
        """
        # total downward load = sum of all vertical loads
        w = self.vertical_loads[:,0]
        x = self.vertical_loads[:,1]
        self.total_weight = w.sum()

        # moment about hold_1: sum(w_i * (x_i - hold_1_coord))
        distances = x - self.hold_1_coord
        M1 = np.dot(w, distances)

        # span between supports
        span = self.hold_2_coord - self.hold_1_coord
        if span <= 0:
            raise ValueError("Support coordinates invalid: hold_2 must be > hold_1.")

        self.R2 = M1 / span
        self.R1 = self.total_weight - self.R2

        print(f"Total load     : {self.total_weight:.2f} N")
        print(f"Moment about 1 : {M1:.2f} N·m (span = {span:.2f} m)")
        print(f"Reaction R2    : {self.R2:.2f} N")
        print(f"Reaction R1    : {self.R1:.2f} N")
        return self.R1, self.R2

    def compute_internal_loads(self, plot=False):
        self.R1, self.R2 = self.compute_reaction_loads()

        # Add reaction loads to the vertical loads array
        self.vertical_loads = np.vstack([
            self.vertical_loads,
            [-self.R1, self.hold_1_coord],
            [-self.R2, self.hold_2_coord]
        ])

        # Sort by x location to ensure loads are processed in order
        self.vertical_loads = self.vertical_loads[self.vertical_loads[:, 1].argsort()]

        shears = []
        x_locations = []
        total_shear = 0

        for load in self.vertical_loads:
            total_shear += load[0]
            shears.append(total_shear)
            x_locations.append(load[1])

        self.shears = np.array(shears)
        self.x_locations = np.array(x_locations)
        moments = []
        for xi in self.x_locations:
            moment = 0
            for force, xf in self.vertical_loads:
                if xf <= xi:
                    moment += force * (xi - xf)
            moments.append(moment)

        self.moments = np.array(moments)

        # Optional plotting
        if plot:
            # Plot shear diagram
            plt.figure(figsize=(10, 4))
            plt.step(self.x_locations, self.shears, where='post', color='blue')
            plt.title("Shear Force Diagram")
            plt.xlabel("x location (m)")
            plt.ylabel("Shear Force (N)")
            plt.grid(True)
            plt.axhline(0, color='black', linestyle='--')
            plt.show()

            # Plot moment diagram
            plt.figure(figsize=(10, 4))
            plt.plot(self.x_locations, self.moments, color='red')
            plt.title("Bending Moment Diagram")
            plt.xlabel("x location (m)")
            plt.ylabel("Moment (N·m)")
            plt.grid(True)
            plt.axhline(0, color='black', linestyle='--')
            plt.show()


        shearsthreed = []
        momentsthreed = []
        for shear in self.shears:
            sheary = shear*math.sin(math.radians(self.angle))
            shearx = shear*math.cos(math.radians(self.angle))
            shearsthreed.append([shearx, sheary])
        for moment in self.moments:
            momentsy = moment * math.sin(math.radians(self.angle))
            momentsx = moment * math.cos(math.radians(self.angle))
            momentsthreed.append([momentsx, momentsy])

        self.threedmoments = np.array(momentsthreed)
        self.threedshears = np.array(shearsthreed)


# BORUI CHECK THIS NUMBERS

angle = 115.116-90

min_segment_length_difference = 0.5 # Space to allow for ribs and storage
segment_4_length = 10
segment_3_length = segment_4_length - min_segment_length_difference
segment_2_length = segment_3_length - min_segment_length_difference
segment_1_length = segment_2_length - min_segment_length_difference

segment_4_mass = 106.46
segment_3_mass = 217.67
segment_2_mass = 401.70
segment_1_mass = 425.43

aeroplatform_extra_mass_per_segment = 114

#All overlaps are the length of the bottom segment except the first one

OL_34 = segment_3_length
OL_23 = segment_2_length
OL_12 = segment_1_length
OL_01 = 2.6

segment_4_added_weights = [[aeroplatform_extra_mass_per_segment*9.81,segment_4_length]]
Segment4 = Segment(segment_4_length, angle, 0, segment_3_length, segment_4_mass, segment_4_added_weights)
S4R1, S4R2 = Segment4.compute_reaction_loads()
Segment4.compute_internal_loads(plot=True)

segment_3_added_weights = [[aeroplatform_extra_mass_per_segment*9.81,segment_3_length], [S4R1, 0], [S4R2, segment_3_length]]
Segment3 = Segment(segment_3_length, angle, 0, segment_2_length, segment_3_mass, segment_3_added_weights)
S3R1, S3R2 = Segment3.compute_reaction_loads()
Segment3.compute_internal_loads(plot=True)

segment_2_added_weights = [[aeroplatform_extra_mass_per_segment*9.81,segment_2_length], [S3R1, 0], [S3R2, segment_2_length]]
Segment2 = Segment(segment_2_length, angle, 0, segment_1_length, segment_2_mass, segment_2_added_weights)
S2R1, S2R2 = Segment2.compute_reaction_loads()
Segment2.compute_internal_loads(plot=True)

segment_1_added_weights = [[aeroplatform_extra_mass_per_segment*9.81,segment_1_length], [S2R1, 0], [S2R2, segment_1_length]]
Segment1 = Segment(segment_1_length, angle, 0, OL_01, segment_1_mass, segment_1_added_weights)
S1R1, S1R2 = Segment1.compute_reaction_loads()
Segment1.compute_internal_loads(plot=True)

output_folder = "Internal_Loads"
export_internals = "stowed_load_arrays"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

Segment1_shears = Segment1.threedshears
Segment2_shears = Segment2.threedshears
Segment3_shears = Segment3.threedshears
Segment4_shears = Segment4.threedshears

Segment1_moments = Segment1.threedmoments
Segment2_moments = Segment2.threedmoments
Segment3_moments = Segment3.threedmoments
Segment4_moments = Segment4.threedmoments
print("Segment 1 shears: ", Segment1_shears)
print("Segment 1 moments: ", Segment1_moments)

np.savez(os.path.join(output_folder, f"{export_internals}.npz"),
         Segment1_shears=Segment1_shears,
         Segment2_shears=Segment2_shears,
         Segment3_shears=Segment3_shears,
         Segment4_shears=Segment4_shears,
         Segment1_moments=Segment1_moments,
         Segment2_moments=Segment2_moments,
         Segment3_moments=Segment3_moments,
         Segment4_moments=Segment4_moments,
         )