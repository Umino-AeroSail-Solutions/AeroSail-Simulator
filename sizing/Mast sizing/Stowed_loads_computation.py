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
        Compute vertical reaction forces at the two pin supports.

        Returns
        -------
        R1, R2 : float
            Vertical reaction forces at support 1 (hold_1_coord) and support 2 (hold_2_coord).
        """
        # First we compute the moment around hold_1 and the total weight

        self.total_weight = np.sum(self.addedweights, axis=0)[0] + self.weight
        M_hold_1 = 0
        for load in self.vertical_loads:
            M_hold_1 += load[0] * (load[1]-self.hold_1_coord)

        # Now this is used to compute R2 and R1
        self.R2 = M_hold_1 / (self.hold_2_coord-self.hold_1_coord)
        self.R1 = self.total_weight - self.R2
        print("R1: ", self.R1)
        print("R2: ", self.R2)
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
            shearx = shear*math.sin(math.radians(self.angle))
            sheary = shear*math.cos(math.radians(self.angle))
            shearsthreed.append([shearx, sheary])
        for moment in self.moments:
            momentsx = moment * math.sin(math.radians(self.angle))
            momentsy = moment * math.cos(math.radians(self.angle))
            momentsthreed.append([momentsx, momentsy])

        self.threedmoments = np.array(momentsthreed)
        self.threedshears = np.array(shearsthreed)


# BORUI CHECK THIS NUMBERS

angle = /Look it up LOL Borui I am in Chimera Worksession/

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

Segment4 = Segment(segment_4_length,)