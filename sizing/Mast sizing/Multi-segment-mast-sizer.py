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
        self.bottom_overlap_top_force = (self.top_overlap_top_force * self.length) + self.applied_force_vector * (0.5*(self.length + self.bottom_overlap)) + self.top_overlap_bottom_force * (self.length - self.top_overlap)
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

        if plot:
            positions = np.array([0, self.bottom_overlap, self.bottom_overlap, self.length-((self.length-self.bottom_overlap)*0.5), self.length-((self.length-self.bottom_overlap)*0.5), self.length-self.top_overlap, self.length-self.top_overlap, self.length])
            shear_mags = np.zeros(8)
            i=0
            for shear in shears:
                shear_mags[i] = np.linalg.norm(shear)
                i+=1
            plt.plot(positions, shear_mags, marker='o')
            plt.fill_between(positions, shear_mags, alpha=0.3)
            plt.xlabel('Position')
            plt.ylabel('Shear Force')
            plt.title('Shear Force Distribution')
            plt.grid(True)
            for idx, (pos, shear_mag) in enumerate(zip(positions, shear_mags)):
                plt.annotate(f'{idx}', xy=(pos, shear_mag), textcoords='offset points', xytext=(0,5), ha='center')
            plt.show()
        if np.linalg.norm(shears[0]-self.bottom_overlap_bottom_force) > 0.1:
            error = "There is no shear equilibrium in the segment :("
            return error


# Example usage
segment = Segment(100, 10, 2, 0, np.array([0,0]), np.array([0,0]))
segment.compute_internal_loads(plot=True)
