
class Rivet():
    def __init__(self, Shear_Yield_Stress, diameter, Tension_Yield, cost):
        self.Shear_Yield = Shear_Yield_Stress
        self.diameter = diameter
        self.Tension_Yield = Tension_Yield
        self.cost = cost

class Panel():
    def __init__(self, thickness, shear_yield, tension_yield, E_mod, max_bearing_stress):
        self.thickness = thickness
        self.min_rivet_d = 3*self.thickness
        self.max_bearing_stress = max_bearing_stress
        self.E_mod = E_mod
        self.shear_yield = shear_yield
        self.tension_yield = tension_yield

    def Get_Min_Diameter_Bearing(self, force):
        min_diameter_bearing = force / (self.max_bearing_stress*self.thickness) # Bearing min diameter comp https://roymech.org/Useful_Tables/Rivets.html
        return min(min_diameter_bearing, self.min_rivet_d)
    def Get_Max_Total_Hole_D(self, force, width):
        max_D = (force/(self.tension_yield*self.thickness)) + width
        return max_D
class Connection():
    def __init__(self, Tension_stress, Panel1, Panel2, width, SF=4.5): # Normal safety factor is huge due to 3x stress concentrations in holes
        self.SF = SF
        self.Tension_stress = Tension_stress*self.SF
        self.Panel1 = Panel1
        self.Panel2 = Panel2
        self.width = width

