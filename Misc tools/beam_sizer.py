import numpy as np

class Material(object):
    def __init__(self, material_name, EModulus=None, YieldStrength=None, Poisson=None, Density=None):
        self.material_name = material_name
        self.EMod = EModulus
        self.Yield = YieldStrength
        self.Poisson = Poisson
        self.Density = Density

class Square_beam(object):
    def __init__(self, length, w=None, h=None, t=None, material=None):
        self.length = length
        self.w = w
        self.h = h
        self.t = t
        self.material = material
        self.EMod = material.EMod
        self.Yield = material.Yield
    def get_Ixx_Iyy_Ixy(self):
        b = self.w - 2*self.t
        d = self.h - 2*self.t
        self.Ixx = (self.w*(self.h**3)/12) - (b*(d**3)/12)
        self.Iyy = (self.h * (self.w ** 3) / 12) - (d * (b ** 3) / 12)
        self.Ixy = 0
    def get_bending_tension(self, moment1, moment2, x, y):
        tension = ((moment1*self.Iyy-moment2*self.Ixy)*y + (moment2*self.Ixx-moment1*self.Ixy)*x)/(self.Ixx*self.Iyy-(self.Ixy**2))
        return tension
    def gettensionSF(self, moment1, moment2, tension=0):
        corner_coordinates = [[-self.w/2, self.h/2],[self.w/2, self.h/2],[self.w/2, -self.h/2],[-self.w/2, -self.h/2]]
        tensions = []
        for coordinates in corner_coordinates:
            tensions.append(self.get_bending_tension(moment1, moment2, *coordinates)+tension)
        tensions = np.array(tensions)
        minSF = self.Yield / np.max(tensions)
        return minSF

    def get_shears(self, shear1, shear2, torsion_shear):
        

    def size_profile(self, shear1, shear2, tension, moment1, moment2):

    def thickness_simple_yield_size(self, loads, SF, initial_thickness, use_self_weight=np.array([0,0,0]),
                                    length_subdivisions=30):
        '''loads are defined as [shear1, shear2, normal tension, position across length], use self weight is an
        acceleration vector [x, y, z] where x is shear 1, y is shear 2 and z is normal tension'''

        if self.w is None or self.h is None or self.material is None or self.EMod is None or self.Yield is None:
            raise ValueError("Not enough fixed parameters, ensure Materials, W, H and EMod are set and the material has a yield")
        self.t = initial_thickness
        self.get_Ixx_Iyy_Ixy()
        last_t = 0
        self.actualSF = 0
        while (self.actualSF/SF) < 1 or (self.actualSF/SF) > 1.1:

