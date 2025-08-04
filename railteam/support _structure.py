import numpy as np
import os
import scipy.integrate as integ
import matplotlib.pyplot as plt

g = 9.80665

safety_factor = 1.5

andres_dic = np.load(os.path.join(os.path.dirname(__file__),"forcevalues.npz"))

R1_values = andres_dic["R1_values"]
R2_values = andres_dic["R2_values"]
l_values = andres_dic["l_values"]
l2_values = andres_dic["l2_values"]
support_forces_bot = andres_dic["support_forces_bot"]
support_forces_top = andres_dic["support_forces_top"]

L_Top = max(l2_values)
L_Bot = max(l_values)

alpha = 0.207026994193 # in radians, from desmos
beta = -0.246726146326 # //


R1_values = [i / 2 * safety_factor for i in R1_values ] # because it's a list.
R2_values = [i / 2 * safety_factor for i in R2_values ] # //
R1_max = np.max(np.abs(np.array(R1_values)))
R2_max = np.max(np.abs(np.array(R2_values)))

support_forces_bot = [i/2*safety_factor for i in support_forces_bot]
support_forces_top = [i/2*safety_factor for i in support_forces_top]
support_bot_1_max, support_bot_1_max_loc = np.max(np.abs(np.array(support_forces_bot[:][0]))), np.argmax(np.abs(np.array(support_forces_bot[:][0]))) # max value, index
support_bot_2_max, support_bot_2_max_loc = np.max(np.abs(np.array(support_forces_bot[:][1]))), np.argmax(np.abs(np.array(support_forces_bot[:][1]))) #
support_top_1_max, support_top_1_max_loc = np.max(support_forces_top[:][0]), np.argmax(support_forces_top[:][0])
support_top_2_max, support_top_2_max_loc = np.max(support_forces_top[:][1]), np.argmax(support_forces_top[:][1])

def get_R_max():
    return [R1_max, R2_max]

# toa compression
def check_buckling(skinside, shearside, sigma, tau, kc, ks, E, v, t):
    
    # skin buckling!!
    crit_skin_buckling = np.pi**2 * kc * E / 12/(1-v**2) * (t/skinside)**2
    print(f"Skin Buckling Safety Mango:\n TOP: {abs(crit_skin_buckling/sigma) - 1}")

    # shear buckling !!
    crit_shear_buckling = np.pi**2 * ks * E / 12/(1-v**2) * (t/shearside)**2
    print(f"Shear Buckling Safety Mangosteen:\n TOP: {abs(crit_shear_buckling/tau) - 1}")



class Material(object):
    def __init__(self, material_name, EModulus=None, YieldStrength=None, Poisson=None, Density=None, max_shear=None):
        self.material_name = material_name
        self.EMod = EModulus
        self.Yield = YieldStrength
        self.Poisson = Poisson
        self.Density = Density
        self.max_shear = max_shear

class Square_beam(object):
    def __init__(self, length, w=None, h=None, t=None, material=None):
        self.length = length
        self.w = w
        self.h = h
        self.t = t
        self.material = material
        self.EMod = material.EMod
        self.Yield = material.Yield
        self.max_shear_yield = material.max_shear
        self.A = w*h - (w-2*t)*(h-2*t) # Area [m²]
        self.m = self.A * self.length * self.material.Density

    def get_Ixx_Iyy_Ixy(self):
        self.b = self.w - 2*self.t
        self.d = self.h - 2*self.t
        self.Ixx = (self.w*(self.h**3)/12) - (self.b*(self.d**3)/12)
        self.Iyy = (self.h * (self.w ** 3) / 12) - (self.d * (self.b ** 3) / 12)
        self.Ixy = 0

    def get_Q(self):
        self.Q = (self.h/2 * self.w)*(self.h/4) - (self.h/2-self.t)*(self.w-2*self.t)*(self.h/4-self.t/2)



check_bar = "toa"

if __name__ == "__main__":
    # Buckling!
    kc = 4 # Buckling coefficient, assume simply supported (conservative)
    ks = 5.5 # Buckling coefficient for shear :D

    if check_bar == "toa":
        steel = Material("mild steel", 207e9, 344e6, 0.3, 7850, (344e6)/2)
        toa = Square_beam(2.2, 100e-3, 100e-3, 10e-3, steel)

        # case for maximum on support_bot
        supp_bot = support_bot_1_max
        supp_top = support_forces_top[support_bot_1_max_loc][0]
        supp_bot_V = supp_bot * np.cos(beta)
        supp_bot_H = supp_bot * np.sin(abs(beta))
        supp_top_V = - supp_top * np.cos(alpha)
        supp_top_H = supp_top * np.sin(alpha)

        toa_base_V = - supp_bot_V - supp_top_V
        toa_base_H = - supp_bot_H - supp_top_H

        toa_axial = np.array([toa_base_V, supp_bot_V, supp_top_V])
        toa_shear = np.array([toa_base_H, supp_bot_H, supp_top_H])
        # print(toa_axial, np.sum(toa_axial), toa_shear, np.sum(toa_shear))
        # do shear and buckling calculations with these values.
        
        # axial loading
        toa_loc = np.array([0, 1.75, 2.1])
        toa_axial = np.vstack((toa_loc, toa_axial))
        # tension-compression graph
        for i, val in enumerate(toa_axial[1]):
            if i != 0:
                toa_tc.append(-val+toa_tc[i-1])
            else: 
                toa_tc = [-val]                
        # just for plotting
        for i, val in enumerate(toa_loc):
            if i != 0:
                toa_loc = np.insert(toa_loc, i+countr, toa_loc[i+countr]-0.0001)

                toa_tc.insert(i+countr, toa_tc[i-1+countr])
                countr += 1
            else:
                countr = 0
        plt.plot(toa_loc, toa_tc)
        plt.show()

        


    # # calculated stresses
    # tau = V * Q / I / t
    # sigma = M * b / 2 / I