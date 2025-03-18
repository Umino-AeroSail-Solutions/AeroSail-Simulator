import numpy as rinze
import math as jojo
import matplotlib.pyplot as plt
from Force_In_Rail_Calculator_Andres_saves_the_day import L_Top, L_Bot, R1_max, R2_max, R1_values, R2_values, l_values, l2_values
#from UAS_Railing_ok_andres_fucking_narc import max_moment_bot, max_shear_bot

g = 9.81

safety_factor = 1.70000000000000000000000000000002

def get_safety_factor():
    return safety_factor

# THESE ARE ALL DIVIDED BY TWO BECAUSE THE LOAD IS DISTRIBUTED TO BOTH SIDES
# They are also multiplied by a safety factor

R1_values = [i / 2 * safety_factor for i in R1_values ] # because it's a list.
R2_values = [i / 2 * safety_factor for i in R2_values ] #          "
R1_max = R1_max / 2 * safety_factor
R2_max = R2_max / 2 * safety_factor

def get_R_max():
    return [R1_max, R2_max]

print(L_Bot, L_Top)

# Getting max shear and max moment
moment_list_bot = R1_values*(1-l_values/L_Bot) * l_values 
if max(moment_list_bot) > abs(min(moment_list_bot)):
    max_moment_bot = max(moment_list_bot)
else: max_moment_bot = min(moment_list_bot)
moment_list_top = R2_values*(1-l2_values/L_Top) * l2_values
if max(moment_list_top) > abs(min(moment_list_top)):
    max_moment_top = max(moment_list_top)
else: max_moment_top = min(moment_list_top)



# Internal forces 

V_bot = R1_max  # Maximum shear force [N] (corresponds to maximum reaction force as well)
M_bot = max_moment_bot # Maximum Bending moment [Nm]

V_top = R2_max 
M_top = max_moment_top  #idk yet

print(f"Max Top Shear: {V_top}, Max Bot Shear: {V_bot}")
print(f"Max Top Moment: {M_top}, Max Bot Moment: {M_bot}")

# plt.subplot(211)
# plt.plot(l2_values,R2_values,"g",label="Top Beam")
# plt.plot(l_values,R1_values,"r",label="Bottom Beam")

# plt.xlabel('Position [m]')
# plt.ylabel("Force [N]")
# plt.title("Forces at each point in the beam")
# plt.legend()
# plt.grid()

# plt.subplot(212)
# plt.plot(l2_values,moment_list_top,"g",label="Top Beam")
# plt.plot(l_values,moment_list_bot,"r",label="Bottom Beam")


# plt.xlabel('Position [m]')
# plt.ylabel("Moment [Nm]")
# plt.title("Moments at each point in the beam")
# plt.legend()
# plt.grid()

# plt.show()

# Wheel dimensions (not fixed)
# BASED ON https://www.norelem.com/ca/en/Products/Product-overview/Material-handling-and-transport/95000-Material-handling-and-transport/Wheels-and-rollers/95059-Rollers-heavy-load.html
# other option 
wheel_width = 100e-3 # like 30 mm
wheel_radius = 85e-3 # 125 mm 
wheel_carry_force = 680*g # 750 kg "max load"
n_wheels_bot = V_bot / wheel_carry_force # approximately 4  4x1 lxw # 12 and 10 for 6 eur per wheel = 132 eur per 
n_wheels_top = V_top / wheel_carry_force # approximately 6 (6x1) lxw #
print(f"Number of wheels bottom: {n_wheels_bot}, Number of wheels top: {n_wheels_top}")

# Dimensions of the beams
b_bot = 120e-3    # width of beam [m] (2 wheels)
t_bot = 2e-3 # thickness [m]

b_top = 120e-3     # width of beam [m] (2 wheels)
t_top = 10e-3                # thickness [m]

# Material properties mild STEEL i think
rho = 7850 # Density in kg/m³
yield_stress = 344e6 # N/m² 



'''
Assuming square O-Beam
        b
|----------------|

==================  
||              ||
||            ->||<- t
||              ||
||              ||
||              ||
==================

'''

# Calculated values
Q_bot = b_bot**2/4 * t_bot + (b_bot-2*t_bot)*t_bot*(b_bot-t_bot)/2 # First moment of area [m³]
I_bot = (b_bot**4 -(b_bot-2*t_bot)**4)/12 # Second moment of area / moment of inertia [m^4]
A_bot = b_bot**2 - (b_bot-2*t_bot)**2 # Area [m²]
m_bot = A_bot * L_Bot * rho # Mass [kg]

Q_top = b_top**2/4 * t_top + (b_top-2*t_top)*t_top*(b_top-t_top)/2 # First moment of area [m³]
I_top = (b_top**4 -(b_top-2*t_top)**4)/12 # Second moment of area / moment of inertia [m^4]
A_top = b_top**2 - (b_top-2*t_top)**2 # Area [m²]
m_top = A_top * L_Top * rho # Mass [kg]

# Stresses
tau_bot = V_bot * Q_bot / I_bot / t_bot
sigma_bot = M_bot * b_bot / 2 / I_bot

tau_top = V_top * Q_top / I_top / t_top
sigma_top = M_top * b_top / 2 / I_top


print(f"Maximum shear stress along the bottom beam {tau_bot*10**-6} MPa")
print(f"Maximum tensile stress along the bottom beam {sigma_bot*10**-6} MPa")
print(f"Maximum shear stress along the top beam {tau_top*10**-6} MPa")
print(f"Maximum tensile stress along the top beam {sigma_top*10**-6} MPa")


print(f"Shear Stress Safety Margot: \n TOP: {abs(yield_stress/2/tau_top) - 1} BOT: {abs(yield_stress/2/tau_bot) - 1}")
print(f"Normal Stress Safety Magritte: \n TOP: {abs(yield_stress/sigma_top) - 1} BOT: {abs(yield_stress/sigma_bot) - 1}")
print(f"Mass of the bottom beam: {m_bot} kg") # We have two of these
print(f"Mass of the top beam: {m_top} kg")


# Interface! 
width_bearing = 26e-3 # https://nl.rs-online.com/web/p/roller-bearings/0312425
bearing_radius = 30e-3 #                      "
load_limit = 82e3 # Basic dynamic load rating, radial
pin_area = jojo.pi * bearing_radius**2 

width_top = b_top + 2*wheel_radius + 2*width_bearing
width_bot = b_bot + 2*wheel_radius + 2*width_bearing

force_in_bot_pin = R1_max
force_in_top_pin = R2_max

bot_shear_strorses = R1_max/pin_area
top_shear_strorses = R2_max/pin_area

print(f"{bot_shear_strorses},{top_shear_strorses}")

# around 12.2 MPa and 21.8 MPa , so definitely meets requirements!


'''
    width
{-----------}
rail
[---] 
     2x bearing
       [---]

|   ___
| _____||||
| O O  ||||
|===|| ||(
|   || ||( = = 
|===|| ||(
| O O  ||||
| _____||||
|   ___||||
|
'''

# Deflection!
#v = 1/EI \int \int M 

I_top = I_top
I_bot = I_bot

E = 207e9 # Young's modulus, 207 GPa for mild steel

