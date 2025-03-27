import numpy as rinze
import math as jojo
import scipy.integrate as suminlee
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

print(f"Bottom Length: {L_Bot}, Top Length: {L_Top}")

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
b_bot = 80e-3    # width of beam [m] 
a_bot = 160e-3 # height [m]
t_bot = 5e-3 # thickness [m]

b_top = 80e-3     # width of beam [m]
a_top = 190e-3    # height [m]
t_top = 7.5e-3     # thickness [m]

# Material properties mild STEEL i think
rho = 7850 # Density in kg/m³
yield_stress = 344e6 # N/m² 



'''
Assuming rectangular O-Beam
        b
|----------------|
                        ________
==================         ^
||              ||         |
||            ->||<- t     | 
||              ||         | 
||              ||         | a
||              ||         | 
||              ||         | 
||              ||         | 
||              ||         | 
||              ||         v 
==================      ________

'''

# Calculated values
Q_bot = (a_bot/2 * b_bot)*(a_bot/4) - (a_bot/2-t_bot)*(b_bot-2*t_bot)*(a_bot/4-t_bot/2) # First moment of area [m³]
I_bot = 1/12 * (a_bot**3 * b_bot - (a_bot-2*t_bot)**3 * (b_bot-2*t_bot)) # Second moment of area / moment of inertia [m^4]
A_bot = b_bot*a_bot - (b_bot-2*t_bot)*(a_bot-2*t_bot) # Area [m²]
m_bot = A_bot * L_Bot * rho # Mass [kg]

Q_top = (a_top/2 * b_top)*(a_top/4) - (a_top/2-t_top)*(b_top-2*t_top)*(a_top/4-t_top/2) # First moment of area [m³]
I_top = 1/12 * (a_top**3 * b_top - (a_top-2*t_top)**3 * (b_top-2*t_top)) # Second moment of area / moment of inertia [m^4]
A_top = b_top*a_top - (b_top-2*t_top)*(a_top-2*t_top) # Area [m²]
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
#deflection = 1/EI \int \int M 

E = 207e9 # Young's modulus, 207 GPa for mild steel


d_top_list = []
d_bot_list = []

def max_deflection(P,b,l,E,I): # Force and deflection
    return P*b*(l**2-b**2)**(3/2) / (9 * jojo.sqrt(3) * l * E * I )
    #return P*b**2 *(3*l-4*b) / 48 / E / I

for i in range(len(R1_values)):
    if l_values[i] > L_Bot/2:
        d_bot_list.append(max_deflection(R1_values[i],L_Bot - l_values[i],L_Bot,E,I_bot))
    else: d_bot_list.append(max_deflection(R1_values[i],l_values[i],L_Bot,E,I_bot))

    if l2_values[i] > L_Top / 2:
        d_top_list.append(max_deflection(R2_values[i],L_Top - l2_values[i],L_Top,E,I_top))
    else: d_top_list.append(max_deflection(R2_values[i],l2_values[i],L_Top,E,I_top))

plt.plot(212)
plt.plot(l_values,rinze.array(d_top_list) * 1e3,"g",label="Top Beam")
plt.plot(l_values,rinze.array(d_bot_list) * 1e3,"r",label="Bottom Beam")


plt.xlabel('Position [m]')
plt.ylabel("Deflection [mm]")
plt.axhline(10,color='black',linestyle='--')
plt.axhline(-10,color='black',linestyle='--')
plt.title("Maximum Deflection at each point in the beam")
plt.legend()
plt.grid()
plt.show()

# Buckling!
kc = 4 # Buckling coefficient, assume simply supported (conservative)
ks = 5.5 # Buckling coefficient for shear :D
v = .33 #poisson ratio



def check_buckling(skinsidetop,shearsidetop,skinsidebot,shearsidebot,sigma_top,sigma_bot,tau_top,tau_bot):
    # skin buckling!!
    crit_skin_buckling_top = jojo.pi**2 * kc * E / 12/(1-v**2) * (t_top/skinsidetop)**2
    crit_skin_buckling_bot = jojo.pi**2 * kc * E / 12/(1-v**2) * (t_bot/skinsidebot)**2
    print(f"Skin Buckling Safety Mango:\n TOP: {abs(crit_skin_buckling_top/sigma_top) - 1} BOT: {abs(crit_skin_buckling_bot/sigma_bot) - 1}")

    # shear buckling !!
    crit_shear_buckling_top = jojo.pi**2 * ks * E / 12/(1-v**2) * (t_top/shearsidetop)**2
    crit_shear_buckling_bot = jojo.pi**2 * ks * E / 12/(1-v**2) * (t_bot/shearsidebot)**2
    print(f"Shear Buckling Safety Mangosteen:\n TOP: {abs(crit_shear_buckling_top/tau_top) - 1} BOT: {abs(crit_shear_buckling_bot/tau_bot) - 1}")


check_buckling(b_top, a_top, b_bot, a_bot,sigma_top,sigma_bot,tau_top,tau_bot)

# Side force calculations
# WE GET THESE FROM ANDRÈS next week
V_bot_side = 22000 * safety_factor #TBD
V_top_side = 2.2e4 *safety_factor #TBD

M_bot_side = M_bot/2 * safety_factor #TBD
M_top_side = M_top/2 *safety_factor #TBD

Qy_bot = (b_bot/2 * a_bot)*(b_bot/4) - (b_bot/2-t_bot)*(a_bot-2*t_bot)*(b_bot/4-t_bot/2) # First moment of area [m³]
Iyy_bot = 1/12 * (b_bot**3 * a_bot - (b_bot-2*t_bot)**3 * (a_bot-2*t_bot)) # Second moment of area / moment of inertia [m^4]


Qy_top = (b_top/2 * a_top)*(b_top/4) - (b_top/2-t_top)*(a_top-2*t_top)*(b_top/4-t_top/2) # First moment of area [m³]
Iyy_top = 1/12 * (b_top**3 * a_top - (b_top-2*t_top)**3 * (a_top-2*t_top)) # Second moment of area / moment of inertia [m^4]


tau_bot_side = V_bot_side * Qy_bot / Iyy_bot / t_bot
sigma_bot_side = M_bot * b_bot / 2 / Iyy_bot

tau_top_side = V_top_side * Qy_top / Iyy_top / t_top
sigma_top_side = M_top * b_top / 2 / Iyy_top


print("\n\nSide forcesssssss")
check_buckling(a_top, b_top, a_bot, b_bot,sigma_top_side,sigma_bot_side,tau_top_side,tau_bot_side)

