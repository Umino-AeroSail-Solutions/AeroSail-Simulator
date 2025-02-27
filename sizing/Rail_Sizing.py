import numpy as rinze
import math as jojo
from Force_In_Rail_Calculator_Andres_saves_the_day import L_Top, L_Bot

g = 9.81


# Internal forces #FIX THIS SO THEY PULL DIRECTLY FROM CALCULATIONS
V_bot = 200_000 # Maximum shear force [N] (corresponds to maximum reaction force as well)
M_bot = 67_000 # Maximum Bending moment [Nm]

V_top = 170_000
M_top = 50_000 #idk yet

# Wheel dimensions (not fixed)
wheel_width = 30e-3 # like 30 mm
wheel_radius = 125e-3 # 125 mm 
wheel_carry_force = 750*g # 750 kg "max load"
n_wheels_bot = V_bot / wheel_carry_force # approximately 28 right now, 4x7 grid should be fine
n_wheels_top = V_top / wheel_carry_force # approximately 24 right now, 3x8 grid should be fine (or 4x6)

# Dimensions of the beams
b_bot = .12    # width of beam [m] (4 wheels)
t_bot = (30)*(10)**(-3) # thickness [m]

b_top = .09     # width of beam [m] (3 wheels)
t_top = 20e-3   # thickness [m]

# Material properties mild STEEL i think
rho = 7850 # Density in kg/m³
yield_stress = 250e6 # N/m² 



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


print(f"Failure due to normal stress: {sigma_bot > yield_stress or sigma_top > yield_stress }")
print(f"Failure due to shear stress: {tau_bot > yield_stress / 2 or tau_top > yield_stress /2 }")
print(f"Mass of the bottom beam: {m_bot} kg")
print(f"Mass of the top beam: {m_top} kg")

