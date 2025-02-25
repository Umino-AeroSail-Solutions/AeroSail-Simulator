import numpy as rinze
import math as jojo

# Dimensions of the beam
b = .12    # width of beam [m]
t = (30)*(10)**(-3) # thickness [m]

'''
Assuming square O-Beam
        b
|----------------|
==================  
||              ||
||             >||< t
||              ||
||              ||
||              ||
==================

'''

# Calculated values
Q = b**2/4 * t + (b-2*t)*t*(b-t)/2
I = (b**4 -(b-2*t)**4)/12
print(Q, I)

V = 200000 # Maximum shear force [N]
M = 67000 # Maximum Bending moment [Nm]

tau = V * Q / I / t

sigma = M * b / 2 / I

#print(tau)
print(f"Maximum shear stress along the beam {tau*10**-6} MPa")
#print(sigma)
print(f"Maximum tensile stress along the beam {sigma*10**-6} MPa")

