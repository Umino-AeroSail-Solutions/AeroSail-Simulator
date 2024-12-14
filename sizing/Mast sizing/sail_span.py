import numpy as np
import math

#we want to find maximum load
max_load = 0
for deflec in flap_deflection:
    for angle in AoA:
        lift = 1/2 * Cl(angle, deflec) * rho * v^2 * c # lift per unit span
        drag = 1/2 * Cd(angle, deflec) * rho * v^2 * c # drag per unit span
        resultant = ( lift ** 2 + drag **2 ) ** (1/2)
        if resultant > max_load:
            max_load = resultant

max_moment = max_load * ( span / 2 + 2.896)

#maximum loads
max_compressive = ( max_moment * container_width / 2 ) / 2 + weight / 4
min_compressive = -1 * ( max_moment * container_width / 2 ) / 2 + weight / 4
max_shear = (max_compressive - min_compressive) + max_load

print(max_compressive, "maximum compressive force")
print(max_shear, "maximum shear in container walls")