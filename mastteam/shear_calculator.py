import numpy as np
import math
import matplotlib.pyplot as plt


#Define dimensions
w = 1440 #mm
h = 700 #mm
Ixx_8080 = Iyy_8080 = 1011255 #mm^4
t_tb = 0.0015 #mm
t_s = 0.0055 #mm
A_8080 = 1612.2 #mm^2

#Define loads
Vx = 100 #N
Vy = 100 #N
Mx = 2376234908000000 #Nm
My = 0 #Nm

#Define Idealizations
n_booms_side = 40
n_booms_top_bottom = 5

#Define Situation
deployed = True



def calcMOI(w, h, Ixx_8080, Iyy_8080, t_tb, t_s, A_8080):

    Ixx_sides = (1/12) * t_s * (h**3)
    Ixx_top_bottom = (1/12) * w * (t_tb**3)
    Iyy_sides = (1/12) * h * (t_s**3)
    Iyy_top_bottom = (1/12) * t_tb * (w**3)

    A_side = t_s * h
    A_top_bottom = t_tb * w
    d_side = w/2
    d_top_bottom = h/2

    Ixx_8080 = Ixx_8080 + A_8080 * (d_top_bottom ** 2)
    Iyy_8080 = Iyy_8080 + A_8080 * (d_side ** 2)
    Ixx_top_bottom = Ixx_top_bottom + A_top_bottom * (d_top_bottom ** 2)
    Iyy_sides = Iyy_sides + A_side * (d_side ** 2)

    Ixx = (Ixx_sides + Ixx_top_bottom) * 2 + (Ixx_8080 * 4)
    Iyy = (Iyy_sides + Iyy_top_bottom) * 2 + (Iyy_8080 * 4)

    return Ixx, Iyy

Ixx, Iyy = calcMOI(w, h, Ixx_8080, Iyy_8080, t_tb, t_s, A_8080)

print(f"Ixx = {Ixx}, Iyy = {Iyy}")


def calcNormStress(w, h, n_booms_side, n_booms_top_bottom):
    normStressArray = []
    coordsArray = []

    # 1. Right edge: top to bottom
    for i in range(n_booms_side):
        x = w/2
        y = h/2 - (h * i / (n_booms_side - 1))
        sigma_z = (Mx * Iyy * y + My * Ixx * x) / (Ixx * Iyy)
        normStressArray.append(sigma_z)
        coordsArray.append((x, y))

    # 2. Bottom edge: right to left
    for i in range(1, n_booms_top_bottom):
        x = w/2 - (w * i / (n_booms_top_bottom - 1))
        y = -h/2
        sigma_z = (Mx * Iyy * y + My * Ixx * x) / (Ixx * Iyy)
        normStressArray.append(sigma_z)
        coordsArray.append((x, y))

    # 3. Left edge: bottom to top
    for i in range(1, n_booms_side):
        x = -w/2
        y = -h/2 + (h * i / (n_booms_side - 1))
        sigma_z = (Mx * Iyy * y + My * Ixx * x) / (Ixx * Iyy)
        normStressArray.append(sigma_z)
        coordsArray.append((x, y))

    # 4. Top edge: left to right
    for i in range(1, n_booms_top_bottom - 1):
        x = -w/2 + (w * i / (n_booms_top_bottom - 1))
        y = h/2
        sigma_z = (Mx * Iyy * y + My * Ixx * x) / (Ixx * Iyy)
        normStressArray.append(sigma_z)
        coordsArray.append((x, y))

    return np.array(normStressArray), np.array(coordsArray)






def calcBoomArea(normStressArray):

    boomAreaArray = []
    b_side = h / (n_booms_side - 1)
    b_top_bottom = w / (n_booms_top_bottom - 1)
    totalBooms = (n_booms_side + n_booms_top_bottom) * 2 - 4
    
    
    for i in range(totalBooms ):
        if i == 0:
            boomArea = t_s * b_side / 6 * (2 + normStressArray[i + 1]/normStressArray[i]) + t_tb * b_top_bottom / 6 * (2 + normStressArray[normStressArray.size - 1]/normStressArray[i])
        elif i == totalBooms - 1:
            boomArea = t_tb * b_top_bottom / 6 * (2 + normStressArray[i - 1] / normStressArray[i]) + t_tb * b_top_bottom / 6 * (2 + normStressArray[0] / normStressArray[i])
        elif i == n_booms_side - 1:
            boomArea = t_tb * b_top_bottom / 6 * (2 + normStressArray[i + 1]/normStressArray[i]) + t_s * b_side / 6 * (2 + normStressArray[i - 1]/normStressArray[i])
        elif i == n_booms_side + n_booms_top_bottom - 2:
            boomArea = t_s * b_side / 6 * (2 + normStressArray[i + 1]/normStressArray[i]) + t_tb * b_top_bottom / 6 * (2 + normStressArray[normStressArray.size - 1]/normStressArray[i])
        elif i == 2 * n_booms_side + n_booms_top_bottom - 3:
            boomArea = t_tb * b_top_bottom / 6 * (2 + normStressArray[i + 1]/normStressArray[i]) + t_s * b_side / 6 * (2 + normStressArray[i - 1]/normStressArray[i])
        elif 0 < i < n_booms_side - 1 or n_booms_side + n_booms_top_bottom - 2 < i < 2 * n_booms_side + n_booms_top_bottom - 3:
            boomArea = t_s * b_side / 6 * (2 + normStressArray[i + 1]/normStressArray[i]) + t_s * b_side / 6 * (2 + normStressArray[i - 1]/normStressArray[i])
        elif n_booms_side - 1 < i < n_booms_side + n_booms_top_bottom - 2 or 2 * n_booms_side + n_booms_top_bottom - 3 < i < totalBooms:
            boomArea = t_tb * b_top_bottom / 6 * (2 + normStressArray[i + 1]/normStressArray[i]) + t_tb * b_top_bottom / 6 * (2 + normStressArray[i - 1]/normStressArray[i])

        #print(f"2 + i + 1/i + 2 + i-1/i = {4 + (normStressArray[i - 1]/normStressArray[i]) + (normStressArray[i + 1]/normStressArray[i]) }")

        boomAreaArray.append(boomArea)


    return np.array(boomAreaArray)

#def calcBoomArea(normStressArray):
    boomAreaArray = []
    b_side = h / (n_booms_side - 1)
    b_top_bottom = w / (n_booms_top_bottom - 1)
    totalBooms = (n_booms_side + n_booms_top_bottom) * 2 - 4

    for i in range(totalBooms):
        if i == 0:
            # Bottom-right corner (Right & Bottom skin)
            s1 = normStressArray[i]
            s2 = normStressArray[i + 1]
            s3 = normStressArray[-1]
            B = (t_s * b_side / 6 * (2 + s2 / s1)) + (t_tb * b_top_bottom / 6 * (2 + s3 / s1))

        elif i == totalBooms - 1:
            # Top-right corner (Right & Top skin)
            s1 = normStressArray[i]
            s2 = normStressArray[0]
            s3 = normStressArray[i - 1]
            B = (t_tb * b_top_bottom / 6 * (2 + s2 / s1)) + (t_s * b_side / 6 * (2 + s3 / s1))

        elif i == n_booms_side - 1:
            # Bottom-left corner (Bottom & Left skin)
            s1 = normStressArray[i]
            s2 = normStressArray[i + 1]
            s3 = normStressArray[i - 1]
            B = (t_tb * b_top_bottom / 6 * (2 + s2 / s1)) + (t_s * b_side / 6 * (2 + s3 / s1))

        elif i == n_booms_side + n_booms_top_bottom - 2:
            # Top-left corner (Left & Top skin)
            s1 = normStressArray[i]
            s2 = normStressArray[i + 1]
            s3 = normStressArray[i - 1]
            B = (t_s * b_side / 6 * (2 + s2 / s1)) + (t_tb * b_top_bottom / 6 * (2 + s3 / s1))

        elif 0 < i < n_booms_side - 1 or n_booms_side + n_booms_top_bottom - 2 < i < 2 * n_booms_side + n_booms_top_bottom - 3:
            # Right or Left edge (vertical skin)
            s1 = normStressArray[i]
            s2 = normStressArray[i + 1]
            s3 = normStressArray[i - 1]
            B = (t_s * b_side / 6 * (2 + s2 / s1)) + (t_s * b_side / 6 * (2 + s3 / s1))

        elif n_booms_side - 1 < i < n_booms_side + n_booms_top_bottom - 2 or 2 * n_booms_side + n_booms_top_bottom - 3 < i < totalBooms:
            # Bottom or Top edge (horizontal skin)
            s1 = normStressArray[i]
            s2 = normStressArray[i + 1]
            s3 = normStressArray[i - 1]
            B = (t_tb * b_top_bottom / 6 * (2 + s2 / s1)) + (t_tb * b_top_bottom / 6 * (2 + s3 / s1))

        boomAreaArray.append(B)

    return np.array(boomAreaArray)




def calcSF_x():

    return SF_arrayX



# Get stress and coordinates
norm_stress, boom_coords = calcNormStress(w, h, n_booms_side, n_booms_top_bottom)


# Plotting
x_vals = boom_coords[:, 0]
y_vals = boom_coords[:, 1]
stress_vals = norm_stress
boom_areas = calcBoomArea(norm_stress)
#print(boom_areas)
#print(stress_vals)

plt.figure(figsize=(10, 5))
sc = plt.scatter(x_vals, y_vals, c=stress_vals, cmap='viridis', s=10)
plt.colorbar(sc, label='Normal Stress (σ_z)')
plt.title('Stress Distribution Around Rectangular Cross Section')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.axis('equal')
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))
sc = plt.scatter(x_vals, y_vals, c=boom_areas, cmap='viridis', s=10)
plt.colorbar(sc, label='Boom Area')
plt.title('Boom Area Around Rectangular Cross Section')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.axis('equal')
plt.grid(True)
plt.show()