import numpy as np
import matplotlib.pyplot as plt


####################################################################
#Inputs: coordinate of P1, P2, P4, mass of mast M, height of mast H#
####################################################################

P1 = np.array([0,1]) # Bottom rail in retracted
P2 = np.array([2,1]) # Top rail in retracted
P3 = np.array([8,2]) # Top rail in extended
P4 = np.array([8,0]) # Bottom rail in extended


def circle_line_intersection(center, diameter, P1, P2):
    center = np.array(center)
    P1, P2 = np.array(P1), np.array(P2)

    radius = diameter / 2  # Convert diameter to radius
    x1, y1 = P1
    x2, y2 = P2
    cx, cy = center

    # Parametric line equation: (x, y) = (x1, y1) + t * ((x2 - x1), (y2 - y1))
    dx, dy = x2 - x1, y2 - y1

    # Quadratic coefficients for intersection
    a = dx ** 2 + dy ** 2
    b = 2 * (dx * (x1 - cx) + dy * (y1 - cy))
    c = (x1 - cx) ** 2 + (y1 - cy) ** 2 - radius ** 2

    # Solve quadratic equation: a*t^2 + b*t + c = 0
    discriminant = b ** 2 - 4 * a * c

    if discriminant < 0:
        print("No intersection")
        return None

    t1 = (-b + np.sqrt(discriminant)) / (2 * a)
    t2 = (-b - np.sqrt(discriminant)) / (2 * a)

    # Compute intersection points
    intersection1 = np.array([x1 + t1 * dx, y1 + t1 * dy])
    intersection2 = np.array([x1 + t2 * dx, y1 + t2 * dy])

    return intersection1, intersection2


def get_reactions(P1, P2, P3, P4, l, m, h):
    # Some useful dimensions
    L_Bot = np.linalg.norm(P4-P1)
    L_Top = np.linalg.norm(P3-P2)
    d = np.linalg.norm(P2-P1)

    if l < 0 or l > L_Bot:
        raise ValueError("l is supposed to be between 0 and {}, to remain within the bottom rail".format(L_Bot))

    # Now we define point A, the point at which the attachment happens in the lower rail

    A = P1 + (l * ((P4-P1)/L_Bot))

    # Now you find the intersections between a circle centered in point A with radius d and the line passing through P3 and P2

    B1, B2 = circle_line_intersection(A, d*2, P3, P2)

    if B1[0] > B2[0]:
        B = B1
    else:
        B = B2

    C_vector = B-A

    phi = np.atan2(C_vector[1], C_vector[0])
    # return phi # for testing

    alpha = np.atan2((P4[1] - P1[1]), (P4[0] - P1[0]))
    beta = np.atan2((P3[1] - P2[1]), (P3[0] - P2[0]))

    R2 = (m * 9.81 * h * np.sin(phi)) / d * np.cos(90 - phi + beta)
    T = (R2 * (np.cos(beta) * np.tan(alpha) + np.sin(beta)) - m * 9.81 * np.tan(alpha)) / (
                np.tan(alpha) * np.sin(alpha) + np.cos(alpha))
    R1 = (m * 9.81 + T * np.sin(alpha) - R2 * np.cos(beta)) / np.cos(alpha)

    return R1, R2, T


# Compute bottom rail length
L_Bot = np.linalg.norm(P4 - P1)

# Generate values of l from 0 to L_Bot
l_values = np.linspace(0, L_Bot, 100)

# Define mass and height parameters
m = 2000  # Example mass in kg
h = 30   # Example height in meters

# Compute reactions for each l
R1_values = []
R2_values = []
T_values = []

for l in l_values:
    R1, R2, T = get_reactions(P1, P2, P3, P4, l, m, h)
    R1_values.append(R1)
    R2_values.append(R2)
    T_values.append(T)

# Plot reaction forces as a function of l
plt.figure(figsize=(8, 5))
plt.plot(l_values, R1_values, label='R1', color='r')
plt.plot(l_values, R2_values, label='R2', color='g')
plt.plot(l_values, T_values, label='T', color='b')
plt.xlabel('$l$ (Attachment Position on Bottom Rail)')
plt.ylabel('Reaction Forces')
plt.title('Reaction Forces as a Function of Attachment Position $l$')
plt.legend()
plt.grid(True)
plt.show()


# Phi testing shenanigans:
# # Compute bottom rail length
# L_Bot = np.linalg.norm(P4 - P1)
#
# # Generate values of l from 0 to L_Bot
# l_values = np.linspace(0, L_Bot, 100)
# phi_values = []
# m,h = 0,0
# for l in l_values:
#     print(l)
#     phi = get_reactions(P1, P2, P3, P4, l, m, h)
#     print(phi)
#     print()
#     phi_values.append(phi)
#
# # Convert phi to degrees for better readability
# phi_values = np.degrees(phi_values)
#
# # Plot phi as a function of l
# plt.figure(figsize=(8, 5))
# plt.plot(l_values, phi_values, label=r'$\phi$ vs. $l$', color='b')
# plt.xlabel('$l$ (Attachment Position on Bottom Rail)')
# plt.ylabel(r'$\phi$ (degrees)')
# plt.title('Angle $\phi$ as a Function of Attachment Position $l$')
# plt.legend()
# plt.grid(True)
# plt.show()