import numpy as np
import matplotlib.pyplot as plt

import pygame

from fontTools.misc.symfont import green
from mpmath import degrees

# Define mass and height parameters
m = 2378  # Example mass in kg
h = 10   # Example height in meters
w = 1 # Example width of the sail typa ribs thingy idk man


# Initialize Pygame
#pygame.init()

# Screen dimensions
WIDTH, HEIGHT = 1000, 700
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Erection Visualizer")


# Font setup
pygame.font.init()
font = pygame.font.SysFont('Arial', 24)

# Colors
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
RED = (255, 0, 0)
BLUE = (0, 0, 255)
GREEN = (0, 255, 0)
YELLOW = (255, 255, 0)
BoruiBlue = (0,122,200)
MateRed = (200,50,50)


####################################################################
#Inputs: coordinate of P1, P2, P4, mass of mast M, height of mast H#
####################################################################

P1 = np.array([0, 1]) # Bottom rail in retracted
P2 = np.array([2.6, 1]) # Top rail in retracted
P3 = np.array([6, 2.65]) # Top rail in extended
P4 = np.array([6, .05]) # Bottom rail in extended

L_Bot = np.linalg.norm(P4-P1)
L_Top = np.linalg.norm(P3-P2)



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


def draw_arrow(screen, origin, angle, length, color):
    end_pos = (origin[0] + length * np.cos(angle), origin[1] - length * np.sin(angle))
    pygame.draw.line(screen, color, origin, end_pos, 3)
    pygame.draw.polygon(screen, color, [
        (end_pos[0], end_pos[1]),
        (end_pos[0] - 10 * np.cos(angle - np.pi / 6), end_pos[1] + 10 * np.sin(angle - np.pi / 6)),
        (end_pos[0] - 10 * np.cos(angle + np.pi / 6), end_pos[1] + 10 * np.sin(angle + np.pi / 6))
    ])

def get_reactions(P1, P2, P3, P4, l, m, h, draw=False, cogloc=h/2):
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
    D_vector = [C_vector[1],-1*C_vector[0]]

    # print(C_vector)
    # print(D_vector)
    phi = np.arctan2(C_vector[1], C_vector[0])
    # return phi # for testing

    alpha = -np.arctan2((P4[1] - P1[1]), (P4[0] - P1[0])) # Bottom rail angle from the horizontal CW
    beta = (np.arctan2((P3[1] - P2[1]), (P3[0] - P2[0]))) # Top rail angle from the horizontal ccw
    # print(np.degrees(alpha), np.degrees(beta))

    # print(np.cos(phi - beta))
    R2 = (m * 9.81 * (cogloc) * np.cos(phi)) / (d * np.cos(phi - beta))
    T = (R2 * (np.cos(beta) * np.tan(alpha) + np.sin(beta)) - m * 9.81 * np.tan(alpha)) / (
                np.tan(alpha) * np.sin(alpha) + np.cos(alpha))
    R1 = (m * 9.81 + T * np.sin(alpha) - R2 * np.cos(beta)) / np.cos(alpha)

    sum_x = R1 * np.sin(alpha) - R2 * np.sin(beta) + T * np.cos(alpha)
    sum_y = R1 * np.cos(alpha) + R2 * np.cos(beta) - T * np.sin(alpha) - m*9.81
    mom_a = (m * 9.81 * (cogloc) * np.cos(phi)) - (R2 * d * np.cos(phi - beta))

    if abs(sum_x)>0.1 or abs(sum_y)>0.1 or abs(mom_a)>0.1:
        print("ERROR ENCOUNTERED: not in equilibrium")
    #
    # print("Sum x is {}".format(sum_x))
    # print("Sum y is {}".format(sum_y))
    # print("Mom a is {}".format(mom_a))
    # print()
    scale = np.array([50,-50])
    offset = np.array([(WIDTH/2) - scale[0]*P3[0], 3*HEIGHT/4])
    if draw:
        pygame.draw.circle(screen, YELLOW, P1 * scale + offset, 4)
        pygame.draw.circle(screen, YELLOW, P2 * scale + offset, 4)
        pygame.draw.circle(screen, YELLOW, P3 * scale + offset, 4)
        pygame.draw.circle(screen, YELLOW, P4 * scale + offset, 4)
        pygame.draw.line(screen, WHITE, P1 * scale + offset, P4 * scale + offset, 1)
        pygame.draw.line(screen, WHITE, P2 * scale + offset, P3 * scale + offset, 1)

        pygame.draw.line(screen, WHITE, A*scale + offset, (A+((h/d)*(C_vector)))*scale + offset, 1)
        pygame.draw.circle(screen, YELLOW, A*scale + offset, 4)
        pygame.draw.circle(screen, YELLOW, B * scale + offset, 4)
        
        #Container Outlines
        pygame.draw.line(screen, BoruiBlue, [0,0] * scale + offset, [0,3]*scale+offset, 2)
        pygame.draw.line(screen, BoruiBlue, [0,0] * scale + offset, [12,0]*scale+offset, 2)
        pygame.draw.line(screen, BoruiBlue, [12,0] * scale + offset, [12,3]*scale+offset, 2)

        #Bottom right corner sail
        pygame.draw.line(screen,MateRed,(A+((h/d)*(C_vector)))*scale + offset, ((A+((h/d)*(C_vector)))+[(w/d)*(D_vector[0]), (w/d)*(D_vector[1])])*scale+offset ,1)
        #Top right corner sail
        pygame.draw.line(screen,MateRed,(A+((h/d)*(C_vector)))*scale + offset, ((A+((h/d)*(C_vector)))-[(.398/.2*w/d)*(D_vector[0]), (.398/.2*w/d)*(D_vector[1])])*scale+offset ,1)


        vector_scale = 1/200
        draw_arrow(screen, A*scale + offset, (np.pi/2-alpha), vector_scale*R1, RED)
        draw_arrow(screen, B * scale + offset, (np.pi / 2 + beta), vector_scale * R2, GREEN)
        draw_arrow(screen, A * scale + offset, (-(alpha)), vector_scale * T, BLUE)

        draw_arrow(screen, (A+((cogloc/d))*(C_vector))*scale + offset, (-np.pi / 2), m * 9.81 *vector_scale , YELLOW)
    return R1, R2, T, B



# Compute bottom rail length
L_Bot = np.linalg.norm(P4 - P1)

# Generate values of l from 0 to L_Bot
l_values = np.linspace(0, L_Bot, 1000)
l2_values = np.zeros(len(l_values))

# Compute reactions for each l
R1_values = []
R2_values = []
T_values = []

delta_t = 0

draw = True

cogloc = 4.12

for l in l_values:
    screen.fill(BLACK)
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
    R1, R2, T, B = get_reactions(P1, P2, P3, P4, l, m, h, draw=draw, cogloc=cogloc)
    R1_values.append(R1) 
    R2_values.append(R2)
    T_values.append(T)
    if draw:
        # Update the display
        pygame.display.flip()
        pygame.time.delay(delta_t)
    i = np.where(l_values == l)
    l2_values[i]=np.linalg.norm(P2-B)

L_Top = max(l2_values)
L_Bot = max(l_values)
R1_max = np.max(np.abs(np.array(R1_values)))
R2_max = np.max(np.abs(np.array(R2_values)))
# print("R1 max: ", np.max(np.abs(np.array(R1_values))))
# print("R2 max: ", np.max(np.abs(np.array(R2_values))))
# print("T max: ", np.max(np.abs(np.array(T_values))))

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

plt.figure(figsize=(8, 5))
plt.plot(l2_values, R2_values, label='R2', color='g')
plt.xlabel('$l2$ (Attachment Position on Top Rail)')
plt.ylabel('Reaction Forces')
plt.title('Reaction Forces as a Function of Attachment Position $l2$')
plt.legend()
plt.grid(True)
plt.show()


# Main loop
running = False
while running:
    screen.fill(BLACK)
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    mx, my = pygame.mouse.get_pos()
    l = (mx) * (L_Bot/WIDTH)
    R1, R2, T = get_reactions(P1, P2, P3, P4, l, m, h, draw=True, cogloc=cogloc)
    pygame.display.flip()

# Quit Pygame
pygame.quit()

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
#     phi = get_reactions(P1, P2, P3, P4, l, m, h, cogloc=cogloc)
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