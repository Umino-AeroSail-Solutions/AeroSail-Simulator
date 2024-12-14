import numpy as np
import matplotlib.pyplot as plt

def compute_Ixx_Iyy_Ixy(w, h, t_top_bottom, t_sides,d, a):
    # See Cross_section.png, x is horizontal and centered pointing right. y is vertical and centered asumes corners are squares
    Ixx_corner = (1/12) * (d**4) + a*(((h-(d/2))/2)**2)
    Iyy_corner = (1/12) * (d**4) + a*(((w-(d/2))/2)**2)
    Ixx_corners = Ixx_corner*4
    Iyy_corners = Iyy_corner*4

    Ixx_horizontal_plates = 2*((t_top_bottom*w)*((h/2)**2))
    Iyy_horizontal_plates = t_top_bottom*(w**3)/6

    Ixx_vertical_plates = t_sides*(h**3)/6
    Iyy_vertical_plates = 2*((t_sides*h)*((w/2)**2))

    Ixx = Ixx_horizontal_plates + Ixx_vertical_plates + Ixx_corners
    Iyy = Iyy_horizontal_plates + Iyy_vertical_plates + Iyy_corners

    Ixy = 0

    return Ixx, Iyy, Ixy

def compute_area(w, h, t_top_bottom, t_sides,d, a):
    return 4*a + 2*t_top_bottom*w + 2*t_sides*h

def compute_bending_tension(Mx, My, Ixx, Iyy, Ixy, x, y):
    a = Mx*Iyy - My*Ixy
    b = My*Ixx - Mx*Ixy
    c = Ixx*Iyy - Ixy*Ixy
    tension = (a*y + b*x)/c
    return tension

def check_bending_ok(w, h, t_top_bottom, t_sides,d, a, Mx, My, tension_max):
    '''Returns true, SF if nothing breaks and false, SF if tension_max is reached'''
    Ixx, Iyy, Ixy = compute_Ixx_Iyy_Ixy(w, h, t_top_bottom, t_sides,d, a)
    test_coordinates = np.array([[w/2, h/2],[-w/2, h/2], [-w/2, -h/2], [w/2, -h/2]])
    tensions = np.array([0,0,0,0])
    for i in range(4):
        tensions[i] = compute_bending_tension(Mx, My, Ixx, Iyy, Ixy, test_coordinates[i,0], test_coordinates[i,1])
    SF = tension_max/tensions.max()
    if SF > 1:
        return True, SF
    else:
        return False, SF

def find_areas(w, h, t_top_bottom, t_sides,d, a, Mx, My, tension_max):
    failure, SF = check_bending_ok(w, h, t_top_bottom, t_sides,d, a, Mx, My, tension_max)
    while SF > 1.01 or SF < 1:
        failure, SF = check_bending_ok(w, h, t_top_bottom, t_sides, d, a, Mx, My, tension_max)
        a = a/SF
    return a

print(find_areas(2, 0.7, 0.001, 0.001, 0.08, 100, 100, 300, 100000))