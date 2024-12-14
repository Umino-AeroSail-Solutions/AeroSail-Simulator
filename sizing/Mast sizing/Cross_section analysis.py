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

# SHEAR ANALYSIS-------------------------------------

def create_areas_and_thicknesses(w, h, d, a, subdivisions, top_thickness, side_thickness):
    up_pos = (h - d) / 2
    right_pos = (w - d) / 2
    areas = np.array([[right_pos, up_pos, a]])
    thicknesses = np.array([side_thickness])
    for i in range(subdivisions - 1):
        areas = np.append(areas, np.array([[right_pos, up_pos - ((i + 1) * up_pos * 2 / subdivisions), 0]]), axis=0)
        thicknesses = np.append(thicknesses, np.array([side_thickness]))
    areas = np.append(areas, np.array([[right_pos, -up_pos, a]]), axis=0)
    thicknesses = np.append(thicknesses, np.array([top_thickness]))
    for i in range(subdivisions - 1):
        areas = np.append(areas, np.array([[right_pos - ((i + 1) * right_pos * 2 / subdivisions), -up_pos, 0]]), axis=0)
        thicknesses = np.append(thicknesses, np.array([top_thickness]))
    areas = np.append(areas, np.array([[-right_pos, -up_pos, a]]), axis=0)
    thicknesses = np.append(thicknesses, np.array([side_thickness]))
    for i in range(subdivisions - 1):
        areas = np.append(areas, np.array([[-right_pos, -up_pos + ((i + 1) * up_pos * 2 / subdivisions), 0]]), axis=0)
        thicknesses = np.append(thicknesses, np.array([side_thickness]))
    areas = np.append(areas, np.array([[-right_pos, up_pos, a]]), axis=0)
    thicknesses = np.append(thicknesses, np.array([top_thickness]))
    for i in range(subdivisions - 1):
        areas = np.append(areas, np.array([[-right_pos + ((i + 1) * right_pos * 2 / subdivisions), up_pos, 0]]), axis=0)
        thicknesses = np.append(thicknesses, np.array([top_thickness]))
    return areas, thicknesses


def plot_areas(areas, thicknesses):
    # Extract x, y coordinates and sizes
    x = areas[:, 0]
    y = areas[:, 1]
    sizes = (areas[:, 2]) * 20 + 10  # Adjust scaling factor for better visualization

    # Plot points
    plt.scatter(x, y, s=sizes, c='blue', alpha=0.5)

    # Plot lines with varying thicknesses
    for i in range(len(x) - 1):
        plt.plot(x[i:i + 2], y[i:i + 2], linestyle='-', color='red', alpha=0.5, linewidth=thicknesses[i] * 1000)

    # Connect the last point to the first one
    plt.plot([x[-1], x[0]], [y[-1], y[0]], linestyle='-', color='red', alpha=0.5, linewidth=thicknesses[-1] * 1000)

    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.title('Areas Plot')
    plt.grid(True)
    plt.show()


# Example parameters
w, h, d, a = 1, 0.7, 0.04, 5
top_thickness = 0.001
sides_thickness = 0.002
subdivisions = 20

areas, thicknesses = create_areas_and_thicknesses(w, h, d, a, subdivisions, top_thickness, sides_thickness)
plot_areas(areas, thicknesses)
