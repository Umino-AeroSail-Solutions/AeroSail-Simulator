import numpy as np
import matplotlib.pyplot as plt



#Checked 15/12/2024
def compute_Ixx_Iyy_Ixy(w, h, t_top_bottom, t_sides,d, a):
    # See Cross_section.png, x is horizontal and centered pointing right. y is vertical and centered asumes corners are squares
    Ixx_corner = (1/12) * (d**4) + a*(((h/2 - d/2))**2)
    Iyy_corner = (1/12) * (d**4) + a*(((w/2 - d/2))**2)
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

#Checked 15/12/2024
def compute_area(w, h, t_top_bottom, t_sides,d, a):
    return 4*a + 2*t_top_bottom*w + 2*t_sides*h

#Checked 15/12/2024
def compute_bending_tension(Mx, My, Ixx, Iyy, Ixy, x, y):
    a = Mx*Iyy - My*Ixy
    b = My*Ixx - Mx*Ixy
    c = Ixx*Iyy - Ixy*Ixy
    tension = (a*y + b*x)/c  #why is there a minus sign here? apparently its going to break
    return tension

#Checked 15/12/2024
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

#Checked 15/12/2024
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
    sizes = (areas[:, 2]) * 1000000 + 10  # Adjust scaling factor for better visualization

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

def get_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

def update_areas_bending(areas, thicknesses, Mx, My, Ixx, Iyy, Ixy):
    i=0
    for area in areas:
        term_2 = 2+ compute_bending_tension(Mx, My, Ixx, Iyy, Ixy, areas[i-1, 0], areas[i-1, 1])/compute_bending_tension(Mx, My, Ixx, Iyy, Ixy, areas[i, 0], areas[i, 1])
        extra_area_1 = term_2*(thicknesses[i-1]*get_distance(areas[i-1, 0], areas[i-1, 1], areas[i,0], areas[i,1]))/6
        if i+1 < np.size(areas, axis=0):
            term_2 = 2+ compute_bending_tension(Mx, My, Ixx, Iyy, Ixy, areas[i+1, 0], areas[i+1, 1])/compute_bending_tension(Mx, My, Ixx, Iyy, Ixy, areas[i, 0], areas[i, 1])
            extra_area_2 = term_2 * (thicknesses[i] * get_distance(areas[i + 1, 0], areas[i + 1, 1], areas[i, 0], areas[i, 1])) / 6
            areas[i, 2] += extra_area_1+extra_area_2
        else:
            term_2 = 2 + compute_bending_tension(Mx, My, Ixx, Iyy, Ixy, areas[0, 0], areas[0, 1]) / compute_bending_tension(Mx, My, Ixx, Iyy, Ixy, areas[i, 0], areas[i, 1])
            extra_area_2 = term_2 * (thicknesses[i] * get_distance(areas[0, 0], areas[0, 1], areas[i, 0], areas[i, 1])) / 6
            areas[i, 2] += extra_area_1 + extra_area_2
        i+=1
def parallel_axis_theorem(Ic, A, d):
    return Ic + A * d**2
def compute_Ixx_Iyy_with_parallel_axis(areas):
    Ixx = 0
    Iyy = 0
    centroid_x = np.mean(areas[:, 0])
    centroid_y = np.mean(areas[:, 1])

    for area in areas:
        x, y, A = area
        d_x = x - centroid_x
        d_y = y - centroid_y
        # Assuming point areas
        Ic = 0  # Moment of inertia for a square about its centroid
        Ixx += parallel_axis_theorem(Ic, A, d_y)
        Iyy += parallel_axis_theorem(Ic, A, d_x)

    return Ixx, Iyy

def compute_shear_flows(areas, Vx, Vy, subdivisions):
    Ixy = 0
    Ixx, Iyy = compute_Ixx_Iyy_with_parallel_axis(areas)
    a = -1* (Vy*Iyy)/(Ixx*Iyy)
    b = -1* (Vx*Ixx)/(Ixx*Iyy)
    qb = np.empty_like(thicknesses)
    for i in range(np.size(areas, axis=0)):
        sum_br_y = np.sum(np.multiply(areas[:(i+1), 2], areas[:(i+1), 1]))
        sum_br_x = np.sum(np.multiply(areas[:(i+1), 2], areas[:(i+1), 0]))
        qb[i] = a*sum_br_y + b*sum_br_x
    totaltorque=0
    moment_arm=(h-d)/2
    length = (w-d)/subdivisions
    index=0
    for i in range(4):
        if moment_arm==(h-d)/2:
            moment_arm=(w-d)/2
            length = (h - d)/subdivisions
        else:
            moment_arm=(h-d)/2
            length = (w - d) / subdivisions
        for j in range(subdivisions):
            totaltorque += qb[index]*moment_arm*length
            index += 1
    qs0 = np.zeros_like(qb)
    qs0[:] = -1*(totaltorque/(2*(h-d)*(w-d)))
    q = np.add(qb, qs0)
    return q, q.max()

def check_shear_ok(shear_flow, t_top_bottom, t_sides, max_allow_shear):
    shear = np.zeros_like(thicknesses)
    SF_shear = 0
    for i in range(np.size(shear)):
        if (i >= 0 and i <= int(np.size(shear)/4)) or (i >= int(np.size(shear)/2) and i <= int(3*np.size(shear)/4)):
            thickness = t_sides
        else:
            thickness = t_top_bottom

        shear[i] = np.abs(np.divide(shear_flow[i], thickness))
    SF_shear = max_allow_shear/shear.max()
    print(SF_shear)
    if SF_shear > 1:
        return True, SF_shear
    else:
        return False, SF_shear



def plot_shear_flows(areas, shear_flows):
    # Extract x, y coordinates and sizes
    x = areas[:, 0]
    y = areas[:, 1]
    sizes = (areas[:, 2]) * 1000000 + 10  # Adjust scaling factor for better visualization

    # Plot points
    plt.scatter(x, y, s=sizes, c='blue', alpha=0.5)

    # Plot lines with varying thicknesses and shear flow values as annotations
    for i in range(len(x) - 1):
        plt.plot(x[i:i + 2], y[i:i + 2], linestyle='-', color='red', alpha=0.5, linewidth=thicknesses[i] * 1000)
        plt.text((x[i] + x[i + 1]) / 2, (y[i] + y[i + 1]) / 2, f'{shear_flows[i]:.2f}', color='green')

    # Connect the last point to the first one
    plt.plot([x[-1], x[0]], [y[-1], y[0]], linestyle='-', color='red', alpha=0.5, linewidth=thicknesses[-1] * 1000)
    plt.text((x[-1] + x[0]) / 2, (y[-1] + y[0]) / 2, f'{shear_flows[-1]:.2f}', color='green')

    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.title('Shear Flows Plot')
    plt.grid(True)
    plt.show()
def plot_shear_flows_with_arrows(areas, shear_flows, Vx, Vy):
    # Extract x, y coordinates and sizes
    x = areas[:, 0]
    y = areas[:, 1]
    sizes = (areas[:, 2]) * 1000000 + 10  # Adjust scaling factor for better visualization

    # Plot points
    plt.scatter(x, y, s=sizes, c='blue', alpha=0.5)

    # Compute arrow directions and lengths
    dx = np.diff(x, append=x[0])
    dy = np.diff(y, append=y[0])

    # Adjust start positions and directions based on the sign of shear flows
    start_x = np.zeros_like(x)
    start_y = np.zeros_like(y)
    for i in range(len(shear_flows)):
        if shear_flows[i] < 0:
            if i+1 < np.size(x):
                start_x[i], start_y[i] = x[i + 1], y[i + 1]
            else:
                start_x[i], start_y[i] = x[0], y[0]
            dx[i], dy[i] = -dx[i], -dy[i]
        else:
            start_x[i], start_y[i] = x[i], y[i]

    # Normalize the lengths of the arrows based on shear flows
    lengths = np.sqrt(dx ** 2 + dy ** 2)
    arrow_lengths = shear_flows / np.max(np.abs(shear_flows)) * lengths

    # Plot arrows for shear flows
    plt.quiver(start_x, start_y, dx, dy, shear_flows, angles='xy', scale_units='xy', scale=1, color='green')

    # Plot Vx and Vy vectors
    centroid_x = np.mean(x)
    centroid_y = np.mean(y)
    plt.quiver(centroid_x, centroid_y, Vx / 700000, Vy / 700000, color='purple', scale=1, scale_units='xy', angles='xy')

    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.title('Shear Flows Plot with Arrows and Vx, Vy Vector')
    plt.grid(True)
    plt.show()


def plot_shear_flows_with_arrows_and_tension(areas, shear_flows, Vx, Vy, Mx, My, Ixx, Iyy, Ixy):
    # Extract x, y coordinates and sizes
    x = areas[:, 0]
    y = areas[:, 1]
    sizes = (areas[:, 2]) * 1000000 + 2  # Adjust scaling factor for better visualization

    # Calculate tension/compression at each point
    tensions = np.array([compute_bending_tension(Mx, My, Ixx, Iyy, Ixy, xi, yi) for xi, yi in zip(x, y)])

    # Define colors based on tension and compression
    colors = ['red' if t > 0 else 'blue' for t in tensions]

    # Plot points with colors based on tension/compression
    plt.scatter(x, y, s=sizes, c=colors, alpha=0.5)

    # Compute arrow directions and lengths
    dx = np.diff(x, append=x[0])
    dy = np.diff(y, append=y[0])

    # Adjust start positions and directions based on the sign of shear flows
    start_x = np.zeros_like(x)
    start_y = np.zeros_like(y)
    for i in range(len(shear_flows)):
        if shear_flows[i] < 0:
            if i + 1 < np.size(x):
                start_x[i], start_y[i] = x[i + 1], y[i + 1]
            else:
                start_x[i], start_y[i] = x[0], y[0]
            dx[i], dy[i] = -dx[i], -dy[i]
        else:
            start_x[i], start_y[i] = x[i], y[i]

    # Normalize the lengths of the arrows based on shear flows
    lengths = np.sqrt(dx ** 2 + dy ** 2)
    arrow_lengths = shear_flows / np.max(np.abs(shear_flows)) * lengths

    # Plot arrows for shear flows
    plt.quiver(start_x, start_y, dx, dy, shear_flows, angles='xy', scale_units='xy', scale=1, color='green')

    # Plot Vx and Vy vectors
    centroid_x = np.mean(x)
    centroid_y = np.mean(y)
    plt.quiver(centroid_x, centroid_y, Vx, Vy, color='purple', scale=1000000, scale_units='xy', angles='xy')
    plt.quiver(centroid_x, centroid_y, Mx, My, color='green', scale=1000, scale_units='xy', angles='xy')
    plt.quiver(centroid_x, centroid_y, Mx, My, color='green', scale=1300, scale_units='xy', angles='xy')

    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.title('Shear Flows Plot with Arrows, Vx, Vy Vector, and Tension/Compression')
    plt.grid(True)
    plt.show()


# Example usage with your current setup
w, h, d = 1, 0.7, 0.04
top_thickness = 0.002
sides_thickness = 0.001
subdivisions = 14

a = find_areas(w, h, top_thickness, sides_thickness, d, 100, 100, 300, 200000)
Ixx, Iyy, Ixy = compute_Ixx_Iyy_Ixy(w, h, top_thickness, sides_thickness, d, a)
areas, thicknesses = create_areas_and_thicknesses(w, h, d, a, subdivisions, top_thickness, sides_thickness)
update_areas_bending(areas, thicknesses, 100, 100, Ixx, Iyy, Ixy)
plot_areas(areas, thicknesses)
qb, qbmax = compute_shear_flows(areas, 100000, -100000, subdivisions)
#print(qb)
#283000000 is the aluminium 2024 shear strength
check_shear_ok(qb, top_thickness, sides_thickness, 283000000)
plot_shear_flows_with_arrows_and_tension(areas, qb, 100000, -100000, 100, 100, Ixx, Iyy, Ixy)