tau_max = 100000000
sigma_max = 200000000
M = 2500  # kg
R = M * 9.81 / 2

# Generate arrays for all possible values of W, D, t
possible_dimensions = []

# Loop through all possible values of W, D, and t
for W in range(100, 475):  # Convert mm to meters
    for h in range(50, 300):
        for t in range(2, 12):
            possible_dimensions.append([W / 1000, h / 1000, t / 1000])  # Convert mm to meters

def check_shear(array):
    passed_array = []
    for dimension in array:
        W, h, t = dimension  # Unpacking for readability
        
        shear_stress = R * ((h ** 2 / 8) + (W * h / 4)) / (t * (h ** 2) * (W / 2 + (h / 6)))
        normal_stress = R * 2.35 / (8 * t * h * (W/2 + h/6))    

        if shear_stress < tau_max and normal_stress < sigma_max:
            passed_array.append(dimension)

    return passed_array

valid_dimensions = check_shear(possible_dimensions)

def smallest_area(array):
    if not array:
        return None, None  # Return None if no valid dimensions exist
    
    # Find the dimension with the smallest area
    min_dimension = min(array, key=lambda dim: 2 * dim[0] * dim[2] + 2 * dim[1] * dim[2])
    min_area = 2 * min_dimension[0] * min_dimension[2] + 2 * min_dimension[1] * min_dimension[2]

    return min_area, min_dimension

# Find the smallest area and the corresponding dimension
min_area, min_dim = smallest_area(valid_dimensions)

print(f"Smallest Area: {min_area} m²")
print(f"Optimal Dimension (W, h, t): {min_dim}")
