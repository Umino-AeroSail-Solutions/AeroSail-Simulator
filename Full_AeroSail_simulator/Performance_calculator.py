from Sail import *
import Profile
import matplotlib.pyplot as plt
import numpy as np


Profile.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')
Sail = Sail_Class(os.path.join(
    os.path.dirname(__file__), ".", "Data", 'E473coordinates.txt')
    , 4.5, 0.4, 30, panels = 20)
Sail.load_interpolation(os.path.join(
        os.path.dirname(__file__), ".", "Data", 'interpolationCR4sail_XFLR5.npz'))

# # Code to plot average thrust as a function of boat speed and wind speed
interpolation = 'yes'
# interpolation = 'Data/interpolationCR4sail_XFLR5.npz'

# # Define the range of boat speeds and wind speeds
boat_speeds = np.array([12, 18, 23])  # Boat speeds
wind_speeds = np.array([15, 20, 25, 30, 35])  # Wind speeds
#
# # Create a matrix to store thrust values
thrust_values = np.zeros((len(boat_speeds), len(wind_speeds)))
#
ship_power = 15000000

# Just 1 test value
# Sail.plot_optimal_values_polar_with_thrust_and_strct(
#             np.arange(np.radians(0), np.radians(180), np.radians(1)),
#             5000, 15000, 25, [30, 150], shipspeed=23
#         , interpolation=interpolation)

# Compute thrust for each (boat_speed, wind_speed) pair
for j, wind_speed in enumerate(wind_speeds):
    for i, boat_speed in enumerate(boat_speeds):
        ship_thrust = ship_power/(boat_speed/1.944)
        thrust_values[i, j] = (Sail.plot_optimal_values_polar_with_thrust_and_strct(
            np.arange(np.radians(0), np.radians(180), np.radians(5)),
            5000, 15000, wind_speed, [30, 150], shipspeed=boat_speed
        , interpolation=interpolation)/ship_thrust)*100

# # Create a heatmap
#
thrust_values = np.array(thrust_values)

plt.figure(figsize=(8, 6))
sns.heatmap(thrust_values, xticklabels=np.round(wind_speeds, 1),
            yticklabels=np.round(boat_speeds, 1), cmap="coolwarm", annot=True, fmt=".2f", annot_kws={"size": 6})
plt.xlabel("Wind Speed (knots)")
plt.ylabel("Boat Speed (knots)")
plt.title("Average Thrust Heatmap")
# plt.show()