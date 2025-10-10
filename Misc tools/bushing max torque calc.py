import numpy as np
import math

def get_momentSF(max_ax_load, d, D, ax_load_applied, moment):
    R = D/2
    r = d/2
    area = math.pi*(R**2 - r**2)
    max_stress = max_ax_load/area
    applied_ax_stress = ax_load_applied/area

    I = (math.pi/64) * (D**4-d**4)
    moment_stress = R*moment/I
    stress = abs(applied_ax_stress) + 2*abs(moment_stress)
    SF = max_stress/stress
    return SF

# print("\nPCMW 629002 E")
# print(get_momentSF(265000, 62/1000, 90/1000, 3000, 5400), "\n") # BAD
#
# print("\nPCMW 386201.5 E")
# print(get_momentSF(150000, 38/1000, 62/1000, 3000, 5400), "\n") # BAD

print("\nBD1B 351890 A")
print(get_momentSF(1330000, 1380/1000, 1535/1000, 3000, 5400), "\n") # GOOD not found

print("\nSKF 51252 M")
print(get_momentSF(488000, 260/1000, 360/1000, 3000, 5400), "\n") # GOOD 793 euro
# https://www.tuli-shop.com/zen-bearing-51252-260x360x79mm?srsltid=AfmBOopNDZLVaOtDnJQ6GkkQ__7xqPt-G-sQoiJf1Ikmw6DZjtC_bkwc

# print("\nSKF 51148 M")
# print(get_momentSF(134000, 240/1000, 300/1000, 3000, 5400), "\n") # BAD

print("\nSKF 51340 M")
print(get_momentSF(624000, 200/1000, 340/1000, 3000, 5400), "\n") # GOOD 1390 euro
# https://www.tuli-shop.com/zen-bearing-51340-m-200x340x110mm?srsltid=AfmBOopjlC2laUOBKm194vDD4-5u-xzXrrkU8xlGLGuZQBle4vt8N0Ka

print("\nNSK 51148X")
print(get_momentSF(229000, 240/1000, 300/1000, 3000, 5400), "\n") # PERFECT 400 euro
# https://www.oss.nsk.com/products/bearings/ball-bearings/thrust-ball-bearings/single-direction-thrust-ball-bearings/51148x-tb-sd.html