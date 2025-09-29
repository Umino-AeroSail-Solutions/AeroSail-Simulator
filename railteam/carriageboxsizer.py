import numpy as np
from matplotlib import pyplot as plt
import math

#         a
# _________________
# |               |
# |               |
# |               |   b
# |               |
# |               |

# t_a is the thickness of block a and so on

alu6061_T6_shear_yield = 207e6
alu6061_T6_tension_yield = 276e6
alu6061_T6_eMod = 68.9e9
alu6061_T6_bearing_yield = 607e6


def get_Ixx(a, b, t_a, t_b, Print=False):
    c_y = (a * t_a * (b - (t_a / 2)) + 2*b * t_b * (b / 2)) / (a * t_a + 2*b * t_b)
    Ixx_top = ( a * (t_a ** 3) / 12 ) + (a*t_a*(b - (t_a / 2) - c_y)**2)
    Ixx_side = (t_b * (b ** 3) / 12) + (b*t_b * ((b / 2) -c_y)**2)
    Ixx = Ixx_top + 2*Ixx_side
    if Print:
        print("Ixx:", Ixx)
        print("C_y", c_y)
    return Ixx, c_y

def get_qb(a, b, t_a, t_b, y, shearup):
    Ixx, c_y = get_Ixx(a, b, t_a, t_b)
    qb = (shearup/Ixx) * (y*t_b) * ((y/2)-c_y)
    return qb


def test_beam(a, b, t_a, t_b, force_up, moment_up ,shear_yield, tension_yield):
    Ixx, c_y = get_Ixx(a, b, t_a, t_b)
    bottom_tension = c_y * moment_up / Ixx
    top_compression = (b-c_y) * moment_up / Ixx
    max_shear = abs(get_qb(a, b, t_a, t_b, c_y, force_up)/t_b)
    tension_SF = tension_yield/max(bottom_tension, top_compression)
    shear_SF = shear_yield/max_shear
    return shear_SF, tension_SF

def optimize_carriage_box(a, b, force_up, width,
                          shear_yield, tension_yield,
                          aim_SF=1.5,
                          t_a=20/1000, t_b=2/1000,
                          tol=1e-2, max_iter=100):

    moment_up = force_up * width / 2
    shear_SF, tension_SF = test_beam(a, b, t_a, t_b,
                                     force_up, moment_up,
                                     shear_yield, tension_yield)

    # Store history for plotting
    history = {
        "t_a": [t_a],
        "t_b": [t_b],
        "shear_SF": [shear_SF],
        "tension_SF": [tension_SF]
    }

    for _ in range(max_iter):
        if (abs(shear_SF - aim_SF) < tol) and (abs(tension_SF - aim_SF) < tol):
            break

        # Update rules (proportional control)
        t_a = t_a / (tension_SF / aim_SF)
        t_b = t_b / (shear_SF / aim_SF)

        print("t_a: ", t_a)
        print("t_b: ", t_b)
        print("shear_SF: ", shear_SF)
        print("tension_SF: ", tension_SF)
        print("\n")

        shear_SF, tension_SF = test_beam(a, b, t_a, t_b,
                                         force_up, moment_up,
                                         shear_yield, tension_yield)

        # Record iteration
        history["t_a"].append(t_a)
        history["t_b"].append(t_b)
        history["shear_SF"].append(shear_SF)
        history["tension_SF"].append(tension_SF)

    # Plot convergence
    fig, axs = plt.subplots(1, 2, figsize=(10,4))

    axs[0].plot(history["t_a"], label="t_a [m]")
    axs[0].plot(history["t_b"], label="t_b [m]")
    axs[0].set_title("Thickness convergence")
    axs[0].set_xlabel("Iteration")
    axs[0].set_ylabel("Thickness [m]")
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(history["shear_SF"], label="Shear SF")
    axs[1].plot(history["tension_SF"], label="Tension SF")
    axs[1].axhline(aim_SF, color="k", linestyle="--", label="Target SF")
    axs[1].set_title("Safety factor convergence")
    axs[1].set_xlabel("Iteration")
    axs[1].set_ylabel("Safety factor")
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.show()

    return t_a, t_b, shear_SF, tension_SF

best_top_carriage_top_rollers = optimize_carriage_box(
    a=530/1000, b=250/1000,
    force_up=30000, width=125/1000,
    shear_yield=alu6061_T6_shear_yield,
    tension_yield=alu6061_T6_tension_yield,
    aim_SF=2
)
print("Optimized design, top carriage, top rollers:", best_top_carriage_top_rollers)



moment_up = 30000*125/1000
top_carriage_top_rollers = test_beam(a=530/1000, b=250/1000, t_a=0.0000000000001, t_b=3/1000,force_up=30000, moment_up=moment_up, shear_yield=alu6061_T6_shear_yield, tension_yield=alu6061_T6_tension_yield)
print("Real design sf:", top_carriage_top_rollers)
print("max shear: ", alu6061_T6_shear_yield/top_carriage_top_rollers[0])

max_shear_panel = alu6061_T6_shear_yield/top_carriage_top_rollers[0]
k=6
buckling_shear = k*(math.pi**2)*alu6061_T6_eMod*((3/1000)**2)/(12*(1-0.33**2)*(250/1000)**2)
print("Buckling shear: ", buckling_shear)
print("Buckling safety factor: ", buckling_shear/max_shear_panel)

# best_top_carriage_bot_rollers = optimize_carriage_box(
#     a=140/1000, b=75/1000,
#     force_up=6000, width=125/1000,
#     shear_yield=alu6061_T6_shear_yield,
#     tension_yield=alu6061_T6_tension_yield,
#     aim_SF=2,
#     t_a = 0.000001,
#     t_b = 0.002
# )
# print("Optimized design, top carriage, bot rollers:", best_top_carriage_bot_rollers)

# TEST
# a = 1
# b = 1
# t_a = 20/1000
# t_b = 1/1000
#
# Ixx, c_y = get_Ixx(a, b, t_a, t_b, Print=True)
#
# bs =  np.arange(0, b, 1/10000)
# qbs = np.zeros(len(bs))
# index = 0
# for y in bs:
#     qbs[index] = get_qb(a, b, t_a, t_b, y, 100)
#     index += 1
#
# plt.plot(bs, qbs)
# plt.show()