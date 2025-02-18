import math
from sympy.solvers import solve
from sympy import Symbol
import sympy
from sympy import symbols, Eq, solve, sqrt
x2 = Symbol('x2')
M = 2000000
H = 30

####################################################################
#Inputs: coordinate of P1, P2, P4, mass of mast M, height of mast H#
####################################################################

P1 = [0,1]
P2 = [2,1]
P3 = [40,11]
P4 = [6,0]

L_Bot = math.sqrt((P4[0]-P1[0])**2 + (P4[1]-P1[1])**2)

D = abs(P2[0] - P1[0])

def y1(l):

    y1 = (P4[1] - P1[1])/(P4[0] - P1[0]) * l * math.cos(alpha) + P1[1]

    return y1


def solve_for_x2(l):
    x2_value = 0


    return x2_value

"""def solve_for_x2(l):
    x2_Value = 0
    equation = P2[0] - P1[0] - (sympy.sqrt(pow((x2 - l), 2) + pow((((P3[1] - P2[1])/(P3[0] - P2[0]) * x2 + P2[1]) - ((P4[1] - P1[1])/(P4[0] - P1[0]) * l * math.cos(alpha) + P1[1])), 2)))
    print(f"Equation to solve at l={l}: {equation}")

    solutions = solve(P2[0] - P1[0] - (sympy.sqrt(pow((x2 - l),2) + pow((((P3[1] - P2[1])/(P3[0] - P2[0]) * x2 + P2[1]) - ((P4[1] - P1[1])/(P4[0] - P1[0]) * l * math.cos(alpha) + P1[1])), 2))), x2)
    print(solutions)
    for solution in solutions:
        if x2_Value < solution:
            x2_Value = solution
    print("Larger value:")
    print(x2_Value)
    return x2_Value"""

def solve_for_x2(l):
    """Solve for x2 ensuring real solutions."""
    eq = Eq(2 - l * sqrt((0.2083 * x2 - 1)**2 + 0.0270 * (0.3335 * x2 + 1)**2), 0)
    solutions = solve(eq, x2)
    real_solutions = [sol.evalf() for sol in solutions if sol.is_real]
    return max(real_solutions) if real_solutions else None

def y2(x2):

    y2 = (P3[1] - P2[1])/(P3[0] - P2[0]) * x2 + P2[1]

    return y2

alpha = math.atan((P4[1] - P1[1])/(P4[0] - P1[0]))
beta = math.atan((P3[1] - P2[1])/(P3[0] - P2[0]))

l = 0
dl = 0.1
F = []

print(L_Bot)

while l < L_Bot:

    x1 = l * math.cos(alpha)

    phi = abs(math.atan((y2(solve_for_x2(l)) - y1(l))/(solve_for_x2(l) - x1)))
    print(phi)

    R2 = (M * 9.81 * H * math.cos(phi)) / D * math.cos(90 - phi + beta)
    T = (R2 * (math.cos(beta) * math.tan(alpha) + math.sin(beta)) - M * 9.81 * math.tan(alpha)) / (math.tan(alpha) * math.sin(alpha) + math.cos(alpha))
    R1 = (M * 9.81 + T * math.sin(alpha) - R2 * math.cos(beta)) / math.cos(alpha)

    F.append(R1)
    l+=dl
    print("L value:")
    print(l)

    print(R1)
    # print(R2)
    # print(T)
#print(F)