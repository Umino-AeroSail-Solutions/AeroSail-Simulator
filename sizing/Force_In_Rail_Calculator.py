import math
from sympy.solvers import solve
from sympy import Symbol
import sympy
from UAS_Railing import l
x2 = Symbol('x2')
M = 2000000
H = 30

####################################################################
#Inputs: coordinate of P1, P2, P4, mass of mast M, height of mast H#
####################################################################

P1 = [0,1]
P2 = [2,1]
P3 = [8,2]
P4 = [8,0]


alpha = math.atan((P4[1] - P1[1])/(P4[0] - P1[0]))
beta = math.atan((P3[1] - P2[1])/(P3[0] - P2[0]))
x1 = l * math.cos(alpha)


D = abs(P2[0] - P1[0])

def y1(l):

    y1 = (P4[1] - P1[1])/(P4[0] - P1[0]) * l * math.cos(alpha) + P1[1]

    return y1

def solve_for_x2(l):
    x2_Value = 0
    solutions = solve(abs(P2[0] - P1[0]) - sympy.sqrt(pow((x2 - l),2) + pow((((P3[1] - P2[1])/(P3[0] - P2[0]) * x2 + P2[1]) - ((P4[1] - P1[1])/(P4[0] - P1[0]) * l * math.cos(alpha) + P1[1])), 2)), x2)
    print(solutions)
    for solution in solutions:
        if x2_Value < solution:
            x2_Value = solution
    return x2_Value

def y2(x2):

    y2 = (P3[1] - P2[1])/(P3[0] - P2[0]) * x2 + P2[1]

    return y2

phi = math.atan((y2(solve_for_x2(l)) - y1(l))/(solve_for_x2(l) - x1))

R2 = (M * 9.81 * H * math.sin(phi)) / D * math.cos(90 - phi + beta)
T = (R2 * (math.cos(beta) * math.tan(alpha) + math.sin(beta)) - M * 9.81 * math.tan(alpha)) / (math.tan(alpha) * math.sin(alpha) + math.cos(alpha))
R1 = (M * 9.81 + T * math.sin(alpha) - R2 * math.cos(beta)) / math.cos(alpha)

print(R1)
F_Bot = R1
print(R2)
F_Top = R2
print(T)