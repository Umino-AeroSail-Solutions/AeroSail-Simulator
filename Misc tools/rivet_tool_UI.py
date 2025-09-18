from rivet_spacing_sizer import *


# Here comes a fastener list (Not only rivets in case you want to it also works with bolts):

Gesipa1433601 = Rivet(2000, 5.1/1000, 2800, 41/250, "Gesipa 1433601")
Gesipa1454074 = Rivet(3400, 6.5/1000, 4600, 60/200, "Gesipa 1454074")
Gesipa1454608 = Rivet(2000, 4.1/1000, 2500, 80/500, "Gesipa 1454608")
Gesipa1456041 = Rivet(1500, 4.9/1000, 2300, 55/250, "Gesipa 1456041")
Gesipa1433932 = Rivet(10000, 6.7/1000, 8000, 50/150, "Gesipa 1433932")
Gesipa1433811 = Rivet(5730, 6.6/1000, 3840, 52/100, "Gesipa 1433811")


rivet_list = [Gesipa1433601, Gesipa1454074, Gesipa1454608, Gesipa1456041, Gesipa1433932, Gesipa1433811]


# https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6063T6

alu6063_T6_shear_yield = 152e6
alu6063_T6_tension_yield = 214e6
alu6063_T6_eMod = 68.9e9
alu6063_T6_bearing_yield = 276e6

top_panel_thickness = 4/1000
side_panel_thickness = 4/1000
column_thickness = 10/1000

top_panel = Panel(top_panel_thickness, alu6063_T6_shear_yield, alu6063_T6_tension_yield, alu6063_T6_eMod, alu6063_T6_bearing_yield)
side_panel = Panel(side_panel_thickness, alu6063_T6_shear_yield, alu6063_T6_tension_yield, alu6063_T6_eMod, alu6063_T6_bearing_yield)
column = Panel(column_thickness, alu6063_T6_shear_yield, alu6063_T6_tension_yield, alu6063_T6_eMod, alu6063_T6_bearing_yield)

column_width = 80/1000
print("--------RIVET SIZER--------\n\n")
print("Select connection: ")
print(" 1. Column top panel")
print(" 2. Column side panel\n")

option = int(input("Enter number of choice: \n"))

print("\n")
stress = float(input("Now enter stress in Pa: "))

if option == 1:
    connection = Connection(stress, column, top_panel, column_width, SF=4.5)
    connection.find_The_Chosen_One(rivet_list, cost_multiplier=1, number_multiplier=0)
    print("\n")

if option == 2:
    connection = Connection(stress, column, side_panel, column_width, SF=4.5)
    connection.find_The_Chosen_One(rivet_list, cost_multiplier=1, number_multiplier=0)
    print("\n")