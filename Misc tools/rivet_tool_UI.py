from rivet_spacing_sizer import *


# Here comes a fastener list (Not only rivets in case you want to it also works with bolts):

# https://en.gesipa.de/products/blind-rivets/
Gesipa1433601 = Rivet(2000, 5.1/1000, 2800, 41/250, "Gesipa 1433601")
Gesipa1454074 = Rivet(3400, 6.5/1000, 4600, 60/200, "Gesipa 1454074")
Gesipa1454608 = Rivet(2000, 4.1/1000, 2500, 80/500, "Gesipa 1454608")
Gesipa1456041 = Rivet(1500, 4.9/1000, 2300, 55/250, "Gesipa 1456041")
Gesipa1433932 = Rivet(11000, 6.7/1000, 8000, 50/150, "Gesipa 1433932")
Gesipa1433811 = Rivet(5730, 6.6/1000, 3840, 52/100, "Gesipa 1433811")
Gesipa1433814 = Rivet(12455, 6.6/1000, 8200, 700/1000,"Gesipa 1433814")
BibusCN02_6420_0700 = Rivet(16000, 6.5/1000, 7500, 26/100, "Bibus CN02.6420.0700") # ONLY FOR THE THIN PANEL
BibusCN03_6419_0700 = Rivet(11700, 6.6/1000, 10400, 28/1000, "Bibus CN03.6419.700")
BN_21443 = Rivet(14500, 6.6/1000, 8330, 50/100, "BN 21443") # ONLY THIN
BN_84012 = Rivet(15600, 6.6/1000, 8500, 50/100, "BN 84012") # ONLY THIN
BN_84013 = Rivet(49400, 10.5/1000, 32300, 50/100, "BN 84013")



rivet_list_THICK = [Gesipa1433601, Gesipa1454074, Gesipa1454608, Gesipa1456041, Gesipa1433932, Gesipa1433811, Gesipa1433814, BibusCN03_6419_0700, BN_84013]
rivet_list_thin = [Gesipa1433601, Gesipa1454074, Gesipa1454608, Gesipa1456041, Gesipa1433932, Gesipa1433811, Gesipa1433814, BibusCN02_6420_0700, BibusCN03_6419_0700, BN_21443, BN_84012, BN_84013]


# https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6063T6

alu6063_T6_shear_yield = 152e6
alu6063_T6_tension_yield = 214e6
alu6063_T6_eMod = 68.9e9
alu6063_T6_bearing_yield = 276e6

top_panel_thickness = 4/1000
side_panel_thickness = 5/1000
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
    connection = Connection(stress, column, top_panel, column_width, SF=1.2)
    connection.find_The_Chosen_One(rivet_list_thin, cost_multiplier=1, number_multiplier=0)
    print("\n")

if option == 2:
    connection = Connection(stress, column, side_panel, column_width, SF=1.2)
    connection.find_The_Chosen_One(rivet_list_THICK, cost_multiplier=1, number_multiplier=0)
    print("\n")