from rivet_spacing_sizer import *


A2yieldshear = 0.75*450e6

# Here comes a fastener list (Not only rivets in case you want to it also works with bolts):


# Blind rivets https://en.gesipa.de/products/blind-rivets/
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
Honsel_10792078225 = Rivet(15700, 7.8/1000, 9100, 50/100, "Honsel 10792078225") # https://www.honsel.de/en/products/productdetails/structural-blind-rivet-fero-bulb-steelsteel-dome-head/
Honsel_10792064205 = Rivet(16500, 6.4/1000, 7800, 50/100, "Honsel 10792064205") # Only thin (techically)
Estrubolt6_4 = Rivet(22600, 6.4/1000, 14500, 50/100, "Estrubolt 6.4")
Estrubolt8_0 = Rivet(35800, 8.0/1000, 23100, 50/100, "Estrubolt 8.0")
Estrubolt9_6 = Rivet(49400, 9.6/1000, 32300, 50/100, "Estrubolt 9.6")
Estrubolt12_7 = Rivet(89600, 12.7/1000, 57800, 50/100, "Estrubolt 12.7")
Estrubolt16_0 = Rivet(126700, 16.0/1000, 91200, 50/100, "Estrubolt 16.0")

# Bolts: https://eurocodeapplied.com/design/en1993/bolt-design-properties
M8_grade_8_8_nonthread = Rivet(19300, 8/1000, 21100, 14/100, "M8 grade 8.8 nonthread")
M10_grade_8_8 = Rivet(22300, 11/1000, 33400, 17/100, "M10 grade 8.8")
M10_grade_8_8_nonthread = Rivet(30200, 11/1000, 33400, 17/100, "M10 grade 8.8 non threaded")
M14_grade_8_8 = Rivet(44200, 15/1000, 66200, 25/100, "M14 grade 8.8")
M14_grade_8_8_nonthread = Rivet(66200, 15/1000, 66200, 25/100, "M14 grade 8.8")
M20_grade_8_8_nonthread = Rivet(90500, 22/1000, 141000, 20/100, "M20 grade 8.8 nonthread")
long25mmrivet = Rivet(4000, 4.8/1000, 5000, 70/250, "Long 25 mm grip rivet ")# https://rivet-expert.com/nl/standaard-popnagel-bolkop-rvs-a2-rvs-a2-4-80x30-00-mm-klembereik-20-00-25-00-mm

# Solid rivets

# A2 rivets

solid4mmA2 = Rivet(A2yieldshear*np.pi*(3.82/2000)**2, 3.82/1000, 450e6*np.pi*(3.82/2000)**2, 80/250, "Solid 4mm rivet ")
solid5mmA2 = Rivet(A2yieldshear*np.pi*(4.82/2000)**2, 4.82/1000, 450e6*np.pi*(4.82/2000)**2, 80/250, "Solid 5mm rivet ")
solid6mmA2 = Rivet(A2yieldshear*np.pi*(5.82/2000)**2, 5.82/1000, 450e6*np.pi*(5.82/2000)**2, 80/250, "Solid 6mm rivet ")
solid8mmA2 = Rivet(A2yieldshear*np.pi*(7.82/2000)**2, 7.82/1000, 450e6*np.pi*(7.82/2000)**2, 80/250, "Solid 8mm rivet ")


rivet_list_THICK = [Gesipa1433601, Gesipa1454074, Gesipa1454608, Gesipa1456041, Gesipa1433932, Gesipa1433811, Gesipa1433814, BibusCN03_6419_0700, BN_84013, M8_grade_8_8_nonthread, M10_grade_8_8, M14_grade_8_8, M20_grade_8_8_nonthread, M14_grade_8_8_nonthread, M10_grade_8_8_nonthread, solid5mmA2, solid6mmA2, solid8mmA2, solid4mmA2, Honsel_10792078225, Honsel_10792064205, Estrubolt6_4, Estrubolt8_0, Estrubolt9_6, Estrubolt12_7, Estrubolt16_0]
rivet_list_thin = [Gesipa1433601, Gesipa1454074, Gesipa1454608, Gesipa1456041, Gesipa1433932, Gesipa1433811, Gesipa1433814, BibusCN02_6420_0700, BibusCN03_6419_0700, BN_21443, BN_84012, BN_84013, M10_grade_8_8, M14_grade_8_8, M20_grade_8_8_nonthread, M14_grade_8_8_nonthread, M10_grade_8_8_nonthread, M14_grade_8_8_nonthread, M8_grade_8_8_nonthread, solid5mmA2, solid6mmA2, solid8mmA2, solid4mmA2, Honsel_10792078225, Honsel_10792064205, Estrubolt6_4, Estrubolt8_0, Estrubolt9_6, Estrubolt12_7, Estrubolt16_0]
rivet_list_SuperTHICK = [long25mmrivet]

rivet_list_ESTRUBOLTS = [ Estrubolt6_4, Estrubolt8_0, Estrubolt9_6, Estrubolt12_7, Estrubolt16_0]

only_bolt_list = [M10_grade_8_8, M14_grade_8_8, M20_grade_8_8_nonthread, M14_grade_8_8_nonthread, M10_grade_8_8_nonthread, M14_grade_8_8_nonthread]

# https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6063T6

alu6063_T6_shear_yield = 152e6
alu6063_T6_tension_yield = 214e6
alu6063_T6_eMod = 68.9e9
alu6063_T6_bearing_yield = 276e6

alu6061_T6_shear_yield = 207e6
alu6061_T6_tension_yield = 276e6
alu6061_T6_eMod = 68.9e9
alu6061_T6_bearing_yield = 607e6

top_panel_thickness = 4/1000
side_panel_thickness = 5/1000
column_thickness = 10/1000

top_panel = Panel(top_panel_thickness, alu6061_T6_shear_yield, alu6061_T6_tension_yield, alu6061_T6_eMod, alu6061_T6_bearing_yield)
side_panel = Panel(side_panel_thickness, alu6061_T6_shear_yield, alu6061_T6_tension_yield, alu6061_T6_eMod, alu6061_T6_bearing_yield)
column = Panel(column_thickness, alu6061_T6_shear_yield, alu6061_T6_tension_yield, alu6061_T6_eMod, alu6061_T6_bearing_yield)

column_width = 80/1000
print("--------RIVET SIZER--------\n\n")
print("Select connection: ")
print(" 1. Column top panel")
print(" 2. Column side panel\n")

option = int(input("Enter number of choice: \n"))

print("\n")
stress = float(input("Now enter stress in Pa: "))

if option == 1:
    connection = Connection(stress, column, top_panel, column_width, SF=1.25)
    # connection.find_The_Chosen_One(rivet_list_thin, cost_multiplier=1, number_multiplier=0)
    # connection.find_The_Chosen_One(rivet_list_ESTRUBOLTS, cost_multiplier=1, number_multiplier=0)
    connection.find_The_Chosen_One([Estrubolt8_0], cost_multiplier=1, number_multiplier=0) # Max 14.7 MPa, rivet limited
    # connection.find_The_Chosen_One(rivet_list_SuperTHICK, cost_multiplier=1, number_multiplier=0)
    # connection.find_The_Chosen_One(only_bolt_list, cost_multiplier=1, number_multiplier=0)
    print("\n")

if option == 2:
    connection = Connection(stress, column, side_panel, column_width, SF=1.25)
    # connection.find_The_Chosen_One(rivet_list_THICK, cost_multiplier=1, number_multiplier=0)
    # connection.find_The_Chosen_One(rivet_list_ESTRUBOLTS, cost_multiplier=1, number_multiplier=0)
    connection.find_The_Chosen_One([Estrubolt8_0], cost_multiplier=1, number_multiplier=0) # Max 14.7 MPa, rivet limited
    # connection.find_The_Chosen_One(rivet_list_SuperTHICK, cost_multiplier=1, number_multiplier=0)
    # connection.find_The_Chosen_One(only_bolt_list, cost_multiplier=1, number_multiplier=0)
    print("\n")