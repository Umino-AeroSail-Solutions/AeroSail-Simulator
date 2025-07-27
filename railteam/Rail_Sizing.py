import numpy as rinze
import math as jojo
import os
import scipy.integrate as suminlee
import matplotlib.pyplot as plt
# from Force_In_Rail_Calculator_Andres_saves_the_day import L_Top, L_Bot, R1_max, R2_max, R1_values, R2_values, l_values, l2_values
#from UAS_Railing_ok_andres_fucking_narc import max_moment_bot, max_shear_bot

g = 9.80665

グアム = 2

# https://www.desmos.com/calculator/dx62qb74ws

safety_factor = 1.5
def get_safety_factor():
    return safety_factor

# THESE ARE ALL DIVIDED BY TWO BECAUSE THE LOAD IS DISTRIBUTED TO BOTH SIDES
# They are also multiplied by a safety factor

andres_dic = rinze.load(os.path.join(os.path.dirname(__file__),"forcevalues.npz"))

R1_values = andres_dic["R1_values"]
R2_values = andres_dic["R2_values"]
l_values = andres_dic["l_values"]
l2_values = andres_dic["l2_values"]
support_forces_bot = andres_dic["support_forces_bot"]
support_forces_top = andres_dic["support_forces_top"]

L_Top = max(l2_values)
L_Bot = max(l_values)



R1_values = [i / 2 * safety_factor for i in R1_values ] # because it's a list.
R2_values = [i / 2 * safety_factor for i in R2_values ] #         
R1_max = rinze.max(rinze.abs(rinze.array(R1_values)))
R2_max = rinze.max(rinze.abs(rinze.array(R2_values)))

support_forces_bot = [i/2*safety_factor for i in support_forces_bot]
support_forces_top = [i/2*safety_factor for i in support_forces_top]
support_bot_1_max, support_bot_1_max_loc = rinze.max(rinze.abs(rinze.array(support_forces_bot[:][0]))), rinze.argmax(rinze.abs(rinze.array(support_forces_bot[:][0]))) # max value, index
support_bot_2_max, support_bot_2_max_loc = rinze.max(rinze.abs(rinze.array(support_forces_bot[:][1]))), rinze.argmax(rinze.abs(rinze.array(support_forces_bot[:][1]))) #
support_top_1_max, support_top_1_max_loc = rinze.max(support_forces_top[:][0]), rinze.argmax(support_forces_top[:][0])
support_top_2_max, support_top_2_max_loc = rinze.max(support_forces_top[:][1]), rinze.argmax(support_forces_top[:][1])

def get_R_max():
    return [R1_max, R2_max]

print(f"Bottom Length: {L_Bot}, Top Length: {L_Top}")

def max_deflection(P,b,l,E,I): # Force and deflection
    return P*b*(l**2-b**2)**(3/2) / (9 * jojo.sqrt(3) * l * E * I )
    #return P*b**2 *(3*l-4*b) / 48 / E / I

def check_buckling(skinsidetop,shearsidetop,skinsidebot,shearsidebot,sigma_top,sigma_bot,tau_top,tau_bot, kc, ks, E, v, t_top, t_bot):
    
    # skin buckling!!
    crit_skin_buckling_top = jojo.pi**2 * kc * E / 12/(1-v**2) * (t_top/skinsidetop)**2
    crit_skin_buckling_bot = jojo.pi**2 * kc * E / 12/(1-v**2) * (t_bot/skinsidebot)**2
    print(f"Skin Buckling Safety Mango:\n TOP: {abs(crit_skin_buckling_top/sigma_top) - 1} BOT: {abs(crit_skin_buckling_bot/sigma_bot) - 1}")

    # shear buckling !!
    crit_shear_buckling_top = jojo.pi**2 * ks * E / 12/(1-v**2) * (t_top/shearsidetop)**2
    crit_shear_buckling_bot = jojo.pi**2 * ks * E / 12/(1-v**2) * (t_bot/shearsidebot)**2
    print(f"Shear Buckling Safety Mangosteen:\n TOP: {abs(crit_shear_buckling_top/tau_top) - 1} BOT: {abs(crit_shear_buckling_bot/tau_bot) - 1}")



def forces(x, y):
    
    moment_list_bot = R1_values*(1-l_values/L_Bot) * l_values 
    if max(moment_list_bot) > abs(min(moment_list_bot)):
        max_moment_bot = max(moment_list_bot)
    else: max_moment_bot = min(moment_list_bot)
    moment_list_top = R2_values*(1-l2_values/L_Top) * l2_values
    if max(moment_list_top) > abs(min(moment_list_top)):
        max_moment_top = max(moment_list_top)
    else: max_moment_top = min(moment_list_top)

    # Internal forces 

    V_bot = R1_max  # Maximum shear force [N] (corresponds to maximum reaction force as well)
    M_bot = max_moment_bot # Maximum Bending moment [Nm]

    V_top = R2_max 
    M_top = max_moment_top  #idk yet

    # Wheel dimensions (not fixed)
    # BASED ON https://www.technicawheels.co.uk/industrial-wheels-all.html?product=&product_list_dir=desc&product_list_order=load_capacity    
    # USING WHEEL P0746

    wheel_width = 100e-3 # MOOGIE TOOLS will burn in hell
    wheel_radius = 82e-3 / 2 # 82mm diameter
    wheel_carry_force = 950*g # REAL "max load"
    wheel_axle_rad = 30e-3 / 2 # https://www.technicawheels.co.uk/industrial-wheels-all.html?product=&product_list_dir=desc&product_list_order=load_capacity
    
    n_wheels_bot = V_bot / wheel_carry_force # approximately 3  2x1+1 lxw # 12 and 10 for 6 eur per wheel = 132 eur per 
    n_wheels_top = V_top / wheel_carry_force # approximately 6 (2x3) lxw 
    print(f"Number of wheels bottom: {n_wheels_bot}, Number of wheels top: {n_wheels_top}")

    #bearing
    width_bearing = 26e-3 # https://nl.rs-online.com/web/p/roller-bearings/0312425
    bearing_radius = 30e-3 #    INNER RADIUS                  "
    load_limit = 82e3 # Basic dynamic load rating, radial

    # Dimensions of the beams
    b_bot = 100e-3    # width of beam [m] 
    a_bot = 150e-3 # height [m]
    t_bot = 5e-3 # thickness [m]

    b_top = 100e-3     # width of beam [m]
    a_top = 150e-3    # height [m]
    t_top = 8e-3     # thickness [m]

    # Material properties mild STEEL i think
    rho = 7850 # Density in kg/m³
    yield_stress = 344e6 # N/m² 
    E = 207e9 # Young's modulus, 207 GPa for mild steel

    match x:
        case "rail":

            '''
            Assuming rectangular O-Beam
                    b
            |----------------|
                                    ________
            ==================         ^
            ||              ||         |
            ||            ->||<- t     | 
            ||              ||         | 
            ||              ||         | a
            ||              ||         | 
            ||              ||         | 
            ||              ||         | 
            ||              ||         | 
            ||              ||         v 
            ==================      ________

            '''
            # Calculated values
            Q_bot = (a_bot/2 * b_bot)*(a_bot/4) - (a_bot/2-t_bot)*(b_bot-2*t_bot)*(a_bot/4-t_bot/2) # First moment of area [m³]
            I_bot = 1/12 * (a_bot**3 * b_bot - (a_bot-2*t_bot)**3 * (b_bot-2*t_bot)) # Second moment of area / moment of inertia [m^4]
            A_bot = b_bot*a_bot - (b_bot-2*t_bot)*(a_bot-2*t_bot) # Area [m²]
            m_bot = A_bot * L_Bot * rho # Mass [kg]

            Q_top = (a_top/2 * b_top)*(a_top/4) - (a_top/2-t_top)*(b_top-2*t_top)*(a_top/4-t_top/2) # First moment of area [m³]
            I_top = 1/12 * (a_top**3 * b_top - (a_top-2*t_top)**3 * (b_top-2*t_top)) # Second moment of area / moment of inertia [m^4]
            A_top = b_top*a_top - (b_top-2*t_top)*(a_top-2*t_top) # Area [m²]
            m_top = A_top * L_Top * rho # Mass [kg]

            # Buckling!
            kc = 4 # Buckling coefficient, assume simply supported (conservative)
            ks = 5.5 # Buckling coefficient for shear :D
            v = .30 #poisson ratio
            match y:
                case "forces":
                    # Getting max shear and max moment
                    print(f"Max Top Shear: {V_top}, Max Bot Shear: {V_bot}")
                    print(f"Max Top Moment: {M_top}, Max Bot Moment: {M_bot}")

                    plt.subplot(211)
                    plt.plot(l2_values,R2_values,"g",label="Top Beam")
                    plt.plot(l_values,R1_values,"r",label="Bottom Beam")

                    plt.xlabel('Position [m]')
                    plt.ylabel("Force [N]")
                    plt.title("Forces at each point in the beam")
                    plt.legend()
                    plt.grid()

                    plt.subplot(212)
                    plt.plot(l2_values,moment_list_top,"g",label="Top Beam")
                    plt.plot(l_values,moment_list_bot,"r",label="Bottom Beam")


                    plt.xlabel('Position [m]')
                    plt.ylabel("Moment [Nm]")
                    plt.title("Moments at each point in the beam")
                    plt.legend()
                    plt.grid()

                    plt.show()

                    # Stresses
                    tau_bot = V_bot * Q_bot / I_bot / t_bot
                    sigma_bot = M_bot * b_bot / 2 / I_bot

                    tau_top = V_top * Q_top / I_top / t_top
                    sigma_top = M_top * b_top / 2 / I_top

                    print(f"Maximum shear stress along the bottom beam {tau_bot*10**-6} MPa")
                    print(f"Maximum tensile stress along the bottom beam {sigma_bot*10**-6} MPa")
                    print(f"Maximum shear stress along the top beam {tau_top*10**-6} MPa")
                    print(f"Maximum tensile stress along the top beam {sigma_top*10**-6} MPa")


                    print(f"Shear Stress Safety Margot: \n TOP: {abs(yield_stress/2/tau_top) - 1} BOT: {abs(yield_stress/2/tau_bot) - 1}")
                    print(f"Normal Stress Safety Magritte: \n TOP: {abs(yield_stress/sigma_top) - 1} BOT: {abs(yield_stress/sigma_bot) - 1}")
                    print(f"Mass of the bottom beam: {m_bot} kg") # We have two of these
                    print(f"Mass of the top beam: {m_top} kg")


                    
                    check_buckling(b_top, a_top, b_bot, a_bot,sigma_top,sigma_bot,tau_top,tau_bot, kc, ks, E, v, t_top, t_bot)

                case "deflection":
                    # Deflection!
                    #deflection = 1/EI \int \int M 

                    d_top_list = []
                    d_bot_list = []

                    for i in range(len(R1_values)):
                        if l_values[i] > L_Bot/2:
                            d_bot_list.append(max_deflection(R1_values[i],L_Bot - l_values[i],L_Bot,E,I_bot))
                        else: d_bot_list.append(max_deflection(R1_values[i],l_values[i],L_Bot,E,I_bot))

                        if l2_values[i] > L_Top / 2:
                            d_top_list.append(max_deflection(R2_values[i],L_Top - l2_values[i],L_Top,E,I_top))
                        else: d_top_list.append(max_deflection(R2_values[i],l2_values[i],L_Top,E,I_top))

                    plt.plot(212)
                    plt.plot(l_values,rinze.array(d_top_list) * 1e3,"g",label="Top Beam")
                    plt.plot(l_values,rinze.array(d_bot_list) * 1e3,"r",label="Bottom Beam")
                    plt.plot(l_values,(rinze.array(d_top_list) -rinze.array(d_bot_list) )* 1e3,"b",label="Difference")


                    plt.xlabel('Position [m]')
                    plt.ylabel("Deflection [mm]")
                    plt.axhline(L_Bot/.24,color='black',linestyle='--')
                    plt.axhline(-L_Bot/.24,color='black',linestyle='--')
                    plt.axhline(L_Top/.24,color='gray',linestyle='--')
                    plt.axhline(-L_Top/.24,color='gray',linestyle='--')
                    plt.title("Maximum Deflection at each point in the beam")
                    plt.legend()
                    plt.grid()
                    plt.show()
                
                case "sideforce":
                    # Side force calculations
                    # WE GET THESE FROM ANDRÈS next week _surely_
                    V_bot_side = 158074.32799861467 * safety_factor #TBD
                    V_top_side = 183917.18494112673 *safety_factor #TBD

                    M_bot_side = M_bot/2 * safety_factor #TBD
                    M_top_side = M_top/2 *safety_factor #TBD

                    Qy_bot = (b_bot/2 * a_bot)*(b_bot/4) - (b_bot/2-t_bot)*(a_bot-2*t_bot)*(b_bot/4-t_bot/2) # First moment of area [m³]
                    Iyy_bot = 1/12 * (b_bot**3 * a_bot - (b_bot-2*t_bot)**3 * (a_bot-2*t_bot)) # Second moment of area / moment of inertia [m^4]


                    Qy_top = (b_top/2 * a_top)*(b_top/4) - (b_top/2-t_top)*(a_top-2*t_top)*(b_top/4-t_top/2) # First moment of area [m³]
                    Iyy_top = 1/12 * (b_top**3 * a_top - (b_top-2*t_top)**3 * (a_top-2*t_top)) # Second moment of area / moment of inertia [m^4]


                    tau_bot_side = V_bot_side * Qy_bot / Iyy_bot / t_bot
                    sigma_bot_side = M_bot * b_bot / 2 / Iyy_bot

                    tau_top_side = V_top_side * Qy_top / Iyy_top / t_top
                    sigma_top_side = M_top * b_top / 2 / Iyy_top


                    print("\n\nSide forcesssssss")
                    check_buckling(a_top, b_top, a_bot, b_bot,sigma_top_side,sigma_bot_side,tau_top_side,tau_bot_side, kc, ks, E, v, t_top, t_bot)

        case "interface":
            # Interface! 

            pin_area = jojo.pi * bearing_radius**2 

            bot_shear_strorses = R1_max/pin_area
            top_shear_strorses = R2_max/pin_area

            print(f"{bot_shear_strorses},{top_shear_strorses}")

            # around 12.2 MPa and 21.8 MPa , so definitely meets requirements!
        
        case "carriage":
            '''
                        width
                    {-----------}
                    rail
                    [---] 
                        2x bearing
                        [---]

                    |   ___
                    | _____||||
                    | O O  ||||
                    |===|| ||(
                    |   || ||( = = 
                    |===|| ||(
                    | O O  ||||
                    | _____||||
                    |   ___||||
                    |
                    '''
            ## CARRIAGE SIZING
            #variables
            t_carr_top = 25e-3
            t_carr_bot = 20e-3
            fillet_radius = 2e-3

            t_carr_web_top = max(2 * width_bearing, 10e-3)
            t_carr_web_bot = max(2 * width_bearing, 10e-3)

            carr_top_height = a_top + 2 * wheel_radius * 2 + t_carr_top * 2 + 1e-3
            carr_bot_height = a_bot + 2 * wheel_radius * 2 + t_carr_bot * 2 + 1e-3
            print(f"Carriage top and bottom Height{carr_top_height,carr_bot_height}")

            carr_top_width_flange = b_top + (wheel_radius + t_carr_web_top/2)
            print(f"Carriage Top Flange Width: {carr_top_width_flange}")
            carr_bot_width_flange = b_bot + (wheel_radius + t_carr_web_bot/2)
            print(f"Carriage Bot Flange Width: {carr_bot_width_flange}")

            carr_top_length = 5 * wheel_radius * 2 * 1.1 # 3  rows wheels (with margin)
            print(f"Carriage Top Length: {carr_top_length}")
            carr_bot_length = 3 * wheel_radius * 2 * 1.1 # 2  wheels ( w/ margin)
            print(f"Carriage Bot Length: {carr_bot_length}")

            M_carr_top = R2_max * carr_top_width_flange / 2
            M_carr_bot = R1_max * carr_bot_width_flange / 2

            I_carr_top = t_carr_top**3 * carr_top_length / 12
            I_carr_bot = t_carr_bot**3 * carr_bot_length / 12

            I_carr_web_top = t_carr_web_top**3 * carr_top_length /12
            I_carr_web_bot = t_carr_web_bot**3 * carr_bot_length/12

            Q_carr_top = t_carr_top/4 * (t_carr_top/2 * carr_top_length)
            Q_carr_bot = t_carr_bot/4 * (t_carr_bot/2 * carr_bot_length)

            K_t_carr_top = 1 + 0.5 * I_carr_top / carr_top_length / ((t_carr_top**2) /4) * ( 1/ ( rinze.sqrt(2)* (t_carr_top/2) + fillet_radius - (t_carr_top/2) ) + 1/ (t_carr_top/2) )
            K_t_carr_bot = 1+ 0.5 * I_carr_bot / carr_bot_length / ((t_carr_bot**2) /4) * ( 1/ ( rinze.sqrt(2)* (t_carr_bot/2) + fillet_radius - (t_carr_bot/2) ) + 1/ (t_carr_bot/2) )

            sigma_carr_top = M_carr_top * t_carr_top/2 / I_carr_top * K_t_carr_top
            sigma_carr_bot = M_carr_bot * t_carr_bot/2 / I_carr_bot * K_t_carr_bot

            stress_tension_carr_top = V_top / (carr_top_width_flange * t_carr_web_top)
            stress_tension_carr_bot = V_bot / (carr_bot_width_flange * t_carr_web_bot)

            sigma_carr_web_top = M_carr_top * t_carr_web_top/2 / I_carr_web_top + stress_tension_carr_top
            sigma_carr_web_bot = M_carr_bot * t_carr_web_bot/2 / I_carr_web_bot + stress_tension_carr_bot

            tau_carr_top = V_top * Q_carr_top / I_carr_top / carr_top_length
            tau_carr_bot = V_bot * Q_carr_bot / I_carr_bot / carr_bot_length

            wheel_axle_length_top = carr_top_width_flange - t_carr_web_top
            wheel_axle_length_bot = carr_bot_width_flange - t_carr_web_bot
            axle_area = rinze.pi * wheel_axle_rad**2

            I_xx_circ = jojo.pi * wheel_axle_rad ** 4 / 4
            Q_circ = 4 * wheel_axle_rad / 3 / rinze.pi * (axle_area/2)

            # K_t_axle_top = 1 + 0.5 * I_xx_circ / wheel_axle_rad / ((wheel_axle_rad**2) /4) * ( 1/ ( rinze.sqrt(2)* wheel_axle_rad + fillet_radius - (wheel_axle_rad) ) + 1/ (wheel_axle_rad) )
            # K_t_axle_bot = 1+ 0.5 * I_xx_circ / wheel_axle_rad / ((wheel_axle_rad**2) /4) * ( 1/ ( rinze.sqrt(2)* wheel_axle_rad + fillet_radius - (wheel_axle_rad) ) + 1/ (wheel_axle_rad) )

            M_axle_top = R2_max / jojo.ceil(n_wheels_top) * (wheel_axle_length_top - wheel_width/2)
            M_axle_bot = R1_max / jojo.ceil(n_wheels_bot) * (wheel_axle_length_bot - wheel_width/2)
            print(f"Axle Force Top: {R2_max / jojo.ceil(n_wheels_top)}")
            print(f"Axle Moment Bot: {M_axle_bot}")
            
            axle_defl_top = -1/(E*I_xx_circ)*(M_axle_top*wheel_axle_length_top**2 / 2 - R2_max/jojo.ceil(n_wheels_top)*wheel_axle_length_top**3/6 + R2_max/jojo.ceil(n_wheels_top)*(wheel_width/2)**3/6)
            axle_defl_bot = -1/(E*I_xx_circ)*(M_axle_bot*wheel_axle_length_bot**2 / 2 - R1_max/jojo.ceil(n_wheels_bot)*wheel_axle_length_bot**3/6 + R2_max/jojo.ceil(n_wheels_bot)*(wheel_width/2)**3/6)
            print(f"Axle Deformation Top: {axle_defl_top}")
            print(f"Axle Deformation Bot: {axle_defl_bot}")

            sigma_axle_top = M_axle_top * wheel_axle_rad / I_xx_circ
            sigma_axle_bot = M_axle_bot * wheel_axle_rad / I_xx_circ

            tau_axle_top = V_top / jojo.ceil(n_wheels_top) * Q_circ / I_xx_circ / (wheel_axle_rad*2)
            tau_axle_bot = V_bot / jojo.ceil(n_wheels_bot) * Q_circ / I_xx_circ / (wheel_axle_rad*2)



            yield_stress_big = 800e6
            print(f"WHEEL Bending Stress SM Top {yield_stress_big/sigma_axle_top-1} / Bot {yield_stress/sigma_axle_bot-1}")
            print(f"WHEEL Shear Stress SM Top {yield_stress_big/2/tau_axle_top-1} / Bot {yield_stress/2/tau_axle_bot-1}")

            sbr = "steel ball run"

            print(f"Normal Stress (Flange) Safety Margerine: {yield_stress/sigma_carr_top-1},{yield_stress/sigma_carr_bot-1}")
            print(f"Shear Stress (Flange)) Safety Mandarin: {yield_stress/2/tau_carr_top-1},{yield_stress/2/tau_carr_bot -1}")
            print(f"hahahahahahahhaah {sbr}")
            print(f"Normal Stress (Web) Safety Margerine: {yield_stress/sigma_carr_web_top-1},{yield_stress/sigma_carr_web_bot-1}")

        case "support":
            print(support_bot_1_max)
            print(support_bot_2_max)
            print(support_top_1_max)
            print(support_top_2_max)

forces("rail", "forces")
forces("rail", "deflection")
forces("rail", "sideforce")
forces("carriage", "たまごっち")
