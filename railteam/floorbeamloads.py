import numpy as np
import os
import matplotlib.pyplot as plt
from railteam.Force_In_Rail_Calculator_Andres_saves_the_day import alpha, beta, P1, P2, P3, P4, P1_ext, P4_ext
CONTAINERLENGTH = 12.024

def get_safety_factor(safety_factor = 1.5):
    return safety_factor

def initialize_rails():
    # processes forces and stuff
    andres_dic = np.load(os.path.join(os.path.dirname(__file__),"forcevalues.npz"))

    safety_factor = get_safety_factor() 

    # R1_values = andres_dic["R1_values"]
    # R2_values = andres_dic["R2_values"]
    l_values = andres_dic["l_values"]
    l2_values = andres_dic["l2_values"]
    support_forces_bot = andres_dic["support_forces_bot"]
    support_forces_top = andres_dic["support_forces_top"]

    L_Top = max(l2_values)
    L_Bot = max(l_values)

    # R1_values = [i / 2 * safety_factor for i in R1_values ] # because it's a list.
    # R2_values = [i / 2 * safety_factor for i in R2_values ] #         
    # R1_max = np.max(np.abs(np.array(R1_values)))
    # R2_max = np.max(np.abs(np.array(R2_values)))

    support_forces_bot = np.transpose(np.array([i/2*safety_factor for i in support_forces_bot]))
    support_forces_top = np.transpose(np.array([i/2*safety_factor for i in support_forces_top]))
    # print("Bottom support forces. Ideally I have Left in one and right in one array")
    # print(support_forces_bot)
    # print("Top Support forces.")
    # print(support_forces_top)

    # Everything is multiplied by -1 because the support forces are now applied forces on the things that hold the rails.
    # Support forces bottom left
    sup_b_l = -1*support_forces_bot[0]
    # Support forces bottom right
    sup_b_r = -1*support_forces_bot[1]
    # Support forces top left
    sup_t_l = -1*support_forces_top[0]
    # Support forces top right
    sup_t_r = -1*support_forces_top[1]
    # print("Bot Left, Should be zero at the end")
    # print(sup_b_l[-1])
    # print("Bot right, should be zero at the beginning")
    # print(sup_b_r[0])
    # print("Top Left, should be zero at the end")
    # print(sup_t_l[-1])
    # print("Top right, should be zero at the beginning")
    # print(sup_t_r[0])

    sup_max_idx  = {
        "bl" : np.argmax(np.abs(sup_b_l)), # max value, index
        "br" : np.argmax(np.abs(sup_b_r)), # max value, index
        "tl" : np.argmax(np.abs(sup_t_l)), # max value, index
        "tr" : np.argmax(np.abs(sup_t_r)), # max value, index
    }

    # print(sup_b_l[sup_max_idx["bl"]])
    # print(sup_b_r[sup_max_idx["br"]])
    # print(sup_t_l[sup_max_idx["tl"]])
    # print(sup_t_r[sup_max_idx["tr"]])
    
    return sup_b_l, sup_b_r, sup_t_l, sup_t_r, sup_max_idx, L_Top, L_Bot 



def frame_elements():
    frame_p_A1 = np.array([0,2])
    frame_p_A2 = P2

    frame_p_B1 = P3
    frame_p_B2 = np.array([P4_ext[0], P3[1]])

    frame_p_C = P1_ext

    frame_p_D = P4_ext

    return frame_p_A1, frame_p_A2, frame_p_B1, frame_p_B2, frame_p_C, frame_p_D



# Assumes Pin and roller at twist locks
def analyticaltest(sup_b_l, sup_b_r, sup_t_l, sup_t_r, sup_max_idx, alpha, beta, condition):
    frame_p_A1, frame_p_A2, frame_p_B1, frame_p_B2, frame_p_C, frame_p_D = frame_elements()

    # Element A
    gamma1 = -np.arctan2((frame_p_A2[1]-frame_p_A1[1]),(frame_p_A2[0]-frame_p_A1[0]))
    l_A = np.linalg.norm(frame_p_A2-frame_p_A1)

    sup_t_l_max = sup_t_l[sup_max_idx[condition]] 
    F1 = -sup_t_l_max*np.array([-np.sin(beta), np.cos(beta)])
    M1 = -sup_t_l_max*np.sin(np.pi/2-beta-gamma1)*l_A

    # A-Wall intesection (Applied load on the left wall)
    F1 = -F1
    M1 = -M1
    # print(F1, M1)

    # Element/Point C (Applied load on the left wall)
    sup_b_l_max = sup_b_l[sup_max_idx[condition]]
    F2 = sup_b_l_max*np.array([np.sin(alpha), np.cos(alpha)])
    # print(F2)

    # Left Bottom (Reaction Loads)
    C1 = -F1-F2
    # print("check moment directions. SHould be neg, neg.")
    # print(-F1[0]*frame_p_A1[1], -F2[0]*frame_p_C[1])
    c1_M = -M1+F1[0]*frame_p_A1[1]+F2[0]*frame_p_C[1]
    print(C1, c1_M)
    
    # Applied loads on the floor beam
    C1_app =  -C1
    c1_M_app = -c1_M


    # Element B
    gamma2 = -np.arctan2((frame_p_B2[1]-frame_p_B1[1]),(frame_p_B2[0]-frame_p_B1[0]))
    l_B = np.linalg.norm(frame_p_B2-frame_p_B1)

    sup_t_r_max = sup_t_r[sup_max_idx[condition]]
    F3 = -sup_t_r_max*np.array([-np.sin(beta), np.cos(beta)])
    M3 = sup_t_r_max*np.sin(np.pi/2-beta-gamma2)*l_B

    # B-Wall interfacing (Applied load on the right wall)
    F3 = -F3
    M3 = -M3
    # print(gamma2, F3, M3)

    # Element/Point D (Applied load on the right wall)
    sup_b_r_max = sup_b_r[sup_max_idx[condition]]
    F4 = sup_b_r_max*np.array([np.sin(alpha), np.cos(alpha)])
    # print(F4)

    # Right Bottom (Reaction Loads)
    C2 = -F3-F4
    # print("check moment directions. SHould be neg, zero.")
    # print(-F3[0]*frame_p_B2[1], -F4[0]*frame_p_D[1])
    # print("Should be zero", frame_p_D[1])
    c2_M = -M3+F3[0]*frame_p_B2[1]+F4[0]*frame_p_D[1]
    print(C2, c2_M)

    # Applied loads on the floor beam
    C2_app =  -C2
    c2_M_app = -c2_M

    
    # Floor Beam Equilibrium
    tw1 = np.array([0.0, 0.0])
    tw2 = np.array([0.0, 0.0])

    # x-direection
    tw1[0] = -C1_app[0]-C2_app[0]

    # y-direction
    # moment around tw1
    tw2[1] = -1/CONTAINERLENGTH*(c1_M_app+c2_M_app+C2_app[1]*frame_p_D[0])
    tw1[1] = -1*(tw2[1] + C1_app[1] + C2_app[1])


    return tw1, tw2, C1_app, c1_M_app, C2_app, c2_M_app
    




def floorbeamintloads(tw1, tw2, C1_app, c1_M_app, C2_app, c2_M_app):
    frame_p_A1, frame_p_A2, frame_p_B1, frame_p_B2, frame_p_C, frame_p_D = frame_elements()
    
    # Grid
    location = np.linspace(0, CONTAINERLENGTH, 100)

    # Axial loads
    N = np.zeros(100)
    # Shear loads (shearing up diagonal is positive in this case)
    V = np.zeros(100)
    # Moments (bending concave is positive)
    M = np.zeros(100)


    for i, val in enumerate(location):
        N[i] = -tw1[0]-C1_app[0]
        V[i] = -tw1[1]-C1_app[1]
        M[i] = c1_M_app-(tw1[1]+C1_app[1])*val
        if val >= frame_p_D[0]:
            N[i] -= C2_app[0]
            V[i] -= C2_app[1]
            M[i] += c2_M_app-C2_app[1]*(val-frame_p_D[0])
    

    fig, axs = plt.subplots(3, 1)
    ax1 = axs[0]
    ax1.plot(location, N, color="red")
    ax1.fill_between(location, 0, N, color="red", alpha=0.3)
    ax1.set_xlabel("Axial")

    ax2 = axs[1]
    ax2.plot(location, V, color="green")
    ax2.fill_between(location, 0, V, color="green", alpha=0.3)
    ax2.set_xlabel("Shear")
    
    ax3 = axs[2]
    ax3.plot(location, M)
    ax3.fill_between(location, 0, M, color="blue", alpha=0.3)
    ax3.set_xlabel("Moment")

    plt.show()

    return N, V, M 


def main():
    sup_b_l, sup_b_r, sup_t_l, sup_t_r, sup_max_idx, L_Top, L_Bot = initialize_rails()
    print(f"Bottom Rail part Length: {L_Bot}, Top Length: {L_Top}")
    print("horizontal to bottom, CW:", alpha, "Horizontal to top, CCW:", beta)
    tw1, tw2, C1_app, c1_M_app, C2_app, c2_M_app = analyticaltest(sup_b_l, sup_b_r, sup_t_l, sup_t_r, sup_max_idx, alpha, beta, "tr")
    print(tw1, tw2, C1_app, C2_app, c1_M_app, c2_M_app)
    floorbeamintloads(tw1, tw2, C1_app, c1_M_app, C2_app, c2_M_app)

if __name__ == "__main__":
    main()