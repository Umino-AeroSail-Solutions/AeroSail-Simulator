import numpy as np
import matplotlib.pyplot as plt

# define variables
g = 9.81 # you already know it [m/s^2]
rho = 1025 # sea water density [kg/m^3]

## THE VARIABLES BELOW ARE TAKEN FROM NEWPANAMAX SHIP ##
L = 366 # Length [m]
T_SC = 15.2 # Draught "draft" [m]
B = 51.25 # breadth / beam [m]
GM = 0.11 * B # metacentric height [m] (since B > 40)
D = 18.3 # Moulded Depth [m] depth of the locks? according to wikipedia
m_max = 120_000_000 # deadweight of new panamax [kg]
m_min = 52_500_000 # deadweight of panamax [kg]
## -------------------------------------------------- ##

k_r = 0.39 * B # roll radius of gyration [m]
L_0 = np.max(110,L) 
L_1 = np.min(250,L)

## f coeffs
f_beta = 1  
f_ps = 0.8 # From figure 3-1, assuming 35 kts
f_p = f_ps # unless fatigue based
f_BK = 1 # bilge keel factor 
f_T = 1 # ratio between draught at a loading condition and scantling draught

T_theta = 2.3 * np.pi * k_r / np.sqrt(g * GM) # Roll period
theta = 9000 * (1.4 - 0.035*T_theta)*f_p*f_BK/(1.15*B+55)/np.pi # Roll angle

T_phi = np.sqrt(2*np.pi*(0.6*L*(1+f_T)) / g)  # Pitch period
phi = 920 * f_p * L**-0.84 * (1+2.57/np.sqrt(g*L)**1.2)   # Pitch angle

C_B = m_max/rho / (L * B * T_SC) # block coefficient = volume displacement / volume of block  BETWEEN 0.6-0.8 
R = np.min(D/4 + T_SC/2,D/2)


# accelerations
a_0 = (1.58 - 0.47 * C_B) *( 2.4/np.sqrt(L) + 34/L + 600 / L**2)
a_surge = 0.2 * (1.6 + 1.5/np.sqrt(g*L))*f_p*a_0*g # acc in x direction
a_sway = 0.3 * (2.25 - 20/np.sqrt(g*L))*f_p*a_0*g # acc in y direction
a_heave = (1.15 - 6.5/np.sqrt(g*L))*f_p*a_0*g # acc in z direction, ONLY FOR L > 150!!

a_roll = f_p * theta * np.pi/180 * (2*np.pi/T_theta)**2
a_pitch = f_p * (1.75 - 22/np.sqrt(g*L)) * phi * np.pi/180 * (2*np.pi/T_phi)**2 # Only for L > 150


C_XG =0
C_XS =0
C_XP = 

C_YG = 
C_YS = 
C_YR = 

C_ZH = 
C_ZR = 
C_ZP = 






def get_acc(x,y,z):
    a_x = f_beta *(-C_XG * g * np.sin(phi) +C_XS * a_surge + C_XP* a_pitch * (z-R))

    a_y = f_beta *(C_YG * g * np.sin(theta) +C_YS * a_sway - C_YR* a_roll * (z-R))

    a_z = f_beta * (C_ZH * a_heave + C_ZR*a_roll * y - C_ZP * a_pitch * (x-0.45*L))
    
    return [a_x,a_y,a_z]


def get_acc_env(x,y,z):
    if L < 90:
        f_L = 1
    elif L < 150:
        f_L = 1.3 - L/300
    else: f_L = 0.8

    a_x_env = 0.7 * f_L * (0.65 + 2*z/7/T_SC)*np.sqrt(a_surge**2 + L_0/325 * (g*np.sin(phi)+a_pitch*(z-R))**2)
   
    a_y_env =  (1-np.exp(-B*L/215*GM))*np.sqrt(a_sway**2+(g*np.sin(theta)+a_roll*(z-R))**2)

    a_z_env = np.sqrt(a_heave**2 + ((0.95+np.exp(-L/15))*a_pitch*(1.08*x-0.45*L))**2 + (1.2*a_roll*y)**2)

    a_z_env_pitch = np.sqrt(a_heave**2 + ((0.95+np.exp(-L/15))*a_pitch*(1.08*x-0.45*L))**2)
    
    a_z_env_roll = np.sqrt(a_heave**2 + (1.2*a_roll*y)**2)

    envelope = []
    # load cases
    envelope.append([a_x_env, 0, a_z_env_pitch])
    envelope.append([a_x_env, 0, -a_z_env_pitch])
    envelope.append([0, a_y_env, a_z_env_roll])
    envelope.append([0, a_y_env, -a_z_env_roll])
    envelope.append([0.6 * a_x_env, 0.6 * a_y_env, a_z_env])
    envelope.append([0.6 * a_x_env, 0.6 * a_y_env, -a_z_env])

    return envelope
    
