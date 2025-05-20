import numpy as np
import matplotlib.pyplot as plt

# define variables
g = 9.81

## THE VARIABLES BELOW ARE TAKEN FROM NEWPANAMAX SHIP ##
L = 366 # Length [m]
T_SC = 15.2 # Draught "draft" [m]
B = 51.25 # breadth / beam [m]
GM = 0.11 * B # metacentric height [m] (since B > 40)
## no ##

k_r = 0.39 * B # roll radius of gyration [m]


f_beta = 1  

C_XG =
C_XS =
C_XP = 

C_YG = 
C_YS = 
C_YR = 

C_ZH = 
C_ZR = 
C_ZP = 

L_0 = np.max(110,L)


phi =  # Pitch angle
theta =  # Roll angle 


def get_acc(x,y,z):
    a_x = f_beta *(-C_XG * g * np.sin(phi) +C_XS * a_surge + C_XP* a_pitch * (z-R))

    a_y = f_beta *(C_YG * g * np.sin(theta) +C_YS * a_sway - C_YR* a_roll * (z-R))

    a_z = f_beta * (C_ZH * a_heave + C_ZR*a_roll * y - C_ZP * a_pitch * (x-0.45*L))


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
    