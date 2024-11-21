from re import match

import numpy as np
import matplotlib.pyplot as plt

import Profile

# This Sail class computes atributes for a finite wing from profile parameters. It also works with flaps :)

class Sail():
    def __init__(self, plainfoil, chord, chordratio, height=None, oswalde = 1):
        self.chord = chord
        self.chordratio = chordratio
        self.height = height
        self.oswalde = oswalde
        self.plainfoil = plainfoil
        self.airfoil = Profile.Profile(self.plainfoil)
        self.ar = self.height/self.chord
        self.re = 1000000
        self.mach = 0
        self.area = self.height * self.chord
        self.l_d_m = [None, None, None]
        self.cl, self.cd, self.cm = 0,0,0

    # Sets a parameter (set in the string) to a value (set as an input)
    def set_p(self, parameter, value):
        match parameter:
            case 'chord':
                self.chord = value
            case 'chordratio':
                self.chordratio = value
            case 'height':
                self.height = value
            case 'oswalde':
                self.oswalde = value
            case 'plainfoil':
                self.plainfoil = value
            case 're':
                self.re = value
            case 'mach':
                self.mach = value

    # Adds a flap to the entire wing at a certain deflection
    def add_flap(self, flapdeflection):
        self.airfoil = Profile.Profile(self.plainfoil)
        self.airfoil.add_flap(self.chordratio, flapdeflection)

    # Returns the whole sail coefficients using the Prandtl approximation: https://webstor.srmist.edu.in/web_assets/srm_mainsite/files/downloads/class4-2012.pdf
    def get_sail_coefficients(self, alpha, flapdeflection):
        self.add_flap(flapdeflection)
        airfoil_coefficients = self.airfoil.get_coefficients(alpha, self.mach, self.re)
        profile_cl = airfoil_coefficients[0]
        profile_cd = airfoil_coefficients[1]
        self.cm = airfoil_coefficients[2]
        self.cl = profile_cl / (1 + (profile_cl / (np.pi * self.oswalde * self.ar)))
        self.cd = profile_cd + ((self.cl ** 2 )/ (np.pi * self.oswalde * self.ar))
        return [self.cl, self.cd, self.cm]

    # Returns the whole sail lift and drag using the area
    def get_l_d_m(self, alpha, flapdeflection, V, rho=1.225):
        self.get_sail_coefficients(alpha, flapdeflection)
        q = 0.5 * rho * (V ** 2)
        self.l_d_m = [q * self.cl * self.area, q * self.cd * self.area, q * self.cm * self.area * self.chord]
        return self.l_d_m

    # Plots a polar of the sail, doesn't work at big angles of attack
    def plot_polar(self, almin, almax, alint, flapdeflection):
        alphas = np.arange(almin, almax, alint)
        cl = np.zeros_like(alphas)
        cd = np.zeros_like(alphas)
        print(alphas)
        for i, alpha in enumerate(alphas):
            coefficients = self.get_sail_coefficients(alpha, flapdeflection)
            cl[i] = coefficients[0]
            cd[i] = coefficients[1]
            print(f'Alpha: {alpha}, Cl: {cl[i]}, Cd: {cd[i]}')
        print(alphas, cl, cd)
        plt.figure()
        plt.plot(cl, cd, label=f'Flap Deflection: {np.degrees(flapdeflection)} degrees')
        plt.ylabel('Drag Coefficient (Cd)')
        plt.xlabel('Lift Coefficient (Cl)')
        plt.title('Sail Polar Plot')
        plt.legend()
        plt.grid(True)
        plt.show()

# TESTING CODE -------------------------------------------------

Profile.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')
Sail = Sail('Data/E473coordinates.txt', 5, 0.4, 30)
printprint(Sail.get_sail_coefficients(10, np.radians(10)))
print(Sail.get_l_d_m(10, np.radians(10), 10))
# print(Sail.get_l_d_m(0, 0, 10))
# Sail.plot_polar(-5, 10, 0.5, np.radians(5))