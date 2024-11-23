from re import match
from unittest import case

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Profile

# This Sail class computes atributes for a finite wing from profile parameters. It also works with flaps :)


class Sail():
    def __init__(self, plainfoil, chord, chordratio, height=None, oswalde = 1, panels=160):
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
        self.panels = panels
        self.airfoil.set_panels(panels)

    # Sets a parameter (set in the string) to a value (set as an input)
    def set_p(self, parameter, value):
        match parameter:
            case 'chord':
                self.chord = value
                self.area = self.height * self.chord
                self.ar = self.height / self.chord
            case 'chordratio':
                self.chordratio = value
            case 'height':
                self.height = value
                self.area = self.height * self.chord
                self.ar = self.height / self.chord
            case 'oswalde':
                self.oswalde = value
            case 'plainfoil':
                self.plainfoil = value
                self.airfoil = Profile.Profile(self.plainfoil)
            case 're':
                self.re = value
            case 'mach':
                self.mach = value
            case 'panels':
                self.panels = value
                self.airfoil.set_panels(self.panels)

    # Returns a specified parameter
    def get_p(self, parameter):
        match parameter:
            case 'chord':
                return self.chord
            case 'chordratio':
                return self.chordratio
            case 'height':
                return self.height
            case 'oswalde':
                return self.oswalde
            case 'plainfoil':
                return self.plainfoil
            case 're':
                return self.re
            case 'mach':
                return self.mach
            case 'panels':
                return self.panels
            case 'airfoil':
                return self.airfoil
            case 'ar':
                return self.ar
            case 'area':
                return self.area

    # Adds a flap to the entire wing at a certain deflection
    def add_flap(self, flapdeflection):
        self.airfoil = Profile.Profile(self.plainfoil)
        self.airfoil.add_flap(self.chordratio, flapdeflection)

    # Returns the whole sail coefficients using the Prandtl approximation: https://webstor.srmist.edu.in/web_assets/srm_mainsite/files/downloads/class4-2012.pdf
    # RETURNS [0,0,0] IF IT COULDN'T GET THE COEFFICIENTS --> PYXFOIL RETURNED NONE
    def get_sail_coefficients(self, alpha, flapdeflection):
        self.add_flap(flapdeflection)
        airfoil_coefficients = self.airfoil.get_coefficients(alpha, self.mach, self.re)
        if (airfoil_coefficients[0] is not None) and (airfoil_coefficients[1] is not None) and (airfoil_coefficients[2] is not None):
            profile_cl = airfoil_coefficients[0]
            profile_cd = airfoil_coefficients[1]
            self.cm = airfoil_coefficients[2]
            self.cl = profile_cl / (1 + (profile_cl / (np.pi * self.oswalde * self.ar)))
            self.cd = profile_cd + ((self.cl ** 2 )/ (np.pi * self.oswalde * self.ar))
            print([self.cl, self.cd, self.cm])
            return [self.cl, self.cd, self.cm]
        else:
            print("FAILURE")
            return [0,0,0]

    # Returns the whole sail lift and drag using the area
    def get_l_d_m(self, alpha, flapdeflection, V, rho=1.225):
        self.get_sail_coefficients(alpha, flapdeflection)
        q = 0.5 * rho * (V ** 2)
        self.l_d_m = [q * self.cl * self.area, q * self.cd * self.area, q * self.cm * self.area * self.chord]
        return self.l_d_m

    # Plots a polar of the sail, doesn't work at big angles of attack. Returns lists with coefficients
    def plot_polar(self, almin, almax, alint, flapdeflection):
        alphas = np.arange(almin, almax, alint)
        cl = np.zeros_like(alphas, dtype=float)
        cd = np.zeros_like(alphas, dtype=float)
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
        return alphas, cl, cd

    # Creates some arrays with polar values and saves it in local lists
    def create_interpolation(self, almin, almax, alstep, flapmin, flapmax, flapstep):
        alphas = np.arange(almin, almax, alstep)
        flaps = np.arange(flapmin, flapmax, flapstep)
        Alphas, Flaps = np.meshgrid(alphas, flaps)
        Cl = np.zeros_like(Alphas)
        Cd = np.zeros_like(Alphas)
        CloCd = np.zeros_like(Alphas)
        for i in range(len(flaps)):
            for j in range(len(alphas)):
                alpha = alphas[j]
                flapdeflection = flaps[i]
                coefficients = self.get_sail_coefficients(alpha, flapdeflection)
                Cl[i, j] = coefficients[0]
                Cd[i, j] = coefficients[1]
                CloCd[i, j] = Cl[i, j]/Cd[i, j]
        self.InterpAlphas = alphas
        self.InterpFlaps = Flaps
        self.InterpCls = Cl
        self.InterpCds = Cd
        self.InterpCloCds = CloCd
        return alphas, flaps, Cl, Cd, CloCd

    # Saves the interpolation arrays in a npz file
    def save_interpolation(self, filename):
        np.savez(filename, interpAlphas=self.InterpAlphas, interpFlaps=self.InterpFlaps, interpCls=self.InterpCls, interpCds=self.InterpCds, interpCloCds=self.InterpCloCds)

    # Loads the interpolation arrays from a npz file
    def load_interpolation(self, filename):
        npzfile = np.load(filename)
        self.InterpAlphas = npzfile['interpAlphas']
        self.InterpFlaps = npzfile['interpFlaps']
        self.InterpCls = npzfile['interpCls']
        self.InterpCds = npzfile['interpCds']
        self.InterpCloCds = npzfile['interpCloCds']

    # Plots the interpolation arrays
    def plot_2d_polar_Interp(self):
        plt.figure()
        fig1 = plt.figure()
        fig2 = plt.figure()
        fig3 = plt.figure()
        ax1 = fig1.add_subplot(111, projection='3d')
        ax2 = fig2.add_subplot(111, projection='3d')
        ax3 = fig3.add_subplot(111, projection='3d')
        surf1 = ax1.plot_surface(self.InterpAlphas, self.InterpFlaps, self.InterpCls, cmap='plasma')
        ax1.set_zlabel('Lift Coefficient (Cl)')
        plt.title('3D Surface Plot of Cl')
        surf2 = ax2.plot_surface(self.InterpAlphas, self.InterpFlaps, self.InterpCds, cmap='plasma')
        ax2.set_zlabel('Drag Coefficient (Cd)')
        plt.title('3D Surface Plot of Cd')
        surf3 = ax3.plot_surface(self.InterpAlphas, self.InterpFlaps, self.InterpCloCds, cmap='plasma')
        ax3.set_zlabel('(CloCd)')
        ax1.set_xlabel('Angle of Attack (radians)')
        ax1.set_ylabel('Flap Deflection (radians)')
        ax2.set_xlabel('Angle of Attack (radians)')
        ax2.set_ylabel('Flap Deflection (radians)')
        ax3.set_xlabel('Angle of Attack (radians)')
        ax3.set_ylabel('Flap Deflection (radians)')
        fig1.colorbar(surf1)
        fig2.colorbar(surf2)
        fig3.colorbar(surf3)
        plt.show()


# TESTING CODE -------------------------------------------------

Profile.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')
Sail = Sail('Data/E473coordinates.txt', 5, 0.4, 30, panels = 20)
# print(Sail.get_sail_coefficients(15, np.radians(10)))
# print(Sail.get_l_d_m(10, np.radians(10), 10))
# print(Sail.get_l_d_m(0, 0, 10))
# Sail.plot_polar(-10, 20, 0.5, np.radians(15))
# Sail.create_interpolation(-10, 20, 0.5, np.radians(0), np.radians(10), np.radians(1))
# Sail.save_interpolation('Data/test_interpolation.npz')
Sail.load_interpolation('Data/test_interpolation.npz')
Sail.plot_2d_polar_Interp()