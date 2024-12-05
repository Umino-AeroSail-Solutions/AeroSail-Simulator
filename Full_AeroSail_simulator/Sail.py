from re import match
from unittest import case

from scipy.interpolate import RegularGridInterpolator
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Profile
import XFLR5_interpolarion_creator as XFLR5_interp
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
        self.airfoil.add_flap(self.chordratio, flapdeflection, reset_foil=True)

    # Returns the whole sail coefficients using the Prandtl approximation: https://webstor.srmist.edu.in/web_assets/srm_mainsite/files/downloads/class4-2012.pdf
    # RETURNS [0,0,0] IF IT COULDN'T GET THE COEFFICIENTS --> PYXFOIL RETURNED NONE
    def get_sail_coefficients(self, alpha, flapdeflection, p_interpolation=None):
        self.add_flap(flapdeflection)
        airfoil_coefficients = self.airfoil.get_coefficients(alpha, self.mach, self.re, interpolate=p_interpolation)
        if (airfoil_coefficients[0] is not None) and (airfoil_coefficients[1] is not None) and (airfoil_coefficients[2] is not None):
            profile_cl = airfoil_coefficients[0]
            profile_cd = airfoil_coefficients[1]
            self.cm = airfoil_coefficients[2]
            self.cl = profile_cl / (1 + (profile_cl / (np.pi * self.oswalde * self.ar)))
            self.cd = profile_cd + ((self.cl ** 2 )/ (np.pi * self.oswalde * self.ar))
            # print([self.cl, self.cd, self.cm])
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
    def create_interpolation(self, almin, almax, alstep, flapmin, flapmax, flapstep, p_interpolation=None):
        alphas = np.arange(almin, almax, alstep, dtype=float)
        flaps = np.arange(flapmin, flapmax, flapstep, dtype=float)
        Alphas, Flaps = np.meshgrid(alphas, flaps)
        Cl = np.zeros_like(Alphas, dtype=float)
        Cd = np.zeros_like(Alphas, dtype=float)
        CloCd = np.zeros_like(Alphas, dtype=float)
        for i in range(len(flaps)):
            for j in range(len(alphas)):
                alpha = alphas[j]
                flapdeflection = flaps[i]
                coefficients = self.get_sail_coefficients(alpha, flapdeflection, p_interpolation=p_interpolation)
                Cl[i, j] = coefficients[0]
                Cd[i, j] = coefficients[1]
                # CloCd[i, j] = Cl[i, j]/Cd[i, j]
                CloCd[i, j] = Cl[i, j] # Just for testing
                print([alpha, flapdeflection, Cl[i, j], Cd[i, j]])
        self.InterpAlphas = alphas
        self.InterpFlaps = flaps
        self.InterpCls = Cl
        self.InterpCds = Cd
        self.InterpCloCds = np.divide(Cl, Cd, out=np.zeros_like(Cl), where=Cd != 0)
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
    def plot_2d_polar_interp(self):
        fig = plt.figure(figsize=(18, 6))
        Alphas, Flaps = np.meshgrid(self.InterpAlphas, self.InterpFlaps)

        # Plot Cl
        ax1 = fig.add_subplot(131, projection='3d')
        ax1.plot_surface(Alphas, Flaps, self.InterpCls, cmap='coolwarm')
        ax1.set_xlabel('Alpha (Degrees)')
        ax1.set_ylabel('Flap Deflection (Radians)')
        ax1.set_zlabel('Cl')
        ax1.set_title('3D Interpolation of Cl')

        # Plot Cd
        ax2 = fig.add_subplot(132, projection='3d')
        ax2.plot_surface(Alphas, Flaps, self.InterpCds, cmap='plasma')
        ax2.set_xlabel('Alpha (Degrees)')
        ax2.set_ylabel('Flap Deflection (Radians)')
        ax2.set_zlabel('Cd')
        ax2.set_title('3D Interpolation of Cd')

        # Plot Cl/Cd
        ax3 = fig.add_subplot(133, projection='3d')
        ax3.plot_surface(Alphas, Flaps, self.InterpCloCds, cmap='viridis')
        ax3.set_xlabel('Alpha (Degrees)')
        ax3.set_ylabel('Flap Deflection (Radians)')
        ax3.set_zlabel('Cl/Cd')
        ax3.set_title('3D Interpolation of Cl/Cd')

        plt.tight_layout()
        plt.show()
    def create_XFLR5_interpolation(self, dir):
        self.InterpAlphas, self.InterpFlaps, self.InterpCls, self.InterpCds, self.InterpCloCds = XFLR5_interp.crt_XFLR5_interpolation(dir)


# TESTING CODE -------------------------------------------------

Profile.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')
Sail = Sail('Data/E473coordinates.txt', 5, 0.4, 30, panels = 20)
# print(Sail.get_sail_coefficients(15, np.radians(10)))
# print(Sail.get_l_d_m(10, np.radians(10), 10))
# print(Sail.get_l_d_m(0, 0, 10))
# Sail.plot_polar(-10, 20, 0.5, np.radians(15))
# Sail.create_interpolation(-10, 20, 1, np.radians(0), np.radians(20), np.radians(1), p_interpolation='Data/interp0.4profile.npz')
Sail.create_XFLR5_interpolation('Data/XFLR5_5_30_0,5_10m_s_INTERPOLATION')
Sail.save_interpolation('Data/interpolationCR4sail_XFLR5.npz')
Sail.load_interpolation('Data/interpolationCR4sail_XFLR5.npz')
print(Sail.InterpAlphas, Sail.InterpFlaps, Sail.InterpCds, Sail.InterpCloCds)
Sail.plot_2d_polar_interp()