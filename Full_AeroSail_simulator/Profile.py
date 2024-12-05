import numpy as np
import matplotlib.pyplot as plt
from numpy import dtype
from pyxfoil import Xfoil, set_workdir, set_xfoilexe
from scipy.interpolate import griddata

# Creates a Profile class to compute Cl Cd Cm coefficients to a manipulable profile started from a standard DAT file

def initializeXfoil(workdir, xfoilexe):
    set_workdir(workdir)
    set_xfoilexe(xfoilexe)

class Profile():
    def __init__(self, PlainDATfile):
        self.PlainDAT = np.genfromtxt(PlainDATfile, delimiter=' ')
        self.x = self.PlainDAT[:,0]
        self.y = self.PlainDAT[:,1]
        self.panels = None
        self.createXfoil_foil()
        self.flapdeflection = 0
        self.chordratio = 0.4
        self.cp = [] # Alpha, mach, re, results
        # print(self.PlainDAT)

    # Defines the number of panels used in simulations, set by defect to 160, the Xfoil standard
    def set_panels(self, panels):
        self.panels = panels

    # Adds a flap at the end of the profile with a certain chord ratio to the total chord and a deflection. It can save it to an optional target file
    # NOTE: If a flap is added and another flap is added over it then a double flap will be created unless reset foil is set to true
    def add_flap(self, chordratio, radiansdeflection, optionaltargetfile=None, reset_foil=False):
        # print(f"Initial deflection: {self.flapdeflection}, New deflection: {radiansdeflection}")

        # WARNING: THIS ONLY WORKS UP TO A DEFLECTION AROUND 30 DEGREES
        self.chordratio = chordratio
        if reset_foil:
            self.x = self.PlainDAT[:, 0]
            self.y = self.PlainDAT[:, 1]
            self.flapdeflection = radiansdeflection
            # print("Foil reset.")
        else:
            self.flapdeflection += radiansdeflection

        # print(f"Flap deflection after update: {self.flapdeflection}")

        self.flappedProfile = self.PlainDAT.copy()
        for i in range(np.size(self.x, axis=0)):
            rotatedpointx = (1 - chordratio) + (
                        np.cos(-self.flapdeflection) * (self.x[i] - (1 - chordratio)) - np.sin(-self.flapdeflection) *
                        self.y[i])
            rotatedpointy = (
                        np.sin(-self.flapdeflection) * (self.x[i] - (1 - chordratio)) + np.cos(-self.flapdeflection) *
                        self.y[i])

            if self.flappedProfile[i, 0] > (1 - chordratio):
                if self.flapdeflection > 0:
                    if rotatedpointy < self.y[i]:
                        self.flappedProfile[i, 0] = rotatedpointx
                        self.flappedProfile[i, 1] = rotatedpointy
                    elif self.flapdeflection > np.radians(36):
                        self.flappedProfile[i, 1] = -self.y[i]
                elif rotatedpointy > self.y[i]:
                    self.flappedProfile[i, 0] = rotatedpointx
                    self.flappedProfile[i, 1] = rotatedpointy
                elif self.flapdeflection < np.radians(-36):
                    self.flappedProfile[i, 1] = -self.y[i]

        if optionaltargetfile is not None:
            np.savetxt(optionaltargetfile, self.flappedProfile, fmt=['%.3f', '%.3f'])

        # Update xs and ys
        self.x = self.flappedProfile[:, 0]
        self.y = self.flappedProfile[:, 1]

        # Debugging: Print updated points
        # print(f"Updated flap profile: \n{self.flappedProfile}")

        return self.flappedProfile

    # Computes the coefficients for a certain condition. Saves them and returns them as [Cl, Cd, Cm]
    def get_coefficients(self, alpha, mach, re, errorstep=0.1, errorange=2, interpolate=None):
        if interpolate is None:
            # Computes the coefficients, returns [Cl, Cd, Cm]
            self.createXfoil_foil()

            self.Cl, self.Cd, self.Cm = None, None, None
            polar = self.xfoil.run_polar(alpha, alpha + 0.2, 0.2, mach=mach, re=re)
            if (len(polar.cd) != 0) and (len(polar.cl) != 0) and (len(polar.cm) != 0):  # Only if it converges
                self.Cl = polar.cl[0]
                self.Cd = polar.cd[0]
                self.Cm = polar.cm[0]
            else:  # Did not converge :( --> Try and linearly interpolate the value with values around it
                solvedleft = False
                alphat = alpha
                Clleft, Cdleft, Cmleft, alphaleft = None, None, None, None

                # Find valid left point
                while not solvedleft:
                    alphat -= errorstep
                    if alphat < alpha - errorange:
                        break
                    polar = self.xfoil.run_polar(alphat, alphat + 0.2, 0.2, mach=mach, re=re)
                    if (len(polar.cd) != 0) and (len(polar.cl) != 0) and (len(polar.cm) != 0):  # Only if it converges
                        Clleft = polar.cl[0]
                        Cdleft = polar.cd[0]
                        Cmleft = polar.cm[0]
                        alphaleft = alphat
                        solvedleft = True

                solvedright = False
                alphat = alpha
                Clright, Cdright, Cmright, alpharight = None, None, None, None

                # Find valid right point
                while not solvedright:
                    alphat += errorstep
                    if alphat > alpha + errorange:
                        break
                    polar = self.xfoil.run_polar(alphat, alphat + 0.2, 0.2, mach=mach, re=re)
                    if (len(polar.cd) != 0) and (len(polar.cl) != 0) and (len(polar.cm) != 0):  # Only if it converges
                        Clright = polar.cl[0]
                        Cdright = polar.cd[0]
                        Cmright = polar.cm[0]
                        alpharight = alphat
                        solvedright = True

                if solvedright and solvedleft:
                    alphas = np.array([alphaleft, alpharight])
                    cls = np.array([Clleft, Clright])
                    cds = np.array([Cdleft, Cdright])
                    cms = np.array([Cmleft, Cmright])
                    self.Cl = np.interp(alpha, alphas, cls)
                    self.Cd = np.interp(alpha, alphas, cds)
                    self.Cm = np.interp(alpha, alphas, cms)
        else:
            self.load_interpolation(interpolate)
            Alphas, Flaps = np.meshgrid(self.interpalphas, self.interpflaps)
            # Check and transpose arrays if their shapes do not match
            if Alphas.shape != self.interpCl.shape:
                self.interpCl = self.interpCl.T
            if Alphas.shape != self.interpCd.shape:
                self.interpCd = self.interpCd.T
            if Alphas.shape != self.interpCloCd.shape:
                self.interpCloCd = self.interpCloCd.T
            points = np.vstack((Alphas.flatten(), Flaps.flatten())).T
            self.Cl = griddata(points, self.interpCl.flatten(), (alpha, self.flapdeflection), method='linear').item()
            self.Cd = griddata(points, self.interpCd.flatten(), (alpha, self.flapdeflection), method='linear').item()
            self.Cm = griddata(points, self.interpCm.flatten(), (alpha, self.flapdeflection), method='linear').item()

        return [self.Cl, self.Cd, self.Cm]

    # Plots the airfoil
    def plot_foil(self):
        self.createXfoil_foil()
        self.xfoil.plot_profile(ls='-')
        plt.show()

    # Plots and stores the Cp function of the profile
    def plot_cp(self, alpha, mach, re):
        rescase = self.xfoil.run_result(alpha, mach=mach, re=re)
        self.cp.append([alpha, mach, re, rescase])
        ax = None
        ax = rescase.plot_result(yaxis='cp', ax=ax, ls='-x')
        _ = ax.legend()
        plt.show()

    # Plotes a polar curve for two variables, 'alpha', 'cl', 'cd', 'clocd', 'cm'...
    def plot_curve(self, almin, almax, alint, mach, re, xaxis, yaxis):
        ax = None
        ax = self.xfoil.run_polar(almin, almax, alint, mach=mach, re=re).plot_polar(ax=ax, xaxis=xaxis, yaxis=yaxis, ls='-o')
        _ = ax.legend()
        plt.show()

    def create_interpolation(self, almin, almax, alint, flapmin, flapmaxs, flapint, mach, re, filename=None):
        self.interpflaps = np.arange(flapmin, flapmaxs, flapint, dtype=float)
        self.interpalphas = np.arange(almin, almax, alint, dtype=float)
        Alphas, Flaps = np.meshgrid(self.interpalphas, self.interpflaps)

        self.interpCl = np.zeros((len(self.interpalphas), len(self.interpflaps)), dtype=float)
        self.interpCd = np.zeros((len(self.interpalphas), len(self.interpflaps)), dtype=float)
        self.interpCloCd = np.zeros((len(self.interpalphas), len(self.interpflaps)), dtype=float)
        self.interpCm = np.zeros((len(self.interpalphas), len(self.interpflaps)), dtype=float)

        for j in range(len(self.interpflaps)):
            self.add_flap(self.chordratio, self.interpflaps[j], reset_foil=True)
            self.plot_foil()
            for i in range(len(self.interpalphas)):
                coefficients = self.get_coefficients(self.interpalphas[i], mach,re)
                self.interpCl[i, j] = coefficients[0]
                self.interpCd[i, j] = coefficients[1]
                self.interpCloCd[i, j] = coefficients[0]/coefficients[1]
                self.interpCm[i, j] = coefficients[2]


        if filename is not None:
            self.save_interpolation(filename)



    # Creates the PyXfoil object
    def createXfoil_foil(self):
        self.xfoil = Xfoil('Flapped E473')
        self.xfoil.set_points(self.x.tolist(), self.y.tolist())
        if self.panels is not None:
            self.xfoil.set_ppar(self.panels)
        else:
            self.xfoil.set_ppar(160)

    def save_interpolation(self, filename):
        np.savez(filename, interpalphas=self.interpalphas, interpflaps=self.interpflaps, interpCl=self.interpCl, interpCd=self.interpCd, interpCloCd=self.interpCloCd, interpCm=self.interpCm)

    # Loads the interpolation arrays from a npz file
    def load_interpolation(self, filename):
        npzfile = np.load(filename)
        self.interpalphas = npzfile['interpalphas']
        self.interpflaps = npzfile['interpflaps']
        self.interpCl = npzfile['interpCl']
        self.interpCd = npzfile['interpCd']
        self.interpCloCd = npzfile['interpCloCd']
        self.interpCm = npzfile['interpCm']

    def plot_2d_polar_interp(self):
        fig = plt.figure(figsize=(18, 6))
        Alphas, Flaps = np.meshgrid(self.interpalphas, self.interpflaps)

        # Check and transpose arrays if their shapes do not match
        if Alphas.shape != self.interpCl.shape:
            self.interpCl = self.interpCl.T
        if Alphas.shape != self.interpCd.shape:
            self.interpCd = self.interpCd.T
        if Alphas.shape != self.interpCloCd.shape:
            self.interpCloCd = self.interpCloCd.T

        # Plot Cl
        ax1 = fig.add_subplot(131, projection='3d')
        ax1.plot_surface(Alphas, Flaps, self.interpCl, cmap='coolwarm')
        ax1.set_xlabel('Alpha (Degrees)')
        ax1.set_ylabel('Flap Deflection (Radians)')
        ax1.set_zlabel('Cl')
        ax1.set_title('3D Interpolation of Cl')

        # Plot Cd
        ax2 = fig.add_subplot(132, projection='3d')
        ax2.plot_surface(Alphas, Flaps, self.interpCd, cmap='plasma')
        ax2.set_xlabel('Alpha (Degrees)')
        ax2.set_ylabel('Flap Deflection (Radians)')
        ax2.set_zlabel('Cd')
        ax2.set_title('3D Interpolation of Cd')

        # Plot Cl/Cd
        ax3 = fig.add_subplot(133, projection='3d')
        ax3.plot_surface(Alphas, Flaps, self.interpCloCd, cmap='viridis')
        ax3.set_xlabel('Alpha (Degrees)')
        ax3.set_ylabel('Flap Deflection (Radians)')
        ax3.set_zlabel('Cl/Cd')
        ax3.set_title('3D Interpolation of Cl/Cd')

        plt.tight_layout()
        plt.show()

    def plot_random_interpolated_points(self, num_points=10):
        # Generate random alpha and flap values within the range
        random_alphas = np.random.uniform(self.interpalphas.min(), self.interpalphas.max(), num_points)
        random_flaps = np.random.uniform(self.interpflaps.min(), self.interpflaps.max(), num_points)

        # Interpolate the values using get_coefficients
        random_Cl, random_Cd, random_CloCd = [], [], []
        for alpha, flap in zip(random_alphas, random_flaps):
            cl, cd, cm = self.get_coefficients(alpha, 0, 1000000, interpolate='Data/interp0.4profile.npz')
            random_Cl.append(cl)
            random_Cd.append(cd)
            random_CloCd.append(cl / cd if cd != 0 else 0)

        random_Cl = np.array(random_Cl)
        random_Cd = np.array(random_Cd)
        random_CloCd = np.array(random_CloCd)

        # Plot the interpolated points
        fig = plt.figure(figsize=(18, 6))

        ax1 = fig.add_subplot(131, projection='3d')
        ax1.scatter(random_alphas, random_flaps, random_Cl, c='r', label='Random Interpolated Cl')
        ax1.set_xlabel('Alpha (Degrees)')
        ax1.set_ylabel('Flap Deflection (Radians)')
        ax1.set_zlabel('Cl')
        ax1.set_title('Random Interpolated Cl Points')
        ax1.legend()

        ax2 = fig.add_subplot(132, projection='3d')
        ax2.scatter(random_alphas, random_flaps, random_Cd, c='b', label='Random Interpolated Cd')
        ax2.set_xlabel('Alpha (Degrees)')
        ax2.set_ylabel('Flap Deflection (Radians)')
        ax2.set_zlabel('Cd')
        ax2.set_title('Random Interpolated Cd Points')
        ax2.legend()

        ax3 = fig.add_subplot(133, projection='3d')
        ax3.scatter(random_alphas, random_flaps, random_CloCd, c='g', label='Random Interpolated Cl/Cd')
        ax3.set_xlabel('Alpha (Degrees)')
        ax3.set_ylabel('Flap Deflection (Radians)')
        ax3.set_zlabel('Cl/Cd')
        ax3.set_title('Random Interpolated Cl/Cd Points')
        ax3.legend()

        plt.tight_layout()
        plt.show()


# TESTING CODE -------------------------------------------------

# initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')
# testProfile = Profile('data/E473coordinates.txt')
# testProfile.plot_foil()
# testProfile.add_flap(0.4, np.radians(5), reset_foil=True)
# testProfile.plot_foil()
# testProfile.add_flap(0.4, np.radians(12), reset_foil=True)
# testProfile.plot_foil()
# testProfile.load_interpolation('Data/interp0.4profile.npz')
# # testProfile.plot_cp(5, 0, 1000000)
# # testProfile.plot_curve(-10, 20, 0.5, 0, 1000000, 'alpha', 'cl')
# # testProfile.create_interpolation(-10, 20, 1, np.radians(0), np.radians(20), np.radians(1), 0, 1000000, 'data/interp0.4profile')
# print(testProfile.get_coefficients(10, 0, 1000000, interpolate='Data/interp0.4profile.npz'))
# testProfile.plot_2d_polar_interp()
# # testProfile.plot_random_interpolated_points(num_points=200)