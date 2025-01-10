from re import match
from unittest import case

from scipy.interpolate import RegularGridInterpolator, griddata
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Profile
import XFLR5_interpolarion_creator as XFLR5_interp
# This Sail_Class class computes atributes for a finite wing from profile parameters. It also works with flaps :)


class Sail_Class():
    def __init__(self, plainfoil, chord, chordratio, height=None, oswalde = 1, panels=160):
        self.chord = chord
        self.chordratio = chordratio
        self.height = height
        self.oswalde = oswalde
        self.plainfoil = plainfoil
        self.airfoil = Profile.Profile(self.plainfoil)
        if height is not None:
            self.ar = self.height/self.chord
            self.area = self.height * self.chord
        self.re = 1000000
        self.mach = 0
        self.l_d_m = [None, None, None]
        self.cl, self.cd, self.cm = 0,0,0
        self.panels = panels
        self.airfoil.set_panels(panels)
        self.yaw = 0.0
    def set_p(self, parameter, value):
        '''Sets a parameter (set in the string) to a value (set as an input)'''
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
    def get_p(self, parameter):
        '''Returns a specified parameter'''
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
    def add_flap(self, flapdeflection):
        '''Adds a flap to the entire wing at a certain deflection'''
        self.airfoil = Profile.Profile(self.plainfoil)
        self.airfoil.add_flap(self.chordratio, flapdeflection, reset_foil=True)
    def get_sail_coefficients(self, alpha, flapdeflection, p_interpolation=None, s_interpolation=None):
        '''Returns the whole sail coefficients using the Prandtl approximation: https://webstor.srmist.edu.in/web_assets/srm_mainsite/files/downloads/class4-2012.pdf
        RETURNS [0,0,0] IF IT COULDN'T GET THE COEFFICIENTS --> PYXFOIL RETURNED NONE'''
        if s_interpolation is None:
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
        else:
            self.load_interpolation(s_interpolation)
            Alphas, Flaps = np.meshgrid(self.InterpAlphas, self.InterpFlaps)
            # Check and transpose arrays if their shapes do not match
            if Alphas.shape != self.InterpCls.shape:
                self.InterpCls = self.InterpCls.T
            if Alphas.shape != self.InterpCds.shape:
                self.InterpCds = self.InterpCds.T
            if Alphas.shape != self.InterpCloCds.shape:
                self.InterpCloCds = self.InterpCloCds.T
            points = np.vstack((Alphas.flatten(), Flaps.flatten())).T
            self.cl = griddata(points, self.InterpCls.flatten(), (alpha, flapdeflection), method='linear').item()
            self.cd = griddata(points, self.InterpCds.flatten(), (alpha, flapdeflection), method='linear').item()
            # self.cm = griddata(points, self.InterpCms.flatten(), (alpha, flapdeflection), method='linear').item()
            self.cm = 0
            return [self.cl, self.cd, self.cm]
    def get_l_d_m(self, alpha, flapdeflection, V, rho=1.225, p_interpolation=None, s_interpolation=None):
        '''Returns the whole sail lift and drag using the area'''
        self.get_sail_coefficients(alpha, flapdeflection, p_interpolation=p_interpolation, s_interpolation=s_interpolation)
        q = 0.5 * rho * (V ** 2)
        self.l_d_m = [q * self.cl * self.area, q * self.cd * self.area, q * self.cm * self.area * self.chord]
        return self.l_d_m
    def plot_polar(self, almin, almax, alint, flapdeflection):
        '''Plots a polar of the sail, doesn't work at big angles of attack. Returns lists with coefficients'''
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
    def create_interpolation(self, almin, almax, alstep, flapmin, flapmax, flapstep, p_interpolation=None):
        '''Creates some arrays with polar values and saves it in local lists'''
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
    def save_interpolation(self, filename):
        '''Saves the interpolation arrays in a npz file'''
        np.savez(filename, interpAlphas=self.InterpAlphas, interpFlaps=self.InterpFlaps, interpCls=self.InterpCls, interpCds=self.InterpCds, interpCloCds=self.InterpCloCds)
    def load_interpolation(self, filename):
        '''Loads the interpolation arrays from a npz file'''
        npzfile = np.load(filename)
        self.InterpAlphas = npzfile['interpAlphas']
        self.InterpFlaps = npzfile['interpFlaps']
        self.InterpCls = npzfile['interpCls']
        self.InterpCds = npzfile['interpCds']
        self.InterpCloCds = npzfile['interpCloCds']
    def plot_2d_polar_interp(self):
        '''Plots the interpolation arrays'''
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
        '''Creates an interpolation using XFLR5 files in a directory, ending in "Flap-XX.txt" where XX is the flap deflection'''
        self.InterpAlphas, self.InterpFlaps, self.InterpCls, self.InterpCds, self.InterpCloCds = XFLR5_interp.crt_XFLR5_interpolation(dir)
    def get_cf(self, plot=False):
        '''Gets the force coefficient arrays'''
        self.cf = np.sqrt(np.square(self.InterpCls) + np.square(self.InterpCds))
        Alphas, Flaps = np.meshgrid(self.InterpAlphas, self.InterpFlaps)
        if plot:
            fig2 = plt.figure()
            cfplot = fig2.add_subplot(111, projection='3d')
            cfplot.plot_surface(Alphas,Flaps,self.cf, cmap='plasma')
            cfplot.set_xlabel('Alpha (Degrees)')
            cfplot.set_ylabel('Flap Deflection (Radians)')
            cfplot.set_zlabel('Cf')
            cfplot.set_title('3D Interpolation of Cf')

            plt.show( )
            print("Max Cf = " + str(np.max(self.cf)))
        return np.max(self.cf)
    def get_cts(self, AWA):
        '''Gets the thrust coefficient array'''
        self.cts = np.add(np.multiply(self.InterpCls, np.sin(AWA)), np.multiply(self.InterpCds, -1*np.cos(AWA)))
        return self.cts
    def get_opt_pos(self, AWA):
        '''Finds the optimum flap deflection and alpha for an Apparent Wind Angle AWA'''
        if AWA > np.pi:
            AWA -= 2*np.pi
        self.get_cts(abs(AWA))
        maxloc = np.unravel_index(np.argmax(self.cts), self.cts.shape)
        if AWA >= 0:
            self.opt_flap = self.InterpFlaps[maxloc[0]]
            self.opt_alpha = self.InterpAlphas[maxloc[1]]
        else:
            self.opt_flap = -self.InterpFlaps[maxloc[0]]
            self.opt_alpha = -self.InterpAlphas[maxloc[1]]
        return self.opt_alpha, self.opt_flap
    def plot_cts_for_AWA(self, AWA):
        '''Plots the thrust coefficient arrays'''
        # Calculate cts
        self.get_cts(AWA)
        opt_alpha, opt_flap = self.get_opt_pos(AWA)

        # Create a meshgrid for plotting
        Alphas, Flaps = np.meshgrid(self.InterpAlphas, self.InterpFlaps)

        # Plotting cts values
        plt.figure(figsize=(10, 6))
        contour = plt.contourf(Alphas, Flaps, self.cts, cmap='viridis')
        plt.colorbar(contour, label='cts')
        plt.xlabel('Alpha (Degrees)')
        plt.ylabel('Flap Deflection (Radians)')
        plt.title(f'cts for AWA = {(np.degrees(AWA)):.2f} radians')

        # Plotting the optimal point
        plt.plot(opt_alpha, opt_flap, 'ro', label='Optimal Solution')
        plt.legend()

        plt.show()
    def plot_optimal_values(self, AWA_range):
        '''Plots the optimum flap deflection and alpha for an Apparent Wind Angle (AWA) range'''
        optimal_alphas = []
        optimal_flaps = []
        for AWA in AWA_range:
            opt_alpha, opt_flap = self.get_opt_pos(AWA)
            optimal_alphas.append(abs(opt_alpha))
            optimal_flaps.append(abs(opt_flap))

        plt.figure(figsize=(10, 6))
        plt.plot(np.degrees(AWA_range), optimal_alphas, label='Optimal Alpha')
        plt.plot(np.degrees(AWA_range), np.degrees(optimal_flaps), label='Optimal Flap Deflection')
        plt.xlabel('AWA (Degrees)')
        plt.ylabel('Value')
        plt.title('Optimal Alpha and Flap Deflection vs AWA')
        plt.legend()
        plt.grid(True)
        plt.show()
    def plot_optimal_values_polar(self, AWA_range):
        '''Plots the optimum flap deflection and alpha for an Apparent Wind Angle (AWA) range in polar graphs'''
        optimal_alphas = []
        optimal_flaps = []
        optimal_cts = []

        for AWA in AWA_range:
            opt_alpha, opt_flap = self.get_opt_pos(AWA)
            optimal_alphas.append(abs(opt_alpha))
            optimal_flaps.append(abs(np.degrees(opt_flap)))  # Convert flap deflection to degrees
            optimal_cts.append(self.cts.max())

        optimal_alphas = np.array(optimal_alphas)
        optimal_flaps = np.array(optimal_flaps)
        optimal_cts = np.array(optimal_cts)

        plt.figure(figsize=(18, 6))

        # Polar plot for optimal alpha
        ax1 = plt.subplot(131, polar=True)
        ax1.plot(AWA_range, optimal_alphas, 'b-', label='Optimal Alpha')
        ax1.set_theta_zero_location('N')
        ax1.set_theta_direction(-1)
        ax1.set_rlabel_position(-22.5)  # Move radial labels to prevent overlap
        ax1.set_title('Optimal Alpha vs AWA', va='bottom')
        ax1.legend(loc='upper right')

        # Polar plot for optimal flap deflection (in degrees)
        ax2 = plt.subplot(132, polar=True)
        ax2.plot(AWA_range, optimal_flaps, 'r-', label='Optimal Flap Deflection (Degrees)')
        ax2.set_theta_zero_location('N')
        ax2.set_theta_direction(-1)
        ax2.set_rlabel_position(-22.5)  # Move radial labels to prevent overlap
        ax2.set_title('Optimal Flap Deflection vs AWA (Degrees)', va='bottom')
        ax2.legend(loc='upper right')

        # Polar plot for optimal ct
        ax3 = plt.subplot(133, polar=True)
        ax3.plot(AWA_range, optimal_cts, 'g-', label='Optimal ct')
        ax3.set_theta_zero_location('N')
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(-22.5)  # Move radial labels to prevent overlap
        ax3.set_title('Optimal ct vs AWA', va='bottom')
        ax3.legend(loc='upper right')

        plt.tight_layout()
        plt.show()
    def get_max_ct(self):
        max_ct = 0
        for AWA in range(180):
            ct = self.get_cts(AWA).max()
            if ct > max_ct:
                max_ct = ct
        return max_ct

    def plot_cf_level_curve(self, cf_value, ax=None, title='Level Curve for Cf', filling=True):
        '''Plots a level curve for a certain value of cf as a function of alpha and flap deflection, with shading for smaller values.'''
        Alphas, Flaps = np.meshgrid(self.InterpAlphas, self.InterpFlaps)
        Cf = self.get_cf()  # Ensure cf is calculated

        if ax is None:
            fig, ax = plt.subplots()

        contour = None

        # Create a filled contour plot for shading
        if filling:
            ax.contourf(Alphas, Flaps, self.cf, levels=np.linspace(np.min(self.cf), cf_value, 100), cmap='Blues',
                        alpha=0.5)

        # Add the contour line for the specified cf_value
        contour = ax.contour(Alphas, Flaps, self.cf, levels=[cf_value], colors='blue')
        ax.clabel(contour, inline=True, fontsize=8, fmt='%.2f')

        ax.set_xlabel('Alpha (Degrees)')
        ax.set_ylabel('Flap Deflection (Radians)')
        ax.set_title(title)
        ax.grid(True)

        return contour

    def get_specific_cf(self, alpha, flap_deflection, interpolation):
        '''Returns the specific force coefficient cf for a given alpha and flap deflection.'''
        self.add_flap(flap_deflection)  # Ensure the flap deflection is set
        coefficients = self.get_sail_coefficients(alpha, flap_deflection, s_interpolation=interpolation)
        cl = coefficients[0]
        cd = coefficients[1]
        cf = np.sqrt(cl ** 2 + cd ** 2)
        return cf

    def plot_optimal_values_polar_with_thrust(self, AWA_range, thrust_min, thrust_max, windspeedkt, radial_angles):
        '''Plots the optimal thrust for an Apparent Wind Angle (AWA) range in a high-resolution polar graph'''
        optimal_thrusts = []
        windspeed = windspeedkt / 1.944
        for AWA in AWA_range:
            opt_alpha, opt_flap = self.get_opt_pos(AWA)
            optimal_thrusts.append(self.cts.max() * 0.5 * 1.225 * (windspeed ** 2) * self.height * self.chord)

        optimal_thrusts = np.array(optimal_thrusts)
        thrust_min_array = np.array([thrust_min] * len(AWA_range))  # Array for thrust_min
        thrust_max_array = np.array([thrust_max] * len(AWA_range))  # Array for thrust_max

        # Polar plot for optimal thrust with thrust_min and thrust_max
        ax = plt.subplot(111, polar=True)
        ax.plot(AWA_range, optimal_thrusts, 'g-', label='Optimal thrust')
        ax.plot(AWA_range, thrust_min_array, 'm--', label='Requirement ASV2-STK-02b')
        ax.plot(AWA_range, thrust_max_array, 'c--', label='Requirement ASV2-STK-02')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_rlabel_position(-22.5)  # Move radial labels to prevent overlap
        ax.set_title('Optimal thrust vs AWA at ' + str(windspeedkt) + ' kt', va='bottom')
        for angle in radial_angles:
            ax.axvline(np.radians(angle), color='k', linestyle='--')  # Adding radial lines

        # Move the legend below the graph
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)

        plt.tight_layout()
        plt.show()


# TESTING CODE -------------------------------------------------
#
Profile.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')
Sail = Sail_Class('Data/E473coordinates.txt', 5, 0.4, 30, panels = 20)
# # # print(Sail.get_sail_coefficients(15, np.radians(10)))
# # # print(Sail.get_l_d_m(10, np.radians(10), 10))
# # # print(Sail.get_l_d_m(0, 0, 10))
# # # Sail.plot_polar(-10, 20, 0.5, np.radians(15))
# # # Sail.create_interpolation(-10, 20, 1, np.radians(0), np.radians(20), np.radians(1), p_interpolation='Data/interp0.4profile.npz')
# Sail.create_XFLR5_interpolation('Data/XFLR5_5_30_0,5_10m_s_INTERPOLATION')
# Sail.save_interpolation('Data/interpolationCR4sail_XFLR5.npz')
Sail.load_interpolation('Data/interpolationCR4sail_XFLR5.npz')
# #
# # # print(Sail.InterpAlphas, Sail.InterpFlaps, Sail.InterpCds, Sail.InterpCloCds)
# Sail.plot_2d_polar_interp()
# Sail.get_cf()
# # Sail.plot_cts_for_AWA(np.radians(5))
# # # print(Sail.get_opt_pos(np.radians(10)))
# # # print(Sail.get_opt_pos(np.radians(120)))
# # # Sail.plot_optimal_values(np.arange(np.radians(5), np.radians(180), np.radians(0.5)))
# # # Sail.plot_optimal_values_polar(np.arange(np.radians(180), np.radians(180), np.radians(0.5)))
# Sail.load_interpolation('Data/interpolationCR4sail_XFLR5.npz')
# Sail.plot_optimal_values(np.arange(np.radians(5), np.radians(180), np.radians(0.01)))
# Sail.plot_optimal_values_polar(np.arange(np.radians(-180), np.radians(180), np.radians(0.01)))
Sail.plot_optimal_values_polar_with_thrust(np.arange(np.radians(-180), np.radians(180), np.radians(0.01)), 5000, 15000, 20, [30,150])