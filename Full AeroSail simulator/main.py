import numpy as np
import matplotlib.pyplot as plt

import Profile
import Sail

Profile.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')

#Class to handle environmental conditions
class Conditions():
    def __init__(self, windspeed, bearing, rho=1.225,re_per_chord=None, viscosity=1.802*(10**(-5))):
        self.windspeed = windspeed
        self.bearing = bearing
        self.rho = rho
        self.viscosity = viscosity
        if re_per_chord is None:
            self.re_per_chord = (self.windspeed * self.rho) / self.viscosity

# class to handle single aerosails
class AeroSail():
    def __init__(self, plainfoil, chord, chordratio, height, oswalde = 1, panels=160):
        self.wing = Sail.Sail(plainfoil, chord, chordratio, height, oswalde, panels)
        self.yaw = 0
        self.flapdeflection = 0.0
        self.chord = chord
        self.height = height
        self.oswalde = oswalde
        self.panels = panels
        self.chordratio = chordratio
        self.plainfoil = plainfoil
    def set_yaw(self, yaw):
        self.yaw = yaw
    def set_flapdeflection(self, flapdeflection):
        self.flapdeflection = flapdeflection
        self.wing.add_flap(self.flapdeflection)
    def set_apparent_wind(self, AWS, bearing, density=1.225, reynolds=None, viscosity=1.802*(10**(-5))):
        self.aws = AWS
        self.rho = density
        self.bearing = bearing
        self.viscosity = viscosity
        if reynolds is None:
            self.reynolds = (self.aws * self.rho * self.chord) / self.viscosity
        else:
            self.reynolds = reynolds
        self.wing.set_p('re', self.reynolds)
    # def get_t_sf_coefficients(self):


# class Ship():