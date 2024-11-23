import numpy as np
import matplotlib.pyplot as plt

import Profile
import Sail

Profile.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')

#Class to handle environmental conditions
class Conditions():
    def __init__(self, speed, bearing, rho=1.225,re_per_chord=None, viscosity=1.802*(10**(-5))):
        self.speed = speed
        self.bearing = bearing
        self.rho = rho
        self.viscosity = viscosity
        if re_per_chord is None:
            self.re_per_chord = (self.speed * self.rho) / self.viscosity

# class to handle single aerosails
class AeroSail():
    def __init__(self, plainfoil, chord, chordratio, height=None, oswalde = 1, panels=160):
        self.wing = Sail.Sail(plainfoil, chord, chordratio, height, oswalde, panels)
        self.yaw = 0
        self.flapdeflection = 0.0
    def set_yaw(self, yaw):
        self.yaw = yaw
    def set_flapdeflection(self, flapdeflection):
        self.flapdeflection = flapdeflection
        self.wing.add_flap(self.flapdeflection)

# class Ship():