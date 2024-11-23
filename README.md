# AeroSail-Simulator
Focused on the development of a simulator that simplifies and simulates the behaviour of an AeroSail system

## Profile
The profile class allows for the creation of custom airfoils using a seed airfoil and adding flaps to it. It uses the PyXfoil (https://github.com/Xero64/pyxfoil) repository to compute aerdynamic coefficient values for these airfoils. It also manages errors in convergance by the use of interpolated results improving stability

## Sail
The sail class allows for computation of finite wing profiles with flap simulation based on the profile class. It accounts for induced drag losses and induced alpha losses in cl losses
