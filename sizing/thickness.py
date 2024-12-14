import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt

#Mast heights
sailHeight = 30 #Mast height [m]
contHeight = 2.896 #Container height [m]
mastLength = sailHeight+contHeight #Length of the mast [m]
nRibs = 8 #Number of ribs [m]

#Mast cross section
c = 5 #Maximum assumed chord length [m] (i made this up)
w = 0.14 * c #Width of the mast cross-section ASSUMING 14% OF CHORD LENGTH (NOT FINAL) 
t = 0.05

Ixx = t*w**3/6 + w**2 * t #Area moment of inertia (for square cross-section)
Qmax = w**2 *t * 3/4 #Maximum First moment of area (for square cross-section)

#Force distribution at ribs
totForce = 2000 #Total aerodynamic forcee [N]


fList = []
hList = []
dh = 0.1
h = 0
while h < mastLength+dh:
    if h < contHeight:
        fList.append(totForce)
    elif h < contHeight + sailHeight / (nRibs-1):
        fList.append(totForce * (nRibs-1)/nRibs)
    elif h < contHeight + 2*sailHeight / (nRibs-1):
        fList.append(totForce * (nRibs-2)/nRibs)
    elif h < contHeight + 3*sailHeight / (nRibs-1):
        fList.append(totForce * (nRibs-3)/nRibs)
    elif h < contHeight + 4*sailHeight / (nRibs-1):
        fList.append(totForce * (nRibs-4)/nRibs)
    elif h < contHeight + 5*sailHeight / (nRibs-1):
        fList.append(totForce * (nRibs-5)/nRibs)
    elif h < contHeight + 6*sailHeight / (nRibs-1):
        fList.append(totForce * (nRibs-6)/nRibs)
    elif h < contHeight + 7*sailHeight / (nRibs-1):
        fList.append(totForce * (nRibs-7)/nRibs)
    elif h > mastLength:
        fList.append(0)
    hList.append(h)
    h+=dh

fV = sp.interpolate.interp1d(hList,fList,kind='linear',fill_value="interpolate")
def fBM(x): #Bending moment
    return -sp.integrate.quad(fV,x,mastLength)[0]
bmList = []
for i in hList:
    bmList.append(fBM(i))


shearMax = fV(0)
momentMax = fBM(0)

normStressMax = momentMax * (w/2) / Ixx
shearStressMax = shearMax * Qmax / Ixx / t

print(f"Max shear: {shearMax}")
print(f"Max moment: {momentMax}")
print(f"Max normal stress: {normStressMax}")
print(f"Max shear stress: {shearStressMax}")

plt.subplot(211)
plt.plot(hList, fV(hList))
plt.xlabel('Mast Length [m]')
plt.ylabel("Force [N]")
plt.title("Force distribution")

plt.subplot(212)
plt.plot(hList, bmList)
plt.xlabel('Mast Length [m]')
plt.ylabel("Bending moment [Nm]")
plt.title("Bending moment distribution")

plt.show()
