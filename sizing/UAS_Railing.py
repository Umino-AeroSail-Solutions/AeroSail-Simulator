import matplotlib.pyplot as plt
import numpy as anhongkudoh
import math as andresblanquerbenito
from Force_In_Rail_Calculator import F

#P1 and P4 are bottom beam
#P2 and P3 are top beam
P1 = [0,1]
P2 = [2,1]
P3 = [8,2]
P4 = [8,0]

L_Bot = andresblanquerbenito.sqrt((P4[0]-P1[0])**2 + (P4[1]-P1[1])**2)
L_Top = andresblanquerbenito.sqrt((P3[0]-P2[0])**2 + (P3[1]-P2[1])**2)

print(len(F))
print(L_Bot)
d = anhongkudoh.linspace(0,L_Bot,len(F))

plt.plot(d,F)
plt.xlabel('Position [m]')
plt.ylabel("Force [N]")
plt.title("Unit Forces Distributed over the beam")
plt.show()

maxShear = []
maxMoment = []

'''
R1     ^      R2
|======|======|           -> z
v      F      v          |
                         v y

|------|
   d
'''

def shearBotFunc(F,d,z):
    if z < d:
        V = F*(1-d/L_Bot)
    elif z > d:
        V  = -F*d/L_Bot
    else: V = 0
    return V

def momentBotFunc(F,d,z):
    if z < d:
        M = F*(1-d/L_Bot) * z
    elif z > d:
        M  = (F*d/L_Bot) * (d-z) + F*(1-d/L_Bot)*d
    else: M = 0
    return M

for i in range(0,len(d)):
    shearBotList=[]
    momentBotList = []
    zList = []
    z = 0
    dz=0.001

    #Close 'em
    shearBotList.append(0)
    momentBotList.append(0)
    zList.append(0)
    while z < L_Bot:
        shearBotList.append(shearBotFunc(F[i],d[i],z))
        momentBotList.append(momentBotFunc(F[i],d[i],z))
        zList.append(z)
        z+=dz
    #Close 'em
    shearBotList.append(0)
    zList.append(L_Bot)
    momentBotList.append(0)
    #print(shearBotList)

    plt.subplot(211)
    plt.plot(zList,shearBotList,linestyle='--',linewidth=1)
    plt.xlabel('Position [m]')
    plt.ylabel("Shear Force [m]")
    plt.title("Shear Force over Bot Beam")
    plt.axhline(0,color = "red",linestyle='--',linewidth=1)

    plt.subplot(212)
    plt.plot(zList,momentBotList,linestyle='--',linewidth=1)
    plt.xlabel('Position [m]')
    plt.ylabel("Moment [Nm]")
    plt.title("Moment over Bot Beam")
    plt.axhline(0,color = "red",linestyle='--',linewidth=1)

    # plt.get_current_fig_manager().window.state('zoomed')

    # plt.show(block=False)
    # plt.pause(.5)
    # plt.close() 

    # Get maximum / minimum
    if max(shearBotList) > abs(min(shearBotList)): 
        maxShear.append(max(shearBotList))
    else: maxShear.append(min(shearBotList))
    maxMoment.append(max(momentBotList))

#print(maxShear)
#print(maxMoment)

plt.subplot(211)
plt.plot(d,maxShear)
plt.xlabel('Position [m]')
plt.ylabel("Maximum Shear Force [N]")
plt.title("Maximum Shear Force along the beam")

plt.subplot(212)
plt.plot(d,maxMoment)
plt.xlabel('Position [m]')
plt.ylabel("Maximum Bending Moment [Nm]")
plt.title("Maximum Bending Moment along the beam")
plt.show()