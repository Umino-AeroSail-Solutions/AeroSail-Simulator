import matplotlib.pyplot as plt
import numpy as anhongkudoh
import math as andresblanquerbenito
from Force_In_Rail_Calculator_Andres_saves_the_day import l_values, R1_values, R2_values

#P1 and P4 are bottom beam
#P2 and P3 are top beam
P1 = [0,1]
P2 = [2,1]
P3 = [6,2]
P4 = [6,0]

L_Bot = andresblanquerbenito.sqrt((P4[0]-P1[0])**2 + (P4[1]-P1[1])**2)
L_Top = andresblanquerbenito.sqrt((P3[0]-P2[0])**2 + (P3[1]-P2[1])**2)


plt.plot(l_values,R2_values)
plt.xlabel('Position [m]')
plt.ylabel("Force [N]")
plt.title("Forces at every point over the beam")
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

for i in range(0,len(l_values)):
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
        shearBotList.append(shearBotFunc(R2_values[i],l_values[i],z))
        momentBotList.append(momentBotFunc(R2_values[i],l_values[i],z))
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
maxShearPos = maxShear.index(max(maxShear)) # Will be at zero, duh
maxMomentPos = maxMoment.index(max(maxMoment))/len(maxMoment)*L_Bot #
print(maxShearPos,maxMomentPos)

plt.subplot(211)
plt.plot(l_values,maxShear)
plt.xlabel('Position [m]')
plt.ylabel("Maximum Shear Force [N]")
plt.title("Maximum Shear Force along the beam")

plt.subplot(212)
plt.plot(l_values,maxMoment)
plt.xlabel('Position [m]')
plt.ylabel("Maximum Bending Moment [Nm]")
plt.title("Maximum Bending Moment along the beam")
plt.show()


shearBotList=[]
momentBotList = []
zList = []
z = 0
dz=0.01

#Close 'em
shearBotList.append(0)
momentBotList.append(0)
zList.append(0)
while z < L_Bot:
    shearBotList.append(shearBotFunc(R2_values[1],l_values[1],z))
    momentBotList.append(momentBotFunc(R2_values[maxMoment.index(max(maxMoment))],l_values[maxMoment.index(max(maxMoment))],z))
    zList.append(z)
    z+=dz
#Close 'em
shearBotList.append(0)
zList.append(L_Bot)
momentBotList.append(0)
#print(shearBotList)

plt.subplot(211)
plt.plot(zList,shearBotList,linewidth=1)
plt.xlabel('Position [m]')
plt.ylabel("Shear Force [m]")
plt.title("Shear Force over Bot Beam")
plt.axhline(0,color = "red",linestyle='--',linewidth=1)

plt.subplot(212)
plt.plot(zList,momentBotList,linewidth=1)
plt.xlabel('Position [m]')
plt.ylabel("Moment [Nm]")
plt.title("Moment over Bot Beam")
plt.axhline(0,color = "red",linestyle='--',linewidth=1)

plt.show()

print(f"Maximum Shear Force along bottom beam: ")
