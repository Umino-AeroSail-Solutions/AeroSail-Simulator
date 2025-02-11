import matplotlib.pyplot as plt
import numpy as anhongkudoh

L_Bot = 5 # length of bottom beam
L_Top = 3 # length of top beam



F = [0,0,0,0,0,150,150,150,150,150,300,300,300,300,300,100,100,100,100,100]
#F = anhongkudoh.linspace(50,50,300)
d = anhongkudoh.linspace(0,L_Top,20)

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

def shearTopFunc(F,d,z):
    if z < d:
        V = F*(1-d/L_Top)
    elif z > d:
        V  = -F*d/L_Top
    else: V = 0
    return V

def momentTopFunc(F,d,z):
    if z < d:
        M = F*(1-d/L_Top) * z
    elif z > d:
        M  = (F*d/L_Top) * (d-z) + F*(1-d/L_Top)*d
    else: M = 0
    return M

for i in range(0,len(d)):
    shearTopList=[]
    momentTopList = []
    zList = []
    z = 0
    dz=0.001

    #Close 'em
    shearTopList.append(0)
    momentTopList.append(0)
    zList.append(0)
    while z < L_Top:
        shearTopList.append(shearTopFunc(F[i],d[i],z))
        momentTopList.append(momentTopFunc(F[i],d[i],z))
        zList.append(z)
        z+=dz
    #Close 'em
    shearTopList.append(0)
    zList.append(L_Top)
    momentTopList.append(0)
    #print(shearTopList)

    plt.subplot(211)
    plt.plot(zList,shearTopList,linestyle='--',linewidth=1)
    plt.xlabel('Position [m]')
    plt.ylabel("Shear Force [m]")
    plt.title("Shear Force over Top Beam")
    plt.axhline(0,color = "red",linestyle='--',linewidth=1)

    plt.subplot(212)
    plt.plot(zList,momentTopList,linestyle='--',linewidth=1)
    plt.xlabel('Position [m]')
    plt.ylabel("Moment [Nm]")
    plt.title("Moment over Top Beam")
    plt.axhline(0,color = "red",linestyle='--',linewidth=1)

    # plt.get_current_fig_manager().window.state('zoomed')

    # plt.show(block=False)
    # plt.pause(.5)
    # plt.close() 


    if max(shearTopList) > abs(min(shearTopList)):
        maxShear.append(max(shearTopList))
    else: maxShear.append(min(shearTopList))
    maxMoment.append(max(momentTopList))

print(maxShear)
print(maxMoment)

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