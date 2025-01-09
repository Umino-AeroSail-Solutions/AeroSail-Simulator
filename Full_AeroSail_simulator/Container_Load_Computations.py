import numpy as np

Force = np.array(([1000, 300]))

def ComputeTWloads(Force, CCLHeight, StackHeight,Containerweight=24390.4, Containerheight=2.59, Containerwidth=2.44, Containerlength=12.19):
    '''Computes the corner loads and outputs [FWPT, FWSTB, BCKWPT, BCKSTB] array (Force[0] is long and Force[1] is lateral)'''
    ApplicationHeight = CCLHeight+(Containerheight*StackHeight)
    StackWeight = Containerweight*StackHeight*9.81

    Forward_moment = Force[0] * ApplicationHeight
    Lateral_moment = Force[1] * ApplicationHeight

    Longitudinal_force_pair = Forward_moment/Containerlength
    Lateral_force_pair = Lateral_moment/Containerwidth

    Forward_port_twlock_force = -Longitudinal_force_pair-Lateral_force_pair-(StackWeight/4) #Negative is compression
    Forward_starboard_twlock_force = -Longitudinal_force_pair+Lateral_force_pair-(StackWeight/4)
    Back_port_twlock_force = +Longitudinal_force_pair - Lateral_force_pair - (StackWeight / 4)  # Negative is compression
    Back_starboard_twlock_force = +Longitudinal_force_pair + Lateral_force_pair - (StackWeight / 4)

    CornerLoads = np.array(([Forward_port_twlock_force, Forward_starboard_twlock_force, Back_port_twlock_force, Back_starboard_twlock_force]))

    return CornerLoads

def CheckCornerloads(Force, CCLHeight, StackHeight,Containerweight=24390.4, Containerheight=2.59, Containerwidth=2.44, Containerlength=12.19, maxTension=250000, maxCompression=848000.000, SF=1.5):
    '''Returns True if there is no faliure and False otherwise (Force[0] is long and Force[1] is lateral)'''
    CornerLoads = ComputeTWloads(Force, CCLHeight, StackHeight,Containerweight=Containerweight, Containerheight=Containerheight, Containerwidth=Containerwidth, Containerlength=Containerlength)
    maxmeasuredtension = SF*max(CornerLoads)
    maxmeasuredcompression = SF*abs(min(CornerLoads))
    if maxmeasuredtension < maxTension and maxmeasuredcompression < maxCompression:
        return True
    else:
        print()
        print("Failure in corner loads")
        if maxmeasuredtension > maxTension:
            print("Maximum tension too high")
        else:
            print("Maximum compression too high")
        return False

def ComputeShears(Force, CCLHeight, StackHeight,Containerweight=24390.4, Containerheight=2.59, Containerwidth=2.44, Containerlength=12.19):
    # CornerLoads = ComputeTWloads(Force, CCLHeight, StackHeight,Containerweight=Containerweight, Containerheight=Containerheight, Containerwidth=Containerwidth, Containerlength=Containerlength)
    # longitudinal_shear = CornerLoads[2]-CornerLoads[0]+Force[0] # Adding differential tension which I think is wrong
    # transverse_shear = CornerLoads[3]-CornerLoads[1]+Force[1] # Adding differential tension which I think is wrong
    longitudinal_shear = Force[0]
    transverse_shear = Force[1]
    return longitudinal_shear, transverse_shear

def CheckShear(Force, CCLHeight, StackHeight, maxtwitlockshear=263000,Containerweight=24390.4, Containerheight=2.59, Containerwidth=2.44, Containerlength=12.19, maxLongShear=150000, maxTransShear=200000, SF=1.5):
    '''Returns True if there is no faliure and False otherwise (Force[0] is long and Force[1] is lateral)'''
    longitudinal_shear, transverse_shear = ComputeShears(Force, CCLHeight, StackHeight,Containerweight=Containerweight, Containerheight=Containerheight, Containerwidth=Containerwidth, Containerlength=Containerlength)

    if SF*np.sqrt((Force[0]**2)+(Force[1]**2))/2 > maxtwitlockshear:
        print()
        print("Failure in twistlock shear: ", maxtwitlockshear/(SF*np.sqrt((Force[0]**2)*(Force[1]**2))/2))
        print(SF*np.sqrt((Force[0]**2)*(Force[1]**2))/2)
        return False

    if SF*abs(longitudinal_shear) < maxLongShear and SF*abs(transverse_shear) < maxTransShear:
        return True
    else:
        print()
        print("Failure in container shear")
        if longitudinal_shear > maxLongShear:
            print("Longitudinal shear too high")
        else:
            print("Transverse shear too high")
        return False
def CheckContainer(Force, CCLHeight, StackHeight,Containerweight=24390.4, Containerheight=2.59, Containerwidth=2.44, Containerlength=12.19, maxLongShear=150000, maxTransShear=200000, maxTension=250000, maxCompression=848000, SF=1.5):
    '''Returns True if there is no faliure and False otherwise (Force[0] is long and Force[1] is lateral)'''
    ShearOK = CheckShear(Force, CCLHeight, StackHeight,Containerweight=Containerweight, Containerheight=Containerheight, Containerwidth=Containerwidth, Containerlength=Containerlength, maxLongShear=maxLongShear, maxTransShear=maxTransShear, SF=SF)
    # ShearOK = True
    CornersOK = CheckCornerloads(Force, CCLHeight, StackHeight,Containerweight=Containerweight, Containerheight=Containerheight, Containerwidth=Containerwidth, Containerlength=Containerlength, maxTension=maxTension, maxCompression=maxCompression, SF=SF)
    # CornersOK = True
    if ShearOK and CornersOK:
        return True
    else:
        return False