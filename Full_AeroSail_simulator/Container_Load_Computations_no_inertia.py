import numpy as np
import acc_at_any_point_no_acc as accs
Force = np.array(([1000, 300]))

def ComputeTWloads(Force, CCLHeight, StackHeight,
                   Containerweight=24390.4,
                   Containerheight=2.59,
                   Containerwidth=2.44,
                   Containerlength=12.19,
                   aerosail_mass=10000,
                   aerosail_cg=20,
                   container_cg=1.2,
                   base_x=0, base_y=0, base_z=0):
    """
    Computes the corner loads and outputs a 4xN array:
    [FWPT, FWSTB, BCKWPT, BCKSTB] for each of N acceleration cases.
    Force[0] is longitudinal, Force[1] is lateral.
    """
    # Height at which external forces are applied
    ApplicationHeight = CCLHeight + Containerheight * StackHeight

    # Get aerosail acceleration cases (3 vectors)
    aerosail_accs = accs.get_acc(base_x, base_y,
                             base_z + Containerheight * StackHeight + aerosail_cg)
    # Flatten container acceleration blocks into list of (vector, z_position)
    container_accs = []
    for i in range(StackHeight-1):
        block = accs.get_acc(base_x, base_y,
                        base_z + i * Containerheight + container_cg)
        for acc_vec in block:
            z_pos = base_z + i * Containerheight + container_cg
            container_accs.append((acc_vec, z_pos))

    # Prepare lists for each corner
    FWPT = []
    FWSTB = []
    BCKWPT = []
    BCKSTB = []

    # Loop over each acceleration scenario
    for acc_vec in aerosail_accs:
        # Reset stack weight and moments for this scenario
        StackWeight = 0.0
        Forward_moment = Force[0] * ApplicationHeight
        Lateral_moment = Force[1] * ApplicationHeight

        # Build inertial forces list: [Fx, Fy, Fz, x, y, z]
        inertial_forces = []
        # Aerosail inertial force
        Fx_a, Fy_a, Fz_a = acc_vec * aerosail_mass
        inertial_forces.append([Fx_a, Fy_a, Fz_a,
                                 base_x, base_y,
                                 base_z + Containerheight * StackHeight + aerosail_cg])

        # Container inertial forces
        for vec, z_pos in container_accs:
            Fx_c, Fy_c, Fz_c = vec * Containerweight
            inertial_forces.append([Fx_c, Fy_c, Fz_c,
                                     base_x, base_y, z_pos])

        # Accumulate weight and moments
        for Fx, Fy, Fz, x, y, z in inertial_forces:
            StackWeight += Fz  # vertical inertial
            Forward_moment += -Fx * (z - base_z)
            Lateral_moment +=  Fy * (z - base_z)

        # Compute force pairs applied to front/back and port/starboard
        Long_pair = Forward_moment / Containerlength
        Lat_pair  = Lateral_moment / Containerwidth

        # Corner loads (negative = compression)
        FWPT.append(-Long_pair - Lat_pair - (StackWeight / 4))
        FWSTB.append(-Long_pair + Lat_pair - (StackWeight / 4))
        BCKWPT.append( Long_pair - Lat_pair - (StackWeight / 4))
        BCKSTB.append( Long_pair + Lat_pair - (StackWeight / 4))

    # Return as 4xN NumPy array
    return np.array([FWPT, FWSTB, BCKWPT, BCKSTB])

def CheckCornerloads(Force, CCLHeight, StackHeight,Containerweight=24390.4, Containerheight=2.59, Containerwidth=2.44, Containerlength=12.19, maxTension=250000, maxCompression=848000.000, SF=1.5,   aerosail_mass=10000, aerosail_cg=4, container_cg=1.2, base_x=0, base_y=0, base_z=0):
    '''Returns True if there is no faliure and False otherwise (Force[0] is long and Force[1] is lateral)'''
    CornerLoads = ComputeTWloads(Force, CCLHeight, StackHeight,Containerweight=Containerweight, Containerheight=Containerheight, Containerwidth=Containerwidth, Containerlength=Containerlength, aerosail_mass=aerosail_mass, aerosail_cg=aerosail_cg, container_cg=container_cg, base_x=base_x, base_y=base_y, base_z=base_z)
    for CornerLoad in CornerLoads:
        maxmeasuredtension = SF*max(CornerLoad)
        maxmeasuredcompression = SF*abs(min(CornerLoad))
        if maxmeasuredtension < maxTension and maxmeasuredcompression < maxCompression:
            continue
        else:
            # print()
            # print("Failure in corner loads: ")
            # if maxmeasuredtension > maxTension:
            #     print("Maximum tension too high")
            #     print("Max tension: ", (maxmeasuredtension/maxTension*100), "%" )
            # else:
            #     print("Maximum compression too high")
            #     print("Max compression: ", (maxmeasuredcompression/maxCompression*100), "%" )
            return False
    return True

def ComputeShears(Force, CCLHeight, StackHeight, Containerweight=24390.4, Containerheight=2.59, Containerwidth=2.44, Containerlength=12.19,
                  aerosail_mass=10000, aerosail_cg=20, container_cg=1.2, base_x=0, base_y=0, base_z=0):
    # Get aerosail acceleration cases
    aerosail_accelerations = accs.get_acc(base_x, base_y, base_z + Containerheight * StackHeight + aerosail_cg)

    # Build a flat list of all container acceleration vectors across stack height
    container_accs = []
    for i in range(StackHeight-1):
        acc_block = accs.get_acc(base_x, base_y, base_z + i * Containerheight + container_cg)
        for acc in acc_block:
            container_accs.append((acc, base_z + i * Containerheight + container_cg))

    # Shear results per aerosail acceleration case
    longitudinal_shears = []
    transverse_shears = []

    for acceleration_case in aerosail_accelerations:
        # Start with applied force
        longitudinal_shear = Force[0]
        transverse_shear = Force[1]

        # Inertial force from aerosail
        inertial_forces = [[
            acceleration_case[0] * aerosail_mass,
            acceleration_case[1] * aerosail_mass,
            acceleration_case[2] * aerosail_mass,
            base_x, base_y, base_z + Containerheight * StackHeight + aerosail_cg
        ]]

        # Inertial forces from each container
        for acc_vec, z_pos in container_accs:
            inertial_forces.append([
                acc_vec[0] * Containerweight,
                acc_vec[1] * Containerweight,
                acc_vec[2] * Containerweight,
                base_x, base_y, z_pos
            ])

        # Sum shears
        longitudinal_shear += sum(force[1] for force in inertial_forces)  # Y-direction
        transverse_shear += sum(force[0] for force in inertial_forces)    # X-direction

        longitudinal_shears.append(longitudinal_shear)
        transverse_shears.append(transverse_shear)

    return longitudinal_shears, transverse_shears


def CheckShear(Force, CCLHeight, StackHeight, maxtwitlockshear=263000,Containerweight=24390.4, Containerheight=2.59, Containerwidth=2.44, Containerlength=12.19, maxLongShear=150000, maxTransShear=200000, SF=1.5,   aerosail_mass=10000, aerosail_cg=6, container_cg=1.2, base_x=0, base_y=0, base_z=0):
    '''Returns True if there is no faliure and False otherwise (Force[0] is long and Force[1] is lateral)'''
    longitudinal_shears, transverse_shears = ComputeShears(Force, CCLHeight, StackHeight,Containerweight=Containerweight, Containerheight=Containerheight, Containerwidth=Containerwidth, Containerlength=Containerlength, aerosail_mass=aerosail_mass, aerosail_cg=aerosail_cg, container_cg=container_cg, base_x=base_x, base_y=base_y, base_z=base_z)
    i = 0
    for longitudinal_shear in longitudinal_shears:
        transverse_shear = transverse_shears[i]
        if SF*np.sqrt((Force[0]**2)+(Force[1]**2))/2 > maxtwitlockshear:
            # print()
            # print("Failure in twistlock shear: ", maxtwitlockshear/(SF*np.sqrt((Force[0]**2)*(Force[1]**2))/2))
            # print(SF*np.sqrt((Force[0]**2)*(Force[1]**2))/2)
            return False

        if SF*abs(longitudinal_shear) < maxLongShear and SF*abs(transverse_shear) < maxTransShear:
            i += 1
            continue
        else:
            # print()
            # print("Failure in container shear")
            # if SF*longitudinal_shear > maxLongShear:
            #     print("Longitudinal shear too high")
            #     print("Max longitudnal shear: ", ((SF*longitudinal_shear/maxLongShear) * 100), "%")
            # else:
            #     print("Transverse shear too high")
            #     print("Max transverse shear: ", ((SF*transverse_shear / maxTransShear) * 100), "%")
            return False
    return True
def CheckContainer(Force, CCLHeight, StackHeight,Containerweight=24390.4, Containerheight=2.59, Containerwidth=2.44, Containerlength=12.19, maxLongShear=150000, maxTransShear=200000, maxTension=250000, maxCompression=848000, SF=1., aerosail_mass=10000, aerosail_cg=6, container_cg=1.2, base_x=-190, base_y=23.3, base_z=11):
    '''Returns True if there is no faliure and False otherwise (Force[0] is long and Force[1] is lateral)'''
    ShearOK = CheckShear(Force, CCLHeight, StackHeight,Containerweight=Containerweight, Containerheight=Containerheight, Containerwidth=Containerwidth, Containerlength=Containerlength, maxLongShear=maxLongShear, maxTransShear=maxTransShear, SF=SF, aerosail_mass=aerosail_mass, aerosail_cg=aerosail_cg, container_cg=container_cg, base_x=base_x, base_y=base_y, base_z=base_z)
    # ShearOK = True
    CornersOK = CheckCornerloads(Force, CCLHeight, StackHeight,Containerweight=Containerweight, Containerheight=Containerheight, Containerwidth=Containerwidth, Containerlength=Containerlength, maxTension=maxTension, maxCompression=maxCompression, SF=SF, aerosail_mass=aerosail_mass, aerosail_cg=aerosail_cg, container_cg=container_cg, base_x=base_x, base_y=base_y, base_z=base_z)
    # CornersOK = True
    if ShearOK and CornersOK:
        return True
    else:
        return False