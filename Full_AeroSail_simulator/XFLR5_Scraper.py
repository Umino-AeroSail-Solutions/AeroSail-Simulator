import numpy as np


def read_polar_data(file_path):
    # Initialize lists to store data
    alphas = []
    CL = []
    CDi = []
    CDv = []
    CD = []
    CY = []
    Cl = []
    Cm = []
    Cn = []
    Cni = []
    QInf = []
    XCP = []

    with open(file_path, 'r') as file:
        # Skip header lines until we find numeric data
        for line in file:
            if line.strip() and line.strip().split()[0].replace('.', '', 1).replace('-', '', 1).isdigit():
                data = line.split()
                alphas.append(float(data[0]))
                CL.append(float(data[2]))
                CDi.append(float(data[3]))
                CDv.append(float(data[4]))
                CD.append(float(data[5]))
                CY.append(float(data[6]))
                Cl.append(float(data[7]))
                Cm.append(float(data[8]))
                Cn.append(float(data[9]))
                Cni.append(float(data[10]))
                QInf.append(float(data[11]))
                XCP.append(float(data[12]))

    # Convert lists to NumPy arrays
    alphas = np.array(alphas)
    CL = np.array(CL)
    CDi = np.array(CDi)
    CDv = np.array(CDv)
    CD = np.array(CD)
    CY = np.array(CY)
    Cl = np.array(Cl)
    Cm = np.array(Cm)
    Cn = np.array(Cn)
    Cni = np.array(Cni)
    QInf = np.array(QInf)
    XCP = np.array(XCP)

    return [alphas, CL, CDi, CDv, CD, CY, Cl, Cm, Cn, Cni, QInf, XCP]


# Example usage:
file_path = 'Data/T1-10_0_m_s-LLT-Flap-10.txt'
alphas = read_polar_data(file_path)[0]
CL = read_polar_data(file_path)[1]

# Print some example data
print('Alphas:', alphas)
print('CL:', CL)
