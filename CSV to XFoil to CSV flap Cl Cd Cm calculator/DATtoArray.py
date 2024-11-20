import numpy as np

def DATtoArray(filename):
    input = open(filename, 'r')
    x = []
    y = []
    for line in input:
        xpart = ""
        ypart = ""
        switch = False
        for character in line:
            if (not switch) and character == ' ':
                switch = True
            elif character.isdigit() or character == '.' or character == '-':
                if not switch:
                    xpart += character
                else:
                    ypart += character
        print(xpart +' ' + ypart)
        x.append(float(xpart))
        y.append(float(ypart))
    x = (np.array(x))
    y = (np.array(y))
    input.close()
    return np.column_stack((x, y))

print(DATtoArray("E473coordinates.txt"))