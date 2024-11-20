import numpy as np

def getExcelAlphaDelta(inputfile):
    input = open(inputfile, "r")
    alphas = []
    deltas = []
    for line in input:
        if line != ' ':
            started = False
            deltabool = False
            alpha = ""
            delta = ""
            for character in line:
                if character == ',':
                    character = '.'
                if character != ' ' and (not started):
                    started = True
                if character == '\t':
                    deltabool = True
                if started and character != ' ' and character != '\n' and character != '\t':
                    if not deltabool:
                        alpha += character
                    else:
                        delta += character
        print(alpha, delta)
        alphas.append(float(alpha))
        deltas.append(float(delta))
    input.close()
    return np.array(alphas), np.array(deltas)

a, d = getExcelAlphaDelta("input.txt")