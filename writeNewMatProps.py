import numpy as np

def writeNewMatProps(filename, callusElementNumbers, materials):

    matLine = 0

    # Extract the first column of materials
    updatedEAll = materials[:, 0]
    updatedNuAll = materials[:, 1]

    with open(filename + '.dat', 'r') as f:
        searchlines = f.readlines()

        for i, line in enumerate(searchlines):
        # find all materials in the .dat file by "isotropic"
            if "isotropic" in line: 
                if (searchlines[i+2].split()[4][1:7] == 'callus') and (searchlines[i+2].split()[4][7] != r'_'):
                    matLine = i+3
                    break

        for i, j in enumerate(callusElementNumbers):

            updatedE = '{:.15e}'.format(updatedEAll[i])
            updatedE = updatedE[:-4] + updatedE[-3] + updatedE[-1]
            updatedNu = '{:.15e}'.format(updatedNuAll[i])
            updatedNu = updatedNu[:-4] + updatedNu[-3] + updatedNu[-1]

            editLine = searchlines[matLine].split()
            editLine[0] = updatedE
            editLine[1] = updatedNu
            editLine = ' ' + ' '.join(editLine) + '\n'

            searchlines[matLine] = editLine

            matLine = matLine + 11

    with open(filename + '.dat', 'w') as f:
        updatedFile = ''.join(searchlines)
        f.write(updatedFile)
