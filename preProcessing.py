import numpy as np
from shutil import copyfile
from scipy.sparse import csr_matrix, identity

#def preProcessing(fileName, newFileName, rigidityFileName, rigidityNewFileName):
def preProcessing(fileName, newFileName, numNodesNeighbor):
    
    # create copy file to preserve original .dat files
    copyfile(fileName + '.dat', newFileName + '.dat')
    #copyfile(rigidityFileName + '.dat', rigidityNewFileName + '.dat')

    numberNodesForANeighbor = numNodesNeighbor

    callusElements = [] # list of all callus element IDs
    callusSurfaceElements = []
    callusMedullarySurfaceElements = []
    callusBoneSurfaceElements = []
    callusConnectivityNodes = [] # 2D array, elementID and node IDs in each row
    callusConnectivityElements = {} # dictionary with element IDs as keys, and a list of neighboring element IDs as the values
    callusConnectivityNodesDict = {} # dictionary with element IDs as keys, and a set of nodes as the values

    # opens the new .dat file and reads in materials data and the element numbers for the callus (the set, not material... can't find element numbers for the material)
    with open(newFileName + '.dat', 'r') as f:
        searchlines = f.readlines()
    for i, line in enumerate(searchlines):

        if "define              element             set                 callus\n" in line:
            callusElementsLowerBound = int(searchlines[i+1].split()[0])
            callusElementsUpperBound = int(searchlines[i+1].split()[2]) 
            callusElements = list(range(callusElementsLowerBound,callusElementsUpperBound + 1))

        #if "define              element             set                 callusSurfaceElements\n" in line:
        #    j = 1
        #    while searchlines[i+j].split()[-1] == "c": # while there is another line (current line ends in "c")
        #        for k in range(len(searchlines[i+j].split()) - 1):
        #            callusSurfaceElements.append(int(searchlines[i+j].split()[k])) # append all the element numbers in the line (but need to leave out the "c"
        #        j += 1
        #    for k in range(len(searchlines[i+j].split())): # then do the last line (doesn't end in "c", even if a full line
        #        callusSurfaceElements.append(int(searchlines[i+j].split()[k]))

        #if "define              element             set                 callusMedullarySurfaceElements\n" in line:
        #    j = 1
        #    while searchlines[i+j].split()[-1] == "c": # while there is another line (current line ends in "c")
        #        for k in range(len(searchlines[i+j].split()) - 1):
        #            callusMedullarySurfaceElements.append(int(searchlines[i+j].split()[k])) # append all the element numbers in the line (but need to leave out the "c"
        #        j += 1
        #    for k in range(len(searchlines[i+j].split())): # then do the last line (doesn't end in "c", even if a full line
        #        callusMedullarySurfaceElements.append(int(searchlines[i+j].split()[k]))

        if "define              element             set                 callusBoneSurfaceElements\n" in line:
            j = 1
            while searchlines[i+j].split()[-1] == "c": # while there is another line (current line ends in "c")
                for k in range(len(searchlines[i+j].split()) - 1):
                    callusBoneSurfaceElements.append(int(searchlines[i+j].split()[k])) # append all the element numbers in the line (but need to leave out the "c"
                j += 1
            for k in range(len(searchlines[i+j].split())): # then do the last line (doesn't end in "c", even if a full line
                callusBoneSurfaceElements.append(int(searchlines[i+j].split()[k]))

    # find where the callus element connectivity starts in the .dat file
    callusElementConnectivityStartLine = 0
    for i, line in enumerate(searchlines):
        if "connectivity" in line:
            if int(searchlines[i+2].split()[0]) == callusElements[0]:
                callusElementConnectivityStartLine = i+2

    # starting from the first callus element, build callusConnectivityNodesDict each element at a time
    for i, elementID in enumerate(callusElements):
        currentLine = [int(x) for x in searchlines[(i*3) + callusElementConnectivityStartLine].split()]
        callusConnectivityNodesDict[currentLine[0]] = set(currentLine[2:])

    # Find the maximum element id and the maximum node ID across all sets
    max_node_id = max(max(node_set, default=0) for node_set in callusConnectivityNodesDict.values())
    max_element_id = max(callusConnectivityNodesDict.keys())

    # create a binary sparse matrix with rows as element IDs and columns as nodes, for each row, put a 1 in each column corresponding to the nodes present in the element
    # this essentially turns callusConnectivityNodesDict into a binary sparse matrix

    row_indices, col_indices = [], []

    for element, nodes in callusConnectivityNodesDict.items():
        row_indices.extend([element - 1] * len(nodes))
        col_indices.extend([node - 1 for node in nodes])

    data = np.ones(len(row_indices), dtype=int)
    connectivityMatrix = csr_matrix((data, (row_indices, col_indices)), shape=(max_element_id, max_node_id), dtype=int)

    # matrix multiplication of the sparse matrix with its transpose. Resultant matrix will be size n x n, where n is the number of callus element IDs (approximately - in reality it will be the size of the highest ID which is not the same, since the callus element IDs do not start from 1)
    # each cell of the matrix corresponds to a pair of elements, and the number in the cell equals the number of shared nodes between the element pair. This matrix is also sparse.
    # this is the fun step

    sharedNodes = connectivityMatrix @ connectivityMatrix.transpose()

    # set a binary mask to the sharedNodes matrix, placing a 1 where the element pair has a number of shared nodes greater or equal to the required number of nodes to define a neighbor, and zeros elsewhere
    # subtract the identity matrix, as each element will have added itself as a neighbor

    sharedNodesBinary = csr_matrix((sharedNodes >= numberNodesForANeighbor).astype(int)) - identity(sharedNodes.shape[0])

    # turn the resulting binary matrix into a dictionary (callusConnectivityElements)

    rows, cols = sharedNodesBinary.nonzero()

    for row, col in zip(rows, cols):
        callusConnectivityElements.setdefault(row, []).append(col)

    # make a set of callusElements to speed up the following step
    callusElementsSet = set(callusElements)

    # add one to each key and value in the dictionary since they will have added indices, but we want to store the actual element IDs
    callusConnectivityElements = {key + 1: [val + 1 for val in values] for key, values in callusConnectivityElements.items() if key + 1 in callusElementsSet}

    return callusElements, callusSurfaceElements, callusMedullarySurfaceElements, callusBoneSurfaceElements, callusConnectivityElements, callusConnectivityNodesDict
