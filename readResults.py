import numpy as np
import numpy.linalg as la
from py_post import *

def readResults(newFileName, callusElementNumbers, matProps, nodeID, IFMDoF_String):

    # create a postpy obj using post_open
    p = post_open(newFileName + ".t16")

    # try opening the results file to check for errors
    try:
        p.moveto(1)
    except:
        print("Error opening post file: " + newFileName + ".t16")

    # move to the last increment
    numIncrements = p.increments() # number of increments in the sim
    p.extrapolation("linear") # set extrapolation to linear
    p.moveto(numIncrements - 1) # move to the final increment

    nns = p.node_scalars() # the number of nodal scalars in results
    nes = p.element_scalars() # the number of element scalars in results

    # read the displacement of the free node (IFM)
    for i in range(0, nns):
        if(p.node_scalar_label(i) == IFMDoF_String):
            IFM = p.node_scalar(nodeID - 1, i)

    minStrainsScalarIdx = 0
    intStrainsScalarIdx = 0
    maxStrainsScalarIdx = 0

    for i in range(0, nes):
        if(p.element_scalar_label(i) == "Minimum Principal Total Strain"):
            minStrainsScalarIdx = i
        elif(p.element_scalar_label(i) == "Intermediate Principal Total Strain"):
            intStrainsScalarIdx = i
        elif(p.element_scalar_label(i) == "Maximum Principal Total Strain"):
            maxStrainsScalarIdx = i

    numElements = len(callusElementNumbers)
    minStrains = np.zeros(numElements)
    intStrains = np.zeros(numElements)
    maxStrains = np.zeros(numElements)

    for i, j in enumerate(callusElementNumbers):
        minStrains[i] = p.element_scalar(j - 1, minStrainsScalarIdx)[0].value
        intStrains[i] = p.element_scalar(j - 1, intStrainsScalarIdx)[0].value
        maxStrains[i] = p.element_scalar(j - 1, maxStrainsScalarIdx)[0].value

    minStrains[minStrains <= -0.5] = -0.4999
    minStrains = (2 * minStrains + 1)**0.5 - 1
    intStrains = (2 * intStrains + 1)**0.5 - 1
    maxStrains = (2 * maxStrains + 1)**0.5 - 1

    hydrostaticStrains = (minStrains + intStrains + maxStrains) / 3
    equivalentStrains = (0.5*((minStrains - intStrains)**2 + (intStrains - maxStrains)**2 + (maxStrains - minStrains)**2))**0.5

    p.close()

    return hydrostaticStrains, equivalentStrains, IFM, minStrains, maxStrains



