import numpy as np
from py_post import *
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator

def readTorsion(rigidityNewFileName, nodeID):

    # create a postpy obj using post_open
    p = post_open(rigidityNewFileName + ".t16")

    # try opening the results file to check for errors
    try:
        p.moveto(1)
    except:
        print("Error opening post file: " + rigidityNewFileName + ".t16")

    # move to the last increment
    numIncrements = p.increments() # number of increments in the sim
    p.extrapolation("linear") # set extrapolation to linear
    p.moveto(numIncrements - 1) # move to the final increment

    nns = p.node_scalars() # the number of nodal scalars in results

    for i in range(0, nns):
        if(p.node_scalar_label(i) == "Reaction Moment Z"):
            reaction = p.node_scalar(nodeID - 1, i) # in Nmm

    p.close()

    return reaction