import numpy as np
import numpy.linalg as la
from py_post import *

def calcVolumes(fileName, callusElementNumbers):

    #referenceVolume = 0.5

    # create a postpy obj using post_open
    p = post_open(fileName + ".t16")

    # try opening the results file to check for errors
    try:
        p.moveto(1)
    except:
        print("Error opening post file: " + fileName + ".t16")

    # move to the last increment
    numIncrements = p.increments() # number of increments in the sim
    p.extrapolation("linear") # set extrapolation to linear
    p.moveto(numIncrements - 1) # move to the final increment

    volumesRelative = [] # contains the volume of each element, divided by the mean volume, the elements are ordered in the same way as callusElementNumbers

    for j in callusElementNumbers: # callusElementNumbers IDs, not indexing...

        nodes = p.element(j-1).items

        points = np.zeros(shape=(4, 4))

        for i in range(4):
            points[0, i] = p.node(nodes[i]-1).x # -1 because INDEXING starts at 0 even though IDs start at 1...
            points[1, i] = p.node(nodes[i]-1).y
            points[2, i] = p.node(nodes[i]-1).z
            points[3, i] = 1

        elVolume = abs((1/6)*la.det(np.array([points[:,0], points[:,1], points[:,2], points[:,3]]))) # volume of a tetrahedral is equal to 1/6th of the determinant of the matrix formed by the point coordinates

        volumesRelative.append(elVolume)


    volumesRelative = volumesRelative/np.mean(volumesRelative)
    #volumesRelative = volumesRelative/np.mean(referenceVolume)

    p.close()

    return volumesRelative

def calcAreas(fileName, callusElementNumbers):

    # create a postpy obj using post_open
    p = post_open(fileName + ".t16")

    # try opening the results file to check for errors
    try:
        p.moveto(1)
    except:
        print("Error opening post file: " + fileName + ".t16")

    # move to the last increment
    numIncrements = p.increments() # number of increments in the sim
    p.extrapolation("linear") # set extrapolation to linear
    p.moveto(numIncrements - 1) # move to the final increment

    areasRelative = [] # contains the volume of each element, divided by the mean volume, the elements are ordered in the same way as callusElementNumbers

    for j in callusElementNumbers: # callusElementNumbers IDs, not indexing...

        nodes = p.element(j-1).items

        points = np.zeros(shape=(4, 3))

        for i in range(3):
            points[0, i] = p.node(nodes[i]-1).x # -1 because INDEXING starts at 0 even though IDs start at 1...
            points[1, i] = p.node(nodes[i]-1).y
            points[2, i] = 1 # just set 1 for all z coords
            points[3, i] = 1

        #elArea = abs((1/6)*la.det(np.array([points[:,0], points[:,1], points[:,2], points[:,3]]))) 
        elArea = abs(0.5*la.det(points[:2, :2].T))

        areasRelative.append(elArea)


    areasRelative = areasRelative/np.mean(areasRelative)
    #volumesRelative = volumesRelative/np.mean(referenceVolume)

    p.close()

    return areasRelative