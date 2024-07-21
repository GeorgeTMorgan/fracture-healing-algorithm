import numpy as np
import numpy.linalg as la
from py_post import *
import builtins

def makeNeighbors(filename, callusElementNumbers):

    centroidDist = 1 # define the distance between element centroids for neighboring elements

    callusConnectivityElements = {}
    centroids = []

    # create a postpy obj using post_open
    p = post_open(filename + ".t16")

    # try opening the results file to check for errors
    try:
        p.moveto(1)
    except:
        print("Error opening post file: " + filename + ".t16")

    # move to the last increment
    numIncrements = p.increments() # number of increments in the sim
    p.extrapolation("linear") # set extrapolation to linear
    p.moveto(numIncrements - 1) # move to the final increment

    for j in callusElementNumbers: # callusElementNumbers IDs, not indexing...

        nodes = p.element(j-1).items

        points = np.zeros(shape=(3, 4))

        for i in range(4):
            points[0, i] = p.node(nodes[i]-1).x # -1 because INDEXING starts at 0 even though IDs start at 1...
            points[1, i] = p.node(nodes[i]-1).y
            points[2, i] = p.node(nodes[i]-1).z

        centroid = np.mean(points, axis=1)
        centroids.append(centroid)

    def calcDist(centroid1, centroid2):
        dist = ((centroid1[0]-centroid2[0])**2 + (centroid1[1]-centroid2[1])**2 + (centroid1[2]-centroid2[2])**2)**0.5
        return dist

    for j in range(len(callusElementNumbers)):
        neighborElIDs = [] # list of neighbor elements

        for i in range(len(callusElementNumbers)):
            if calcDist(centroids[j], centroids[i]) <= centroidDist:
                neighborElIDs.append(callusElementNumbers[i])

        neighborElIDs.remove(callusElementNumbers[j]) # each element will have added itself as a neighbor
        callusConnectivityElements[callusElementNumbers[j]] = neighborElIDs # add the list of neighbors to the dictionary

    return callusConnectivityElements