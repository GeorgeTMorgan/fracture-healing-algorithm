import numpy as np
import numpy.linalg as la
from py_post import *

def calcProximity(fileName, callusElementNumbers):

    proximityWeightsIO = []
    proximityWeightsChon = []

    # Proximity function parameters
    a = 8 # gap width function 'speed'
    b = 1.07 # gap width / 2
    c = 1.5 # radial function 'speed'
    d = 8 # radial boundary

    u = c/2 # axial function decay 'speed'
    v = d # axial boundary

    j = 2 # chonProx 'speed'
    k = 1.5 # chonProx boundary dist

    def getProxIO(centroid):

        radialDist = (centroid[2]**2 + ((centroid[0]**2 + centroid[1]**2)**0.5 - 8)**2)**0.5 # distance to OD of bone fragments on plane z = 0, circle described by x^2 + y^2 = 64
        axialDist = abs(centroid[2]) # distance to the plane z = 0
        
        # Radial proximity function
        radialProx = 1 - (1/(1 + np.exp(-c*(radialDist - d))))
        # Axial proximity function
        axialProx = (1/(1 + np.exp(-a*(axialDist - b)))) - (1/(1 + np.exp(-u*(axialDist - v))))

        return radialProx * axialProx

    def getProxChon(centroid):

        R = (((centroid[0]**2 + centroid[1]**2)**0.5 - 8)**2 + (abs(centroid[2]) - b)**2)**0.5 # distance from OD of bone fragment tips

        # proximity function
        prox = 1 - (1/(1 + np.exp(-j*(R - k))))

        return prox





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

    for j in callusElementNumbers: # callusElementNumbers IDs, not indexing...

        nodes = p.element(j-1).items

        points = np.zeros(shape=(3, 4))

        for i in range(4):
            points[0, i] = p.node(nodes[i]-1).x # -1 because INDEXING starts at 0 even though IDs start at 1...
            points[1, i] = p.node(nodes[i]-1).y
            points[2, i] = p.node(nodes[i]-1).z

        centroid = np.mean(points, axis=1)

        proximityWeightsIO.append(getProxIO(centroid))
        proximityWeightsChon.append(getProxChon(centroid))

    return proximityWeightsIO, proximityWeightsChon