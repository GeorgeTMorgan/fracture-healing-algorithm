import numpy as np
import numpy.linalg as la
from py_post import *
import builtins

def make2Dmesh(filename, callusElementNumbers):

    # create a postpy obj using post_open
    p = post_open(filename + ".t16")

    # try opening the results file to check for errors
    try:
        p.moveto(1)
    except:
        print("Error opening post file: " + fileName + ".t16")

    # move to the last increment
    numIncrements = p.increments() # number of increments in the sim
    p.extrapolation("linear") # set extrapolation to linear
    p.moveto(numIncrements - 1) # move to the final increment

    nodes2D = [] # nodes (as coords like: [x, y]) in the 2D plane
    elements2D = {} # dict of elements connectivity in the 2D plane

    for j in callusElementNumbers: # callusElementNumbers IDs, not indexing...

        nodes = p.element(j-1).items

        points = np.zeros(shape=(4, 4))

        for i in range(4):
            points[0, i] = p.node(nodes[i]-1).x # -1 because INDEXING starts at 0 even though IDs start at 1...
            points[1, i] = p.node(nodes[i]-1).y
            points[2, i] = p.node(nodes[i]-1).z
            points[3, i] = 1

        # check if the element has at least one point (only checking x coord) < x=0 AND one > x=0, add element number to dict CUTPLANEELEMENTNUMBERS as key
        if any(x < 0 for x in points[0, :]) and any(x > 0 for x in points[0, :]):
            edges = [[points[0:3, 0], points[0:3, 1]],
                     [points[0:3, 0], points[0:3, 2]],
                     [points[0:3, 0], points[0:3, 3]],
                     [points[0:3, 1], points[0:3, 2]],
                     [points[0:3, 1], points[0:3, 3]],
                     [points[0:3, 2], points[0:3, 3]]]

            edges_int = [] # a list of edges which intersect the plane x = 0

            for edge in edges:
                if edge[0][0] * edge[1][0] < 0: # if the product of the x_coords is negative, then the edge crosses the plane x = 0
                    edges_int.append(edge)

            nodes_int = [] # the intersection points of the edges which intersect x = 0, these are also the nodes for the 2D plot on x = 0

            # calc the point each edge intersects x = 0
            for edge in edges_int:
                scale = edge[0][0] / (edge[1][0] - edge[0][0])
                y_new = edge[0][1] - (scale * (edge[1][1] - edge[0][1]))
                z_new = edge[0][2] - (scale * (edge[1][2] - edge[0][2]))
                nodes_int.append([y_new, z_new])

            nodeIDs = [] # the 2D nodeIDs of the current element
            tolerance = 0.00001

            for node in nodes_int:
                added = 0
                for coords in nodes2D:
                    dist = ((node[0] - coords[0])**2 + (node[1] - coords[1])**2)**0.5
                    if dist < tolerance: # if node to be added is close (within tolerance) to an existing 2D node, then add the existing nodeID to nodeIDs
                        nodeIDs.append(nodes2D.index(coords))
                        added = 1
                        break
                if added == 0: # if node to be added is not close to any existing nodes, then add it to the 2D global node list and add to nodeIDs
                    nodes2D.append(node)
                    nodeIDs.append(len(nodes2D)-1)

            elements2D[j] = builtins.set(nodeIDs)

    p.close()

    def orderQuad(nodes):

        # line segment intersection checker from https://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
        # Returns true is the 3 nodes A, B, and C are counter clockwise
        def ccw(A,B,C):
            return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0]) # if slope of AB is less than slope of AC (negative slopes corrected by multiplying out denominators), then points are ccw

        # Return true if line segments AB and CD intersect
        def intersect(A,B,C,D):
            return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

        nodesOrdered = []
        for node in nodes:
            nodesOrdered.append(node)

        # if not ordered correctly, swap 1st and 2nd nodes and check again
        if not intersect(nodes2D[nodesOrdered[0]], nodes2D[nodesOrdered[2]], nodes2D[nodesOrdered[1]], nodes2D[nodesOrdered[3]]): # nodes listed by diagonal segments (1-3, 2-4)
            nodesOrdered[1], nodesOrdered[0] = nodesOrdered[0], nodesOrdered[1]
            # if still not ordered corectly, swap back 1st and 2nd, then swap 2nd and 3rd, then order must be correct
            if not intersect(nodes2D[nodesOrdered[0]], nodes2D[nodesOrdered[2]], nodes2D[nodesOrdered[1]], nodes2D[nodesOrdered[3]]):
                nodesOrdered[2], nodesOrdered[1], nodesOrdered[0] = nodesOrdered[0], nodesOrdered[2], nodesOrdered[1]

        return nodesOrdered

    elIDsTri2D = []
    elIDsQuad2D = []
    elConnectivityTri2D = []
    elConnectivityQuad2D = []

    for elID, nodes in elements2D.items():
        if len(nodes) == 3:
            nodesToAdd = []
            for node in nodes:
                nodesToAdd.append(node)
            elIDsTri2D.append(elID)
            elConnectivityTri2D.append(np.array(nodesToAdd))
        elif len(nodes) == 4:
            nodesToAdd = orderQuad(nodes)
            elIDsQuad2D.append(elID)
            elConnectivityQuad2D.append(np.array(nodesToAdd))

    return nodes2D, elIDsTri2D, elIDsQuad2D, elConnectivityTri2D, elConnectivityQuad2D



def makeAxisymmetricMesh(filename, callusElementNumbers, callusConnectivityNodesDict):

    # create a postpy obj using post_open
    p = post_open(filename + ".t16")

    # try opening the results file to check for errors
    try:
        p.moveto(1)
    except:
        print("Error opening post file: " + fileName + ".t16")

    # move to the last increment
    numIncrements = p.increments() # number of increments in the sim
    p.extrapolation("linear") # set extrapolation to linear
    p.moveto(numIncrements - 1) # move to the final increment

    numNodes = p.nodes() # the number of nodes in the increment

    tolerance = 0.0000001

    nodes2D = [] # nodes (as coords like: [x, y]) in the 2D plane
    elements2D = {} # dict of elements connectivity in the 2D plane
    elIDsTri2D = []
    elIDsQuad2D = [] # empty for 2D model
    elConnectivityTri2D = []
    elConnectivityQuad2D = [] # empty for 2D model

    for j in callusElementNumbers:

        nodes_int = []
        nodeIDs = []

        nodes = p.element(j-1).items

        for node in nodes:
            nodes_int.append([p.node(node - 1).x, p.node(node - 1).y])

        for node in nodes_int:
            added = 0
            for coords in nodes2D:
                dist = ((node[0] - coords[0])**2 + (node[1] - coords[1])**2)**0.5
                if dist < tolerance: # if node to be added is close (within tolerance) to an existing 2D node, then add the existing nodeID to nodeIDs
                    nodeIDs.append(nodes2D.index(coords))
                    added = 1
                    break
            if added == 0: # if node to be added is not close to any existing nodes, then add it to the 2D global node list and add to nodeIDs
                nodes2D.append(node)
                nodeIDs.append(len(nodes2D)-1)

        elements2D[j] = builtins.set(nodeIDs)

    p.close()

    for elID, nodes in elements2D.items():
        nodesToAdd = []
        for node in nodes:
            nodesToAdd.append(node)
        elIDsTri2D.append(elID)
        elConnectivityTri2D.append(np.array(nodesToAdd))

    return nodes2D, elIDsTri2D, elIDsQuad2D, elConnectivityTri2D, elConnectivityQuad2D


