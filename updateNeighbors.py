import numpy as np
import time

def updateNeighbors(stateVariables, neighbors, callusElementNumbers, callusBoneSurfaceElements, callusMedullarySurfaceElements, callusConnectivity, currentIteration):

    def calculateMaxNeighborConc(properties_array, neighbor_indices, idx0):
    
        # Create an array to store the properties of neighboring elements
        neighbor_properties_array = properties_array[np.concatenate(list(neighbor_indices.values()))-idx0] # -idx0 is a quick hack to get around the dictionary keys not actually being the indexes for callusElementNumbers

        # Use a list to store arrays of neighboring properties for each element
        neighbor_properties_list = np.split(neighbor_properties_array, np.cumsum([len(neighbors) for neighbors in neighbor_indices.values()])[:-1])

        max_values = np.maximum.reduceat(neighbor_properties_array, np.concatenate([[0], np.cumsum([len(neighbors) for neighbors in neighbor_indices.values()])[:-1]]))

        # Create the final dictionary X using NumPy functions
        #X = dict(zip(neighbor_indices.keys(), max_values))

        return max_values

    properties_array_bone = stateVariables[:,0]
    properties_array_cart = stateVariables[:,1]

    neighbors[:,0] = calculateMaxNeighborConc(properties_array_bone, callusConnectivity, callusElementNumbers[0])
    neighbors[:,1] = calculateMaxNeighborConc(properties_array_cart, callusConnectivity, callusElementNumbers[0])

    # re-enforce BCs here
    # need to set nBone of callusBoneSurfaceElements to 1
    indices = np.where(np.isin(callusElementNumbers, callusBoneSurfaceElements))[0]
    neighbors[indices,0] = 1

    return neighbors
