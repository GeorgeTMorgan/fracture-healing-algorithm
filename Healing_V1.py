import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
from subprocess import Popen
from shutil import copyfile
import pathlib
import ctypes, sys
import time

import preProcessing
import calcVolumes
import readResults
import readTorsion
import doFuzzyTest
import updateNeighbors
import calcNewMatProps
import writeNewMatProps
import make2Dmesh
import show2Dmesh
import makeNeighbors
import calcProximity
import bendingStiffness
import updateLoad
np.set_printoptions(threshold=sys.maxsize)



############### Notes ###############

# Edited controlsystem.py in SKFuzzy library - in CrispValueCalculator.defuzz() - 
# comment out Value error - "crisp output cannot be calculated", 
# replace with: return np.zeros(self.sim._array_shape, dtype=np.float64) (ie. return zero)
# because default output should be "no change"

# Visual Studio must be run as administrator

# Must use same version of Python as is used by the version of Marc where py_post and py_mentat are taken froms




############### Inputs and Set-up ###############

simName = '1High'

currentIteration = 0
numIterations = 150 # Number of algorithm iterations
#reactions = [] # Reaction moments from torsional rigidity test
IFMs = []

newLoad = 500

numNodesNeighbor = 1 # define the number of nodes which adjacent elements must share in order to be defined as neighboring elements (2 nodes for 2D linear tri3 elements, 6 nodes for 3D tet10 elements, or can change to 1 (includes 'diagonal' elements) or 2 (adjacent elements which share an edge) or 3 (adjacent elements which share a face) for 3D linear tet4 elements)
meshSizeCoefficient = 0.57
IFMDoF_String = "Displacement X" # For 2D
#IFMDoF_String = "Displacement Z" # For 3D


# Input Marc file name, and temp copy file name
fileName = '1High_2D_init'
#fileName = 'gap2lowstrainmesh07_torsion_job1' # select the .dat file, make sure it is in this directory (for now)
newFileName = 'FractureHealing_Temp'
#bendingFileName = 'bendingStiffness_mesh07_domainBig_job1'
#bendingNewFileName = 'bendingStiff_Temp'
#rigidityFileName = 'TibiaTorsion_job1'
#rigidityNewFileName = 'Torsion_Temp'

# Input folder to store results
currentSim = r'1High'

RBE2nodeID = 1

pathlib.Path(currentSim + r'\eachStep').mkdir(parents=True, exist_ok=True)  # create results directories
pathlib.Path(currentSim + r'\codeUsed').mkdir(parents=True, exist_ok=True)  # create code copy directory
pathlib.Path(currentSim + r'\FE').mkdir(parents=True, exist_ok=True)
pathlib.Path(currentSim + r'\gifs').mkdir(parents=True, exist_ok=True)
pathlib.Path(currentSim + r'\bending').mkdir(parents=True, exist_ok=True)



############### Pre-processing ###############

# Get callus element numbers and subgroups, as well as callus element connectivity dictionary
print("start")
print(simName)

callusElementNumbers, callusSurfaceElements, callusMedullarySurfaceElements, callusBoneSurfaceElements, callusConnectivity, callusConnectivityNodesDict = preProcessing.preProcessing(fileName, newFileName, numNodesNeighbor)

#don't need relative volumes if neighbor def'n is distance based
#volumesRelative = calcVolumes.calcVolumes(fileName, callusElementNumbers)
volumesRelative = calcVolumes.calcAreas(fileName, callusElementNumbers) # use same variable name for area (2D elements) but use calcArea function

np.save(currentSim + r'\callusElements.npy', callusElementNumbers)
np.save(currentSim + r'\callusSurfaceElements.npy', callusSurfaceElements)
np.save(currentSim + r'\callusMedullarySurfaceElements.npy', callusMedullarySurfaceElements)
np.save(currentSim + r'\callusBoneSurfaceElements.npy', callusBoneSurfaceElements)
np.save(currentSim + r'\callusConnectivityElements.npy', callusConnectivity)
np.save(currentSim + r'\volumesRelative.npy', volumesRelative)
#np.save(currentSim + r'\proximityWeightsIO.npy', proximityWeightsIO)
#np.save(currentSim + r'\proximityWeightsChon.npy', proximityWeightsChon)


## Restart Code
#copyfile(fileName + '.dat', newFileName + '.dat')

#callusElementNumbers = list(np.load(currentSim + r'\callusElements.npy', allow_pickle='TRUE'))
#callusSurfaceElements = list(np.load(currentSim + r'\callusSurfaceElements.npy', allow_pickle='TRUE'))
#callusMedullarySurfaceElements = list(np.load(currentSim + r'\callusMedullarySurfaceElements.npy', allow_pickle='TRUE'))
#callusBoneSurfaceElements = list(np.load(currentSim + r'\callusBoneSurfaceElements.npy', allow_pickle='TRUE'))
#callusConnectivity = np.load(currentSim + r'\callusConnectivityElements.npy', allow_pickle='TRUE')[()] # formatted as a dict with keys as each callus element ID, and the value as a list of neighboring elements
#volumesRelative = list(np.load(currentSim + r'\volumesRelative.npy', allow_pickle='TRUE'))
#proximityWeights = list(np.load(r'proximityWeights.npy', allow_pickle='TRUE'))
#IFMs = list(np.load(currentSim + r'\IFMs.npy', allow_pickle='TRUE'))

proximityWeightsIO = [1] * len(callusElementNumbers)
proximityWeightsChon = [1] * len(callusElementNumbers)
proximityWeightsChonInitial = [1] * len(callusElementNumbers)

# first row is E, second in nu
materials = np.zeros((len(callusElementNumbers),2))
materials[:,0] = 3.0
materials[:,1] = 0.3

print("finished loading")
# Initialize healing state variables: woven bone, fibrocartilage, and vascularity (soft tissue concentration is inferred becase Cwb + Cc + Cs = 1)
# Entries are for each element in the same order as callusElementNumbers, each state variable is between 0 and 1
# Initial values are given in Simon et al 2011 - 100% soft tissue (0 woven and 0 cartilage)
stateVariables = np.zeros((len(callusElementNumbers),3))
neighbors = np.zeros((len(callusElementNumbers),3)) # the maximum for each elements neighbors

# initialize matrix for storing temporal smoothing statevariable history, smoothed for the past trailingIterations iterations
trailingIterations = 7
smoothingSV = np.zeros((len(callusElementNumbers), 3, trailingIterations))

## Restart Code
#stateVariables = np.load('combined_method3\stateVariables49.npy', allow_pickle='TRUE')
#materials = np.load('combined_method3\materials49.npy', allow_pickle='TRUE')
#equivalentStrains = np.load('combined_method3\equivalentStrains49.npy', allow_pickle='TRUE')
#hydrostaticStrains = np.load('combined_method3\hydrostaticStrains49.npy', allow_pickle='TRUE')

neighbors = updateNeighbors.updateNeighbors(stateVariables, neighbors, callusElementNumbers, callusBoneSurfaceElements, callusMedullarySurfaceElements, callusConnectivity, currentIteration)


# Perform first iteration using already completed FE sim

# read callus element strains
# the elements are ordered in the same way as callusElementNumbers
hydrostaticStrains, equivalentStrains, IFM, minStrains, maxStrains = readResults.readResults(fileName, callusElementNumbers, materials, RBE2nodeID, IFMDoF_String)
print("IFM: " + str(IFM))
IFMs.append(IFM)
print("results read - current iteration " + str(currentIteration))

# do Fuzzy logics

stateVariables, activeRules, smoothingSV = doFuzzyTest.doFuzzyTest(callusElementNumbers, hydrostaticStrains, equivalentStrains, stateVariables, neighbors, volumesRelative, proximityWeightsIO, proximityWeightsChonInitial, minStrains, maxStrains, currentIteration, smoothingSV, meshSizeCoefficient)
print("fuzzy done - current iteration " + str(currentIteration))

neighbors = updateNeighbors.updateNeighbors(stateVariables, neighbors, callusElementNumbers, callusBoneSurfaceElements, callusMedullarySurfaceElements, callusConnectivity, currentIteration)
print("neighbors done - current iteration " + str(currentIteration))

# calc new mat props
materials = calcNewMatProps.calcNewMatProps(stateVariables, materials)
print("mat props calc'd - current iteration " + str(currentIteration))

# write new mat props to both FE files
writeNewMatProps.writeNewMatProps(newFileName, callusElementNumbers, materials)
#writeNewMatProps.writeNewMatProps(rigidityNewFileName, callusElementNumbers, materials)
print("mat props written - current iteration " + str(currentIteration))

# save state of each iteration
np.save(currentSim + r'\IFMs_' + currentSim + r'.npy', IFMs)
np.save(currentSim + r'\eachStep\activeRules' + str(currentIteration) + '.npy', activeRules)
np.save(currentSim + r'\eachStep\hydrostaticStrains' + str(currentIteration) + '.npy', hydrostaticStrains)
np.save(currentSim + r'\eachStep\equivalentStrains' + str(currentIteration) + '.npy', equivalentStrains)
np.save(currentSim + r'\eachStep\minStrains' + str(currentIteration) + '.npy', minStrains)
np.save(currentSim + r'\eachStep\maxStrains' + str(currentIteration) + '.npy', maxStrains)
np.save(currentSim + r'\eachStep\stateVariables' + str(currentIteration) + '.npy', stateVariables)
np.save(currentSim + r'\eachStep\materials' + str(currentIteration) + '.npy', materials)
np.save(currentSim + r'\eachStep\neighbors' + str(currentIteration) + '.npy', neighbors)

currentIteration += 1





################ Processing - Iterative loop ###############

while currentIteration < numIterations:
    print(simName)
    print("current iteration " + str(currentIteration))


    #updateLoad.updateLoad(newFileName, newLoad)
    ## increase load after 2 weeks
    #if currentIteration == 14:
    #    updateLoad.updateLoad(newFileName, newLoad)

    # run FE
    # VISUAL STUDIO MUST BE RUN AS ADMINISTRATOR FOR THIS TO WORK WITHOUT NEEDING UAC CONFIRMATION
    # ie, if not run as admin, this line will require a manual prompt response every iteration
    subprocess.call([r"C:\Program Files\MSC.Software\Marc\2021.2.0\marc2021.2\tools\run_marc.bat", r"-jid", r"FractureHealing_Temp.dat", r"-back", r"yes", r"-nts", r"2", r"-nte", r"2"])
    print("FE done - current iteration " + str(currentIteration))

    ## run virtual torsion test
    #subprocess.call([r"C:\Program Files\MSC.Software\Marc\2021.2.0\marc2021.2\tools\run_marc.bat", r"-jid", r"Torsion_Temp.dat", r"-back", r"yes", r"-nts", r"2", r"-nte", r"2"])
    #print("virtual torsion test complete - current iteration " + str(currentIteration))

    ## read and save torsion results in radians
    #reaction = readTorsion.readTorsion(rigidityNewFileName, RBE2nodeID)
    #print(reaction)
    #reactions.append(reaction)
    #print("torsional stability results read - current iteration " + str(currentIteration))

    # read callus element strains
    # the elements are ordered in the same way as callusElementNumbers
    hydrostaticStrains, equivalentStrains, IFM, minStrains, maxStrains = readResults.readResults(newFileName, callusElementNumbers, materials, RBE2nodeID, IFMDoF_String)
    print("IFM: " + str(IFM))
    IFMs.append(IFM)
    print("results read - current iteration " + str(currentIteration))

    # do Fuzzy logics
    stateVariables, activeRules, smoothingSV = doFuzzyTest.doFuzzyTest(callusElementNumbers, hydrostaticStrains, equivalentStrains, stateVariables, neighbors, volumesRelative, proximityWeightsIO, proximityWeightsChon, minStrains, maxStrains, currentIteration, smoothingSV, meshSizeCoefficient)
    print("fuzzy done - current iteration " + str(currentIteration))

    neighbors = updateNeighbors.updateNeighbors(stateVariables, neighbors, callusElementNumbers, callusBoneSurfaceElements, callusMedullarySurfaceElements, callusConnectivity, currentIteration)
    print("neighbors done - current iteration " + str(currentIteration))

    # calc new mat props
    materials = calcNewMatProps.calcNewMatProps(stateVariables, materials)
    print("mat props calc'd - current iteration " + str(currentIteration))

    # write new mat props to both FE files
    writeNewMatProps.writeNewMatProps(newFileName, callusElementNumbers, materials)
    #writeNewMatProps.writeNewMatProps(rigidityNewFileName, callusElementNumbers, materials)
    print("mat props written - current iteration " + str(currentIteration))

    # save state of each iteration
    np.save(currentSim + r'\IFMs_' + currentSim + r'.npy', IFMs)
    np.save(currentSim + r'\eachStep\activeRules' + str(currentIteration) + '.npy', activeRules)
    np.save(currentSim + r'\eachStep\hydrostaticStrains' + str(currentIteration) + '.npy', hydrostaticStrains)
    np.save(currentSim + r'\eachStep\equivalentStrains' + str(currentIteration) + '.npy', equivalentStrains)
    np.save(currentSim + r'\eachStep\minStrains' + str(currentIteration) + '.npy', minStrains)
    np.save(currentSim + r'\eachStep\maxStrains' + str(currentIteration) + '.npy', maxStrains)
    np.save(currentSim + r'\eachStep\stateVariables' + str(currentIteration) + '.npy', stateVariables)
    np.save(currentSim + r'\eachStep\materials' + str(currentIteration) + '.npy', materials)
    np.save(currentSim + r'\eachStep\neighbors' + str(currentIteration) + '.npy', neighbors)

    # save a copy of the .dat and .t16 files every 7 iterations
    if currentIteration % 7 == 6:
        copyfile(newFileName + '.dat', os.path.join(currentSim, 'FE', newFileName + str(currentIteration) + '.dat'))
        copyfile(newFileName + '.t16', os.path.join(currentSim, 'FE', newFileName + str(currentIteration) + '.t16'))

    currentIteration += 1

print(IFMs)

######## POST PROCESSING #######

# create gifs of callus iterations for all relevant calculated values

#nodes2D, elIDsTri2D, elIDsQuad2D, elConnectivityTri2D, elConnectivityQuad2D = make2Dmesh.make2Dmesh(fileName, callusElementNumbers)
nodes2D, elIDsTri2D, elIDsQuad2D, elConnectivityTri2D, elConnectivityQuad2D = make2Dmesh.makeAxisymmetricMesh(fileName, callusElementNumbers, callusConnectivityNodesDict)
np.save(currentSim + r'\gifs\nodes2D.npy', nodes2D)
np.save(currentSim + r'\gifs\elIDsTri2D.npy', elIDsTri2D)
np.save(currentSim + r'\gifs\elIDsQuad2D.npy', elIDsQuad2D)
np.save(currentSim + r'\gifs\elConnectivityTri2D.npy', elConnectivityTri2D)
np.save(currentSim + r'\gifs\elConnectivityQuad2D.npy', elConnectivityQuad2D)

show2Dmesh.show2Dmesh(callusElementNumbers, nodes2D, elIDsTri2D, elIDsQuad2D, elConnectivityTri2D, elConnectivityQuad2D, numIterations, currentSim, currentSim + r'\gifs')

# calc the bending stiffness of the final callus iteration

#bendingStiffness.bendingStiffness(bendingFileName, bendingNewFileName, callusElementNumbers, materials)

# delete FractureHealing_Temp files
os.remove(newFileName + '.dat')
os.remove(newFileName + '.log')
os.remove(newFileName + '.out')
os.remove(newFileName + '.sts')
os.remove(newFileName + '.t16')

# copy all input files used to codeUsed folder - do this at the start?