import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
from subprocess import Popen
from shutil import copyfile
import pathlib
import ctypes, sys
import time
from scipy.interpolate import griddata
from py_post import *

newFileName = 'FractureHealing_Temp146'
bendingFileName = '1high_bending_job1'
bendingNewFileName = 'bendingStiff_Temp'
callusElementNumbers = list(np.load(r'callusElements.npy', allow_pickle='TRUE'))
materials =  np.load('materials149.npy', allow_pickle='TRUE')


startTime = time.time()


##### Load the 2D mesh data (using py_post and points should be centroids already)

print(time.time() - startTime)

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

points = np.zeros((3, 2))
centroids2D = []

for idx, j in enumerate(callusElementNumbers):

    el = p.element(j-1)
    elNodes = el.items

    for i, nodeID in enumerate(elNodes):

        points[i,0] = p.node(nodeID-1).x # -1 because INDEXING starts at 0 even though IDs start at 1...
        points[i,1] = p.node(nodeID-1).y

    centroids2D.append(np.mean(points, axis=0))

centroids2D = np.array(centroids2D)
p.close()


##### Load 3D mesh connectivity data

print(time.time() - startTime)

numMaterials = 0

with open(bendingFileName + '.dat', 'r') as f:
    searchlines = f.readlines()

for i, line in enumerate(searchlines):

    if "define              element             set                 callus\n" in line:
        callusElementsLowerBound = int(searchlines[i+1].split()[0])
        callusElementsUpperBound = int(searchlines[i+1].split()[2]) 
        callusElements3D = list(range(callusElementsLowerBound,callusElementsUpperBound + 1))

for i, line in enumerate(searchlines):

    if "connectivity" in line:
        if int(searchlines[i+2].split()[0]) == callusElements3D[0]:
            callusConnectivityLine = i+1

    if "isotropic" in line:
        lastIsotropicLine = i+1
        numMaterials = numMaterials + 1


##### Calculate 3D mesh data for interpolation (centroid)

print(time.time() - startTime)

p = post_open(bendingFileName + ".t16")

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

points = np.zeros((4, 3))
centroids3DCartesian = []

for idx, j in enumerate(callusElements3D):

    el = p.element(j-1)
    elNodes = el.items

    for i, nodeID in enumerate(elNodes):

        points[i,0] = p.node(nodeID-1).x # -1 because INDEXING starts at 0 even though IDs start at 1...
        points[i,1] = p.node(nodeID-1).y
        points[i,2] = p.node(nodeID-1).z

    centroids3DCartesian.append(np.mean(points, axis=0))

centroids3DCartesian = np.array(centroids3DCartesian)
p.close()

# Map the 3D centroids to an axisymmetric x-y coordinate system matching the 2D system
radialXCoords = centroids3DCartesian[:,0]
radialYCoords = np.sqrt((centroids3DCartesian[:,1])**2 + (centroids3DCartesian[:,2])**2)
centroids3D = np.column_stack((radialXCoords, radialYCoords))


##### Interpolate the 3D mesh data material values

print(time.time() - startTime)

# Linear interpolation, and use nearest for any points falling outside the interpolation region
youngsMods3D_nearest = griddata(centroids2D, materials[:,0], centroids3D, method='nearest')
youngsMods3D = griddata(centroids2D, materials[:,0], centroids3D, method='linear')
youngsMods3D = np.where(np.isnan(youngsMods3D), youngsMods3D_nearest, youngsMods3D)

# Linear interpolation, and use nearest for any points falling outside the interpolation region
poissonsRatios3D_nearest = griddata(centroids2D, materials[:,1], centroids3D, method='nearest')
poissonsRatios3D = griddata(centroids2D, materials[:,1], centroids3D, method='linear')
poissonsRatios3D = np.where(np.isnan(poissonsRatios3D), poissonsRatios3D_nearest, poissonsRatios3D)


##### Initialize and simultaneously write the new mat props for the 3D model using DAT editing

print(time.time() - startTime)

copyfile(bendingFileName + '.dat', bendingNewFileName + '.dat')

with open(bendingNewFileName + '.dat', 'r') as f:
    searchlines = f.readlines()

# start with mat prop stuff to not mess up the line numbering for the connectivity
isotropicChunk = searchlines[lastIsotropicLine - 1: lastIsotropicLine + 10] # example chunk of isotropic structure
outputChunks = []
lineToChange = isotropicChunk[2].split()
matLineToChange = isotropicChunk[3].split()
materialsCounter = numMaterials + 1
spacings = [17, 43, 10, 10, 23]

for i, elID in enumerate(callusElements3D):

    tempLine = lineToChange
    tempLine[0] = str(materialsCounter) + 'elastic'
    materialsCounter += 1
    tempLine[4] = '0callus' + '{:07}'.format(elID)
    formattedLine = ''.join(number.rjust(spacing) for number, spacing in zip(tempLine, spacings)) + '\n'

    updatedE = '{:.15e}'.format(youngsMods3D[i])
    updatedE = updatedE[:-4] + updatedE[-3] + updatedE[-1]
    updatedNu = '{:.15e}'.format(poissonsRatios3D[i])
    updatedNu = updatedNu[:-4] + updatedNu[-3] + updatedNu[-1]
    tempMatLine = matLineToChange
    tempMatLine[0] = updatedE
    tempMatLine[1] = updatedNu
    tempMatLine = ' ' + ' '.join(tempMatLine) + '\n'


    outputChunks.extend(isotropicChunk[0:2] + [formattedLine] + [tempMatLine] + isotropicChunk[4:])

searchlines[lastIsotropicLine + 10:lastIsotropicLine + 10] = outputChunks



# Connectivity stuff

connectivityLine = searchlines[callusConnectivityLine - 1] # funky line number indexing
initialCallusConn = searchlines[callusConnectivityLine]

outputLines = []
tempIDLine = initialCallusConn.split()
matID = int(searchlines[callusConnectivityLine].split()[5])
matID +=1

numCallusElements = len(callusElements3D)

for line in searchlines[(callusConnectivityLine + 1):(callusConnectivityLine + numCallusElements + 1)]:
    outputLines.append(connectivityLine)

    matIDLine = tempIDLine
    matIDLine[5] = str(matID)
    formatted_line = ''.join(number.rjust(10) for number in matIDLine) + '\n'
    matID += 1
    outputLines.extend([formatted_line, line])

del searchlines[callusConnectivityLine - 1:callusConnectivityLine + numCallusElements + 1]
searchlines[callusConnectivityLine - 1:callusConnectivityLine - 1] = outputLines

# rewrite the .dat file

with open(bendingNewFileName + '.dat', 'w') as f:
    f.writelines(searchlines)



##### Run the new bending stiffness file

print(time.time() - startTime)

subprocess.call([r"C:\Program Files\MSC.Software\Marc\2021.2.0\marc2021.2\tools\run_marc.bat", r"-jid", r"bendingStiff_Temp.dat", r"-back", r"no", r"-nts", r"2", r"-nte", r"2"])

print(time.time() - startTime)