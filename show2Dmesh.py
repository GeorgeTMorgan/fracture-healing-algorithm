import meshpy.triangle as triangle
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.collections
from scipy.interpolate import griddata
import builtins
import imageio

def show2Dmesh(callusElementNumbers, nodes2D, elIDsTri2D, elIDsQuad2D, elConnectivityTri2D, elConnectivityQuad2D, numIterations, currentSim, saveDir):

    resultsFolder = currentSim + r'\eachStep'

    width = 12
    height = 8

    #plotValsTri = []
    #plotValsQuad = []

    #for elID in elIDsTri2D:
    #    plotValsTri.append(materials[callusElementNumbers.index(elID), 0])

    #for elID in elIDsQuad2D:
    #    plotValsQuad.append(materials[callusElementNumbers.index(elID), 0])

    mesh_points = np.array(nodes2D)
    mesh_tris = np.asarray(elConnectivityTri2D)
    mesh_quads = np.asarray(elConnectivityQuad2D)

    y = np.array(mesh_points[:, 0])
    z = np.array(mesh_points[:, 1])

    def quatplot(y, z, quatrangles, values, clim, ax=None, **kwargs):

        if not ax: ax=plt.gca()
        yz = np.c_[y, z]
        verts= yz[quatrangles]
        pc = matplotlib.collections.PolyCollection(verts, **kwargs)
        pc.set_array(values)
        pc.set_clim(clim[0],clim[1]) # set the range of the colorbar
        ax.add_collection(pc)
        ax.autoscale()
        return pc

    ##fig, ax = plt.subplots()
    #fig = plt.figure(figsize=(3.5, 4))
    #ax = fig.add_subplot(1, 1, 1)

    #plt.rc('font', family='serif')
    ##plt.rc('xtick', labelsize='small')
    ##plt.rc('ytick', labelsize='small')

    #pc3 = quatplot(y, z, np.asarray(mesh_tris), np.array(plotValsTri), ax=ax, 
    #         edgecolor="black", cmap="rainbow")
    #pc4 = quatplot(y, z, np.asarray(mesh_quads), np.array(plotValsQuad), ax=ax, 
    #         edgecolor="black", cmap="rainbow")
    #fig.colorbar(pc3, ax=ax) 
    ##fig.colorbar(pc4, ax=ax)

    #ax.set_xlim([4, 16])
    #ax.set_xticks([])
    #ax.set_yticks([])

    #plt.tight_layout()
    #plt.show()





    #fig, ax = plt.subplots()
    #plt.rc('font', family='serif')

    #fig = plt.figure(figsize=(3.5, 4))
    #ax.set_xlim([4, 16])
    #ax.set_xticks([])
    #ax.set_yticks([])

    #plt.rc('font', family='serif')
    #plt.rc('xtick', labelsize='small')
    #plt.rc('ytick', labelsize='small')
             

    def plotIteration(iteration, data, clim):

        #materials = np.load(resultsFolder + r'\materials' + str(iteration) + '.npy', allow_pickle='TRUE')
        #SVs = np.load(resultsFolder + r'\stateVariables' + str(iteration) + '.npy', allow_pickle='TRUE')
        #activeRules = np.load(resultsFolder + r'\activeRules' + str(iteration) + '.npy', allow_pickle='TRUE')
        #equivStrains = np.load(resultsFolder + r'\equivalentStrains' + str(iteration) + '.npy', allow_pickle='TRUE')
        #hydroStrains = np.load(resultsFolder + r'\hydrostaticStrains' + str(iteration) + '.npy', allow_pickle='TRUE')

        #data = np.load(resultsFolder + dataStr + str(iteration) + '.npy', allow_pickle='TRUE')

        plotValsTri = []
        plotValsQuad = []

        for elID in elIDsTri2D:
            #plotValsTri.append(hydroStrains[callusElementNumbers.index(elID)])
            plotValsTri.append(data[callusElementNumbers.index(elID)])

        for elID in elIDsQuad2D:
            #plotValsQuad.append(hydroStrains[callusElementNumbers.index(elID)])
            plotValsQuad.append(data[callusElementNumbers.index(elID)])

        pc3 = quatplot(y, z, np.asarray(mesh_tris), np.array(plotValsTri), clim, ax=ax,
                       edgecolor="black", cmap="rainbow")
        if len(plotValsQuad) != 0:
            pc4 = quatplot(y, z, np.asarray(mesh_quads), np.array(plotValsQuad), clim, ax=ax, 
                           edgecolor="black", cmap="rainbow")

        if(iteration == 0):
            fig.colorbar(pc3, ax=ax) 
        fig.canvas.draw()

        img = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))

        return img

    SVs_bone_all = []
    SVs_cart_all = []
    activeRules_io_all = []
    activeRules_chon_all = []
    activeRules_cc_all = []
    activeRules_eo_all = []
    equivStrains_all = []
    hydroStrains_all = []
    minStrains_all = []
    maxStrains_all = []
    youngsMods_all = []

    for iteration in range(numIterations):

        SVs = np.load(resultsFolder + r'\stateVariables' + str(iteration) + '.npy', allow_pickle='TRUE')
        activeRules = np.load(resultsFolder + r'\activeRules' + str(iteration) + '.npy', allow_pickle='TRUE')
        equivStrains = np.load(resultsFolder + r'\equivalentStrains' + str(iteration) + '.npy', allow_pickle='TRUE')
        hydroStrains = np.load(resultsFolder + r'\hydrostaticStrains' + str(iteration) + '.npy', allow_pickle='TRUE')
        minStrains = np.load(resultsFolder + r'\minStrains' + str(iteration) + '.npy', allow_pickle='TRUE')
        maxStrains = np.load(resultsFolder + r'\maxStrains' + str(iteration) + '.npy', allow_pickle='TRUE')
        materials = np.load(resultsFolder + r'\materials' + str(iteration) + '.npy', allow_pickle='TRUE')

        SVs_bone_all.append(SVs[:,0])
        SVs_cart_all.append(SVs[:,1])
        activeRules_io_all.append(activeRules[0,:])
        activeRules_chon_all.append(activeRules[1,:])
        activeRules_cc_all.append(activeRules[2,:])
        activeRules_eo_all.append(activeRules[3,:])
        equivStrains_all.append(equivStrains)
        hydroStrains_all.append(hydroStrains)
        minStrains_all.append(minStrains)
        maxStrains_all.append(maxStrains)
        youngsMods_all.append(materials[:,0])

    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\bone.gif', [plotIteration(n, SVs_bone_all[n], [0, 1]) for n in range(numIterations)], fps=4)
    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\cart.gif', [plotIteration(n, SVs_cart_all[n], [0, 1]) for n in range(numIterations)], fps=4)
    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\eo.gif', [plotIteration(n, activeRules_eo_all[n], [0, 1]) for n in range(numIterations)], fps=4)
    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\equivStrains.gif', [plotIteration(n, equivStrains_all[n], [0, 0.05]) for n in range(numIterations)], fps=4)
    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\hydroStrains.gif', [plotIteration(n, hydroStrains_all[n], [-0.01, 0.01]) for n in range(numIterations)], fps=4)
    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\minStrains.gif', [plotIteration(n, minStrains_all[n], [-0.1, 0]) for n in range(numIterations)], fps=4)
    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\maxStrains.gif', [plotIteration(n, maxStrains_all[n], [0, 0.1]) for n in range(numIterations)], fps=4)
    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\youngsMods.gif', [plotIteration(n, youngsMods_all[n], [0, 4000]) for n in range(numIterations)], fps=4)
    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\io.gif', [plotIteration(n, activeRules_io_all[n], [0, 1]) for n in range(numIterations)], fps=4)
    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\chon.gif', [plotIteration(n, activeRules_chon_all[n], [0, 1]) for n in range(numIterations)], fps=4)
    fig, ax = plt.subplots(figsize=(width,height))
    plt.rc('font', family='serif')
    imageio.mimsave(saveDir + r'\cc.gif', [plotIteration(n, activeRules_cc_all[n], [0, 1]) for n in range(numIterations)], fps=4)