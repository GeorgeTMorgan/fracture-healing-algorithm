import numpy as np
import skfuzzy as fuzz
from skfuzzy import control as ctrl

def doFuzzyTest(callusElementNumbers, hydrostaticStrains, equivalentStrains, stateVariablesInput, neighborsInput, volumesRelative, proximityWeightsIO, proximityWeightsChon, minStrains, maxStrains, currentIteration, smoothingSVInput, meshSizeCoefficientInput):

    stateVariables = stateVariablesInput
    neighbors = neighborsInput
    smoothingSV = smoothingSVInput
    meshSizeCoefficient = meshSizeCoefficientInput

    # shift all in smoothingSV (this will delete oldest iteration in the process)
    for i in range(smoothingSV.shape[2] - 1):
        smoothingSV[:,0,i] = smoothingSV[:,0,i+1] # c_bone
        smoothingSV[:,1,i] = smoothingSV[:,1,i+1] # c_cart
        smoothingSV[:,2,i] = smoothingSV[:,2,i+1] # c_perfusion

    # Define antecedent fuzzy sets and membership functions
    CartilageUniv = np.arange(0, 101, 1)
    nCartUniv = np.arange(0, 101, 1)
    BoneUniv = np.arange(0, 101, 1)
    nBoneUniv = np.arange(0, 101, 1)

    Cartilage_low = nCart_low = Bone_low = nBone_low = fuzz.trapmf(BoneUniv, [0, 0, 20, 40])
    Cartilage_med = nCart_med = Bone_med = nBone_med = fuzz.trapmf(BoneUniv, [20, 40, 60, 80])
    Cartilage_high = nCart_high = Bone_high = nBone_high = fuzz.trapmf(BoneUniv, [60, 80, 100, 100])

    DistortionalUniv = np.arange(0, 20.01, 0.01)
    Distortional_zero = fuzz.trapmf(DistortionalUniv, [0, 0, 0.01, 0.05])
    Distortional_low = fuzz.trapmf(DistortionalUniv, [0.01, 0.05, 5, 7])
    Distortional_med = fuzz.trapmf(DistortionalUniv, [5, 7, 10, 10])
    Distortional_dest = fuzz.trapmf(DistortionalUniv, [10, 10, 20, 20])
    #CHANGED HERE AND IN INPUTS SECTION
    #Distortional_med = fuzz.trapmf(DistortionalUniv, [5, 9, 14, 100])
    #Distortional_dest = fuzz.trapmf(DistortionalUniv, [14, 100, 100, 100])

    DilatationalUniv = np.arange(-8, 8.01, 0.01)
    Dilatational_negDes = fuzz.trapmf(DilatationalUniv, [-8, -8, -6, -4])
    Dilatational_negMed = fuzz.trapmf(DilatationalUniv, [-6, -4, -0.9, -0.8])
    Dilatational_negLow = fuzz.trapmf(DilatationalUniv, [-0.9, -0.8, -0.03, -0.01])
    Dilatational_zero = fuzz.trapmf(DilatationalUniv, [-0.03, -0.01, 0.01, 0.03])
    Dilatational_posLow = fuzz.trapmf(DilatationalUniv, [0.01, 0.03, 0.8, 0.9])
    Dilatational_posMed = fuzz.trapmf(DilatationalUniv, [0.8, 0.9, 4, 6])
    Dilatational_posDes = fuzz.trapmf(DilatationalUniv, [4, 6, 8, 8])

    maxUniv = np.arange(0, 50.01, 0.01)
    #max_zero = fuzz.trapmf(maxUniv, [0, 0, 0.5, 1])
    max_low = fuzz.trapmf(maxUniv, [0, 0, 0.5, 5])
    max_med = fuzz.trapmf(maxUniv, [0.5, 5, 5, 50])
    max_high = fuzz.trapmf(maxUniv, [5, 50, 50, 50])
    #max_dest = fuzz.trapmf(maxUniv, [1, 5, 50, 50])

    minUniv = np.arange(0, 50.01, 0.01)
    min_zero = fuzz.trapmf(minUniv, [0, 0, 0.5, 1])
    min_low = fuzz.trapmf(minUniv, [0.5, 1, 2, 2.5])
    min_med = fuzz.trapmf(minUniv, [2, 2.5, 7.5, 10])
    #min_low = fuzz.trapmf(minUniv, [0.5, 1, 7.5, 10])
    #min_med = fuzz.trapmf(minUniv, [7.5, 10, 10, 12.5])
    #min_high = fuzz.trapmf(minUniv, [10, 12.5, 50, 50])
    min_high = fuzz.trapmf(minUniv, [7.5, 10, 50, 50])
    #min_dest = fuzz.trapmf(minUniv, [14, 16, 50, 50])
    
    # Set rule rates
    #R_io_bone = 21.6
    R_io_bone = 100
    R_chon_bone = 0
    R_cc_bone = 100
    R_eo_bone = 100
    R_at1_bone = -100
    R_at2_bone = -100
    R_remod_bone = -100

    R_io_cart = 0
    R_chon_cart = 100
    R_cc_cart = -100
    R_eo_cart = -100
    R_at1_cart = -100
    R_at2_cart = -100
    R_remod_cart = 0

    ruleRatesBone = np.array([R_io_bone, R_chon_bone, R_cc_bone, R_eo_bone, R_at1_bone, R_at2_bone, R_remod_bone])
    ruleRatesCart = np.array([R_io_cart, R_chon_cart, R_cc_cart, R_eo_cart, R_at1_cart, R_at2_cart, R_remod_cart])
    

    ########## VECTORIZED SECTION ##########

    # Set antecedents values

    cartilage_in = stateVariables[:,1]*100
    nCart_in = neighbors[:,1]*100
    bone_in = stateVariables[:,0]*100
    nBone_in = neighbors[:,0]*100
    distortional_in = np.minimum(equivalentStrains[:]*100,20) # bring dilatational strain into the universe of discourse
    dilatational_in = np.maximum(np.minimum(hydrostaticStrains[:]*100,8),-8) # bring dilatational strain into the universe of discourse
    max_in = np.minimum(maxStrains[:]*100, 50)
    min_in = np.minimum(abs(minStrains[:])*100, 50)

    # Calculate antecedent membership values

    cartMem_low = fuzz.interp_membership(CartilageUniv, Cartilage_low, cartilage_in)
    cartMem_med = fuzz.interp_membership(CartilageUniv, Cartilage_med, cartilage_in)
    cartMem_high = fuzz.interp_membership(CartilageUniv, Cartilage_high, cartilage_in)

    nCartMem_low = fuzz.interp_membership(nCartUniv, nCart_low, nCart_in)
    nCartMem_med = fuzz.interp_membership(nCartUniv, nCart_med, nCart_in)
    nCartMem_high = fuzz.interp_membership(nCartUniv, nCart_high, nCart_in)

    boneMem_low = fuzz.interp_membership(BoneUniv, Bone_low, bone_in)
    boneMem_med = fuzz.interp_membership(BoneUniv, Bone_med, bone_in)
    boneMem_high = fuzz.interp_membership(BoneUniv, Bone_high, bone_in)

    nBoneMem_low = fuzz.interp_membership(nBoneUniv, nBone_low, nBone_in)
    nBoneMem_med = fuzz.interp_membership(nBoneUniv, nBone_med, nBone_in)
    nBoneMem_high = fuzz.interp_membership(nBoneUniv, nBone_high, nBone_in)

    distMem_zero = fuzz.interp_membership(DistortionalUniv, Distortional_zero, distortional_in)
    distMem_low = fuzz.interp_membership(DistortionalUniv, Distortional_low, distortional_in)
    distMem_med = fuzz.interp_membership(DistortionalUniv, Distortional_med, distortional_in)
    distMem_dest = fuzz.interp_membership(DistortionalUniv, Distortional_dest, distortional_in)

    dilMem_negDes = fuzz.interp_membership(DilatationalUniv, Dilatational_negDes, dilatational_in)
    dilMem_negMed = fuzz.interp_membership(DilatationalUniv, Dilatational_negMed, dilatational_in)
    dilMem_negLow = fuzz.interp_membership(DilatationalUniv, Dilatational_negLow, dilatational_in)
    dilMem_zero = fuzz.interp_membership(DilatationalUniv, Dilatational_zero, dilatational_in)
    dilMem_posLow = fuzz.interp_membership(DilatationalUniv, Dilatational_posLow, dilatational_in)
    dilMem_posMed = fuzz.interp_membership(DilatationalUniv, Dilatational_posMed, dilatational_in)
    dilMem_posDes = fuzz.interp_membership(DilatationalUniv, Dilatational_posDes, dilatational_in)

    #maxMem_zero = fuzz.interp_membership(maxUniv, max_zero, max_in)
    maxMem_low = fuzz.interp_membership(maxUniv, max_low, max_in)
    maxMem_med = fuzz.interp_membership(maxUniv, max_med, max_in)
    maxMem_high = fuzz.interp_membership(maxUniv, max_high, max_in)
    #maxMem_dest = fuzz.interp_membership(maxUniv, max_dest, max_in)

    minMem_zero = fuzz.interp_membership(minUniv, min_zero, min_in)
    #minMem_veryLow = fuzz.interp_membership(minUniv, min_veryLow, min_in)
    minMem_low = fuzz.interp_membership(minUniv, min_low, min_in)
    minMem_med = fuzz.interp_membership(minUniv, min_med, min_in)
    minMem_high = fuzz.interp_membership(minUniv, min_high, min_in)
    #minMem_dest = fuzz.interp_membership(minUniv, min_dest, min_in)

    # Calculate rule activation proportions

    # intramembraneous ossification
    activeR3 = np.fmin((nBoneMem_med + nBoneMem_high), np.fmin(maxMem_med, np.fmin(cartMem_low, (distMem_zero + distMem_low + distMem_med))))

    # chondrogenesis
    activeR4 = np.fmin((nBoneMem_high + nBoneMem_med + nCartMem_high + nCartMem_med), np.fmin(boneMem_low, np.fmin((1 - distMem_dest), np.fmin(maxMem_low, (minMem_med)))))

    if currentIteration < 7:
        activeR4 = np.zeros_like(activeR3)

    # cartilage calcification
    activeR5 = np.fmin(boneMem_high, np.fmin(nBoneMem_high, cartMem_low))

    # endochondral ossification
    #activeR6 = np.fmin((nBoneMem_high + nBoneMem_med), np.fmin((1 - cartMem_low), minMem_veryLow))
    activeR6 = np.fmin((nBoneMem_high + nBoneMem_med), np.fmin((1 - cartMem_low), minMem_low))


    # atrophy1
    activeR7 = np.zeros_like(activeR3)
    #activeR7 = np.fmax(dilMem_negDes, dilMem_posDes)

    # atrophy2
    activeR8 = np.zeros_like(activeR3)
    #activeR8 = distMem_dest

    # remodelling
    activeR9 = np.zeros_like(activeR3)
    #activeR9 = np.fmin(dilMem_zero, (distMem_zero + distMem_low))


    activeRules = np.vstack((activeR3, activeR4, activeR5, activeR6, activeR7, activeR8, activeR9))
    changeInBone = ((ruleRatesBone @ activeRules).T).astype(float)
    changeInCart = ((ruleRatesCart @ activeRules).T).astype(float)


    # update smoothingSV
    smoothingSV[:,0,-1] = smoothingSV[:,0,-1] + 0.01 * changeInBone * meshSizeCoefficient
    smoothingSV[:,1,-1] = smoothingSV[:,1,-1] + 0.01 * changeInCart * meshSizeCoefficient
    
    # Enforce 0 to 1 range for state variables
    smoothingSV = np.where(smoothingSV <= 1, smoothingSV, 1)
    smoothingSV = np.where(smoothingSV >= 0, smoothingSV, 0)

    # check that total proportion !>1
    mask = (smoothingSV[:,0,-1] + smoothingSV[:,1,-1]) > 1
    totalProp = np.where(mask, (smoothingSV[:,0,-1] + smoothingSV[:,1,-1]), 1)
    smoothingSV[:,0,-1] = smoothingSV[:,0,-1] / totalProp
    smoothingSV[:,1,-1] = smoothingSV[:,1,-1] / totalProp

    # calculate stateVariables from the updated smoothingSV trailing iterations
    stateVariables = np.mean(smoothingSV, axis=2)

    return stateVariables, activeRules, smoothingSV

