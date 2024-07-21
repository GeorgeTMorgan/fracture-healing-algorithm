import numpy as np

def calcNewMatProps(stateVariables, materials):

    E_wb = 4000 
    E_c = 200
    E_s = 3

    Nu_wb = 0.36
    Nu_c = 0.45
    Nu_s = 0.3 

    materials[:, 0] = E_wb * (stateVariables[:, 0]**3) + E_c * (stateVariables[:, 1]**3) + E_s * ((1 - stateVariables[:, 0] - stateVariables[:, 1])**3)
    materials[:, 1] = Nu_wb * stateVariables[:, 0] + Nu_c * stateVariables[:, 1] + Nu_s * (1 - stateVariables[:, 0] - stateVariables[:, 1])

    return materials

