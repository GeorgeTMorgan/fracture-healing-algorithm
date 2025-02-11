# fracture-healing-algorithm

Description:
A bone fracture-healing algorithm as part of the journal article “A novel strain-based bone-fracture healing algorithm is able to predict a range of healing outcomes” published in Frontiers of Bioengineering and Biotechnology, 2024. The algorithm works in conjunction with a finite-element model in order to update the tissue type and material properties of fracture-callus finite-elements in response to the strains in the callus region. 

Workflow:
Healing_V1
	The main function (a main() is not used). Input parameters such as the FE model filename and number of nodes which define a neighbouring element. Calls all other files/functions.

preProcessing
	Reads the FE .dat text file to pull-out relevant information such as the element IDs of required element sets, and the element connectivity of the callus mesh. This function is specific for the MSC.Software Marc Mentat FE software, and will need to be rewritten for use with another FE software.

readResults
	Reads the strains and IFMs from the FE model results file. This function is specific for the MSC.Software Marc Mentat FE software, and will need to be rewritten for use with another FE software.

doFuzzyTest
	Implements the biological controller described in the journal article using the sci-kit fuzzy logics toolbox. Updates the stateVariables matrix and stores the biological rule activation rates for post-processing.

updateNeighbors
	Updates the values to be used for ‘nBone’ and ‘nCart’ in the biological controller by searching through the stateVariables matrix and connectivity matrix to determine the maximum bone and cartilage concentration of a neighboring element for each element.

calcNewMatProps
	Calculates each callus finite element’s updated material properties based on the updated tissue concentrations within the element. See published journal article for further details.

writeNewMatProps
	Writes the new callus finite element material properties to the FE model. This function is specific for the MSC.Software Marc Mentat FE software, and will need to be rewritten for use with another FE software.

bendingStiffness
	Creates a 3D model from the axisymmetric FE model, and performs a 4 point bending test on the callus. See published journal article for further details. This function is specific for the MSC.Software Marc Mentat FE software, and will need to be rewritten for use with another FE software.
	

Contact:
George Morgan, george.morgan16@imperial.ac.uk
