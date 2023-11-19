# oudelohuis-et-al-2024-natneurosci
Code accompanying Oude Lohuis et al. 'Triple dissociation of visual, auditory and motor processing in mouse primary visual cortex'
Published in Nature Neuroscience 2024

Code Files: The code is organized in directories per main figure. Code associated with producing data presented in the supplementary information is found within the folder of the main figure that the supplementary figure belongs to. All files are named by their corresponding Figures in the paper:

•	Behavior
•	Neural
•	GLM
•	AUC
•	Optogenetics
•	VisuotactileTask
•	Noise correlations

Within each directory the preprocessed data necessary for the reproduction of the figures is present:

•	Data1_1: 
•	Data1_2: 
•	Data2_1: 
•	Data3_1: 
•	Data4_1: 
•	Data4_2: 
•	Data5_1: 
•	Data5_2: 
•	Data5_3: 
•	Data6_1: 
•	Data6_2: 
•	Data6_3:

There are two additional directories:

Utils: set of helper functions

CreateData scripts: original scripts that load the data, filter out the relevant conditions and save the data accordingly. To execute these scripts the total preprocessed dataset is necessary (which is not included). These scripts merely serve to illustrate how the datasets were created.

This repository was tested on MATLAB R2016a, MATLAB R2020a with the following toolboxes: • Optimization toolbox • Curve fitting toolbox • Statistics and machine learning toolbox • Signal processing toolbox • Image processing toolbox • GUI layout toolbox

How to cite this material: Companion code to Oude Lohuis et al., 2022 “Multisensory task demands temporally extend the causal requirement for visual cortex in perception”, doi:10.5281/zenodo.6451623, available at Https://gitlab.com/csnlab/olcese-lab/modid-project/2nd-bump
