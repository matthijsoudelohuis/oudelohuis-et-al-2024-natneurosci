# oudelohuis-et-al-2024-natneurosci
Code accompanying Oude Lohuis et al. 'Triple dissociation of visual, auditory and motor processing in mouse primary visual cortex'
Published in Nature Neuroscience 2024

Code Description:

The code is organized in directories per figure, and all scripts are named by the corresponding figure panel(s) they produce in the paper. Some code associated with producing data presented in the extended data is found within the folder and under the name of the main figure that the extended data figure relates to. Some code produces figures that did not end up in the publication. Most scripts start with code that loads the raw data from the total set of data (which is not included in this form) and filters the relevant conditions and experiments which is then saved and can be loaded. This part of the code is not functional, but nonetheless included for illustrative purposes. The HelperFuncs directory contains a set of helper functions for plotting the data, statistical tests, etc.

Data description: 

Because of the collection of different experiments and the size of the data, the dataset was split in manageable chunks that were saved separately to facilitate loading the relevant data for different analyses. The following presents the list of separate dataset files that were kept as small as possible for that specific question and analysis: 

Dataset1_1.mat -  Behavioral data of performance across the three task versions:

Dataset1_2.mat -  V1 neural data during task performance in the three cohorts without any cortical manipulations, trial manipulations, etc. (no video data included)

Dataset2_1.mat -  Neural, behavioral, spiking and video data during task performance in the three cohorts. Used for analysis of relationship video motion to V1 sensory tuning.
•	Dataset2_2.mat -  Video motion in an example session
•	Dataset2_3.mat -  Cumulative variance explained by video PCs
•	Dataset3_1.mat -  Filtered neural data in V1 and AC to compare response timing
•	Dataset5_1.mat -  Average spiking data across the session with and without muscimol
•	Dataset5_2.mat -  Behavioral data in sessions with and without muscimol
•	Dataset5_3.mat -  Single neuron sensory responses with and without muscimol
•	Dataset5_4.mat -  Video motion data in sessions with and without muscimol
•	Dataset6_1.mat -  Neural data during conflict trials
•	Dataset6_2.mat -  Behavioral data in conflict trials for reaction times
•	Dataset6_3.mat -  Video motion data during A, V, AV trials
•	Dataset6_4.mat -  Behavioral data with stimulus onset asynchronies in the AV trials
•	Dataset6_5.mat -  Behavioral data for dominance computation in conflict trials
•	DatasetS4_1.mat -  Example session in audiovisual detection task.
•	DatasetS4_2.mat – Behavioral data of animals (n=3) in audiovisual detection task.
•	DatasetS4_3.mat – V1 neural data during audiovisual detection task.
•	DatasetS9_1.mat -  Single session local field potential data in V1 with optogenetic stimulation of auditory cortical fibers to V1.
•	DatasetS9_2.mat -  Single session local field potential data in V1 after checkerboard stimulation.
•	DatasetS9_3.mat – Local field potential data in V1 for all sessions with good LFP signal.
•	GLMfits – fits of generalized linear model with sensory and behavioral data to spiking data, per session during task performance across the three cohorts
•	GLMfits_Musc – GLM fits of sessions with muscimol and control
•	GLMfits_motor – GLM fits with 500 video PCs as predictors to assess V1 firing rate predictions based on video dimensionality


This repository was tested on MATLAB R2016a, MATLAB R2020a with the following toolboxes: Optimization toolbox, Curve fitting toolbox, Statistics and machine learning toolbox, Signal processing toolbox, Image processing toolbox, GUI layout toolbox

How to cite this material: Companion code to Oude Lohuis et al., 2024 “Triple dissociation of visual, auditory and motor processing in mouse primary visual cortex”, XXXXXXXXX (zenodo link), available at https://github.com/matthijsoudelohuis/oudelohuis-et-al-2024-natneurosci/
