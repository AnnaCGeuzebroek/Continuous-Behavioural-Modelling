# Neurally-Informed-Modelling
General function needed to analysis behavioural and EEG of decision-making tasks (dataAnalysis.m), as well as creating neurally-informed models of the behavioural data (dataModelling.m). 

Code is depending on:
  1) [EEGLAB](https://eeglab.org/others/How_to_download_EEGLAB.html) (including Biosig extention).
  2) [CSD toolbox](https://psychophysiology.cpmc.columbia.edu/software/csdtoolbox/) and lay-out.
  3) [findNoisyChannels](https://github.com/VisLab/EEG-Clean-Tools)
  4) [Brewermap](https://github.com/DrosteEffect/BrewerMap) (to get colors for plot. Can be easily replaced with just choosing colours)
  5) [panels](https://nl.mathworks.com/matlabcentral/fileexchange/20003-panel)
  
## mainDataAnalysis
** mainDataAnalysis.m ** gives an tutorial on how to use dataAnalysis.m as object. This includes the following steps:
1) Raw data folder structures.
2) Presetting condition parameters and where to find them in the participants trial files.
3) Presetting preprocessing EEG parameters.
4) Presetting figure lay-outing parameters.
5) Preprocessing and epoching of EEG data and extracting behavioural data.
6) EEG waveform plotting and statistics. 

## mainDataModelling
** mainDataModelling.m ** gives an tutorial on how to use dataModelling.m as object. Not this requires you to use * dataAnalysis.m * to
extract behavioural data and when required neural data to constrain the models. 
