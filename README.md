# Neurally-Informed-Modelling
General function needed to analysis behavioural and EEG of decision-making tasks (dataAnalysis.m), as well as creating neurally-informed models of the behavioural data (dataModelling.m). 

Code is depending on:
  1) EEGLAB (including Biosig extention).
  2) [CSD toolbox](https://psychophysiology.cpmc.columbia.edu/software/csdtoolbox/) and lay-out.
  3) [findNoisyChannels](see https://github.com/VisLab/EEG-Clean-Tools)
  4) [Brewermap](https://github.com/DrosteEffect/BrewerMap) (to get colors for plot. Can be easily replaced with just choosing colours)
  5) [panels](https://nl.mathworks.com/matlabcentral/fileexchange/20003-panel)
  
% costum-made code:
%   1) dataAnalysis (object to access all the functions, OR get data
%   structure set-up as see below) --> here the data is cut to get the
%   conditions/epochs and will allow and get the EEG data needed for
%   constraining the NI models. 
    
