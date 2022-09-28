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


<sub>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. If you use the Software for your own research, cite the paper.</sub>

<sub>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.</sub>

Anna Geuzebroek and Simon Kelly, 2022 anna.geuzebroek@ucd.ie / simon.kelly@ucd.ie
