%% Behaviour, EEG data and Eye movement data analysis.
% Continue work on the data set collected in by Craddock, 2016:
% 'Electrophysiological and behavioural indices of decision criterion
% adjustments across context of weak and strong evidence.
% The High-density EEG data are recorded in sampling rate of 512Hz
% using a 128-channel EEG BioSemi in University College Dublin.
%
% Depending on:
%   1) EEGLAB (including Biosig extention).
%   2) dataAnalysis (object to access all the functions).
%   3) CSD layout.
%   5) findNoisyChannels (see zip file)
%   6) Brewermap (to get colors for plot. Can be easily replaced with just choosing colours)


clear
clc
clear global

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------     Input parameters    -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set all paths.
if strcmp(computer, 'PCWIN64')
    
    inputFolder = 'C:\Users\eng\Documents\2. Postdoc_Kelly\Project 5 - Neural correlates of static and dynamic decision-bound adjustments\Project 5.1 - Signatures of bound adjustment\Data';
    addpath(genpath('C:\Users\eng\Documents\MATLAB\dataAnalysis'))
    addpath(genpath('C:\Users\eng\Documents\2. Postdoc_Kelly\Project 5 - Neural correlates of static and dynamic decision-bound adjustments\Project 5.1 - Signatures of bound adjustment\Matlab\Analysis'))

    if ~exist('eeglab', 'file')
        run('C:\Users\eng\Documents\MATLAB\EEGLAB\eeglab'); close all;
    end
else
    % inputFolder = '/Volumes/Keesje/2. Postdoc_Kelly/Project 5.2 - Signatures of bound adjustment/Data';
    % addpath(genpath('/Volumes/Keesje/MATLAB/dataAnalysis'))
end

% some behavioral task that can be check. looking at time on task
% (Within-block effects) or time doing the experiments
blockRandomized = 1; % yes, as you have young and old groups.
numRTBins  = 2; TimeOnTask = 0; TimeOnExperiment = 0;

% Set conditions, here this is left or right tilted as well as inter-trial
% duration. For initial purposes we do not have to look at these, except
% that we bin the reaction times for each conditions to account of possible
% effects.
% IMPORTANT condNames should refer to the way that it was saved in the
% experiment. This can be checked by loading the trial information file and
% see how it will be set in the workspace. Same is happening with
% conditions, e.g. for the left vs. right tilded contrast change, it is
% saved as trialLR as condNames and 1 for left and 2 for right in the so-called
% trialMatrix (matrix that will track the experimental conditions).
conditions{1}  = [1 2]; % easy difficult mixed
condNames{1}   = 'condition';
conditions{2}  = [25 70];
condNames{2}   = 'coherence';

% get your color scheme set here
tmp =  flipud(brewermap(4,'RdBu')); ColoursCond12  = [128 205 193; 1 133 113; 223 194 125; 166 97 26]/255;% [brewermap(2, 'Greens'); tmp([2 1],:)];
ColoursCond3 =  brewermap(4,'Greys'); ColoursCond3 = ColoursCond3-0.16;


% preset the signal processing object.
parameters  = dataAnalysis.getInstance();

% set a couple of standard parameters within this object.
parameters.system = 'bioSemi';   % Object will choice the approtiated loading function of EEGLAB
parameters.analysisEEG      = 1; % Analyse EEG (obviously we want to look at this)
parameters.analysisEyelink  = 0;
parameters.analysisEOG      = 1; % Use additionally EOG electrodes and Frontal electrodes
                                 % to look for blinks and
                                 % exclude trials with blinks (We do this to prevent
                                 % any possible effects of the blink on the sensory perception)
parameters.analysisEMG = 0;
parameters.numTrials   = 24; % number of trials per block.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------     GET EXPERIMENT INFORMATION    ------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function looping through the data files of each participants. Identifying
% the behavioural information. including numBlocks, numBlocks per
% condition, numStimulus in each conditions and their specific condition (
% here initialy all stimulus in each block are 1 condition.) However, it is
% possible that the conditions information are in the EEG files. So here
% simply put all the par and other parameters of each participant away.
parameters.getInformation(0, inputFolder, conditions, condNames, blockRandomized);
clear cond* *Folder


% recode the experiment.trialCond from 1:4 as 1:2 is no gap, 4:6 is gap. As
% well as 1,2,3 being [2 4 6] ITI respectively.
parameters.conditions{3}   = [2 4 6 8];
parameters.condNames{3}    = 'ITI';

parameters.stim.lengthITI  = [2 4 6 8];
parameters.stim.namesITI   = 'ITI';

% NOTE below was needed to created the right trialMatrix construction, this
% has already been preformed in the epoched data
%{
for indPP = 1:length(parameters.ppNames)
    for indBlock = 1:parameters.numBlocks(indPP)
        
        tmpITI = zeros(size(parameters.experiment{indPP}{indBlock}.trialCond));

        if length(parameters.experiment{indPP}{indBlock}.par.correctCoh) == 1
            parameters.experiment{indPP}{indBlock}.coherence  = parameters.experiment{indPP}{indBlock}.par.correctCoh{:};
            parameters.experiment{indPP}{indBlock}.condition  = 1;%find(parameters.conditions{2} == parameters.experiment{indPP}{indBlock}.par.correctCoh{:});
            
            for indITI = 1:length(parameters.conditions{3})
                tmpITI(parameters.experiment{indPP}{indBlock}.trialCond == indITI) = parameters.conditions{3}(indITI);
            end
        else
            tmpCoh = zeros(size(parameters.experiment{indPP}{indBlock}.trialCond));
            tmpCoh(ismember(parameters.experiment{indPP}{indBlock}.trialCond,1:2:8)) = 70;
            tmpCoh(ismember(parameters.experiment{indPP}{indBlock}.trialCond,2:2:8)) = 25;
            
            parameters.experiment{indPP}{indBlock}.coherence  = tmpCoh;
            parameters.experiment{indPP}{indBlock}.condition  = 2;
            indexITI = 1:2:length(parameters.conditions{3})*2;
            for indITI = 1:length(indexITI)
                tmpITI(parameters.experiment{indPP}{indBlock}.trialCond == indexITI(indITI))  = parameters.conditions{3}(indITI);
                tmpITI(parameters.experiment{indPP}{indBlock}.trialCond == indexITI(indITI)+1) = parameters.conditions{3}(indITI);
            end
        end
        parameters.experiment{indPP}{indBlock}.ITI = tmpITI;
    end
end
%}

% assuming that parameters are the same for all participants and with all
% blocks
parameters.DetectOrDisc = 0;

parameters.stim.refreshRate  = 60;   % parameters.experiment{1}{1}.par.videoFrate;	 % Screens refreshRate
parameters.stim.duration     = 1000; % parameters.experiment{1}{1}.par.targetDur.*1000;    % duration of the stimulus
parameters.stim.freqSSVEP    = 15;
parameters.stim.epoch        = [-1800 2500];

parameters.stim.RTdeadLine = [-1.8 1.65];
parameters.stim.RTCutOff   = 0.25; % only include reaction times after the gap offset.
parameters.stim.FACutOff   = -1.5;
parameters.stim.FA         = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------    EEG PROCESSING PARAMETERS    -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All parameters needed to set-up the pre-processing. We use a most
% simplest pre-processing protocol with the following steps:
%   1) detrending for all (1) or only for the channels that it would
%   benefit (2).
%   2) Low-pass filtered using a filter designed by Simon. (see obj.setFilters)
%   3) High-pass filtered using filtfilt.
%   4) Additionally, the findNoisyChannels within AutoMagic.

%---------------     set pre-processing parameters   ----------------------
parameters.eeg.NumberOfChannels = 128;
% load in a structure 'chanlocs' containing the x,y,z locations of each of the 128 scalp channels in the cap.
load chanlocsBioSemi128;
elab = strings(parameters.eeg.NumberOfChannels,1);
for indChan = 1:parameters.eeg.NumberOfChannels
    elab(indChan) = chanlocs(indChan).labels;
end
elab(end+1)='FCz'; clear ind*

parameters.eeg.ChannelsName   	= elab;
parameters.eeg.SampleRate   	= 512;	% EEG sample rate

parameters.eeg.chanlocs         = chanlocs;	% electrode locations
parameters.eeg.transChanlocs	= load('CSD_coords.mat');	% electrode locations

parameters.eeg.ampSaturation	= 3276;	% What's the maximum limiting amplitude of the amplifiers, to detect saturation?

parameters.eeg.postDCC  = 700; % how many samples after DCC is EEG rubbish?
parameters.eeg.applyCSD = 1;

parameters.eeg.applyFindDCcorrect   = 0;

parameters.eeg.badChannel   = 1;	% set on one if you want to apply bad channel selection and interpolation before epoching.
parameters.eeg.applyLPF 	= 1;    % attention here that the lowpass filter is different in New York as in Ireland.
parameters.eeg.cuttoffLPF	= 38;	% 38 cut-off threshold for low-pass filter
% low-pass filter will be different as NY has 60 hz electricity rate.
parameters.eeg.LHamming     = 77;   % 77
parameters.eeg.applyHPF     = 0;    % look at raw data if is nescessary.
parameters.eeg.simonHPF     = 0;
parameters.eeg.HPFcutoff    = 0;
parameters.eeg.applyDetrend	= 1;    % option 1(all)/ 2(only that will benefit from it)

parameters.eeg.artifactThres = 80;
parameters.eeg.artifactEOG   = 200;
parameters.eeg.intArtifact   = 1;

parameters.eeg.channelVEOG   = 1:2;
parameters.eeg.channelHEOG   = 3:4;

parameters.eeg.timing = 'past';

%--------------- set triggers and EPOCHS definition -----------------------
% The trigger codes we used for the different events were as follows:
parameters.triggers.start       = 1;            % Start of the trial
parameters.triggers.stimulusOFF = 5;            % Target off
parameters.triggers.stimulusON  = [170 125];    % 170 for 70% and 125 for 25%
parameters.triggers.response    = [12 13];      % left vs. right

% EPOCH DEFINITION
parameters.eeg.epochLock = parameters.triggers.stimulusON;
parameters.eeg.baseline  = 0 + [-1 1]*1000/parameters.stim.freqSSVEP;  

% Initially we set a rather larger epoch based on the stimulus definitition
% 'parameters.stim.epoch'. This is used to cut out a larger epoch to apply
% all the ERP pre-processing on. Afterwards smaller epochs are created to
% get the target- and response-locked epochs. This is usefull especially
% for the ERD (event related desynchronization) and the SSVEP later
% onwards.

% these range are used for the initial plots to explore the full range.
% We later zoomed in to the new ranges to get a better view. 
parameters.eeg.targetEpoch   = -ceil(500/(1000/parameters.eeg.SampleRate))*(1000/parameters.eeg.SampleRate):1000/parameters.eeg.SampleRate:ceil(parameters.stim.RTdeadLine(2)*1000/(1000/parameters.eeg.SampleRate))*(1000/parameters.eeg.SampleRate);
parameters.eeg.responseEpoch = -ceil(800/(1000/parameters.eeg.SampleRate))*(1000/parameters.eeg.SampleRate):1000/parameters.eeg.SampleRate:ceil(400/(1000/parameters.eeg.SampleRate))*(1000/parameters.eeg.SampleRate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------    Perform eeg prep  -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE usually this needs to be run to get the preprocessed data in order
% to actually created the epoched data. 
% parameters.applyPreprocessing;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- Reset behaviour data    -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE below was needed to created the right trialMatrix construction, this
% has already been preformed in the epoched data
%{
for indPP = 1:length(parameters.ppNames)
    for indBlock = 1:parameters.numBlocks(indPP)
        Responses = find(ismember(parameters.event{indPP}{indBlock}.Number, parameters.triggers.response));
        StimOn    = find(ismember(parameters.event{indPP}{indBlock}.Number, parameters.triggers.stimulusON));
        if length(StimOn) > parameters.numTrials
            parameters.event{indPP}{indBlock}.Number(StimOn(end):end) = [];
            parameters.event{indPP}{indBlock}.Times(StimOn(end):end) = [];
            
            StimOn(end) = [];
            Responses(Responses > StimOn(end)) = [];
        end
        
        parameters.experiment{indPP}{indBlock}.RespT   = parameters.event{indPP}{indBlock}.Times(Responses)*(1/parameters.eeg.SampleRate);
        parameters.experiment{indPP}{indBlock}.TargOnT = parameters.event{indPP}{indBlock}.Times(StimOn)*(1/parameters.eeg.SampleRate);
        parameters.experiment{indPP}{indBlock}.RespLR = ones(size(parameters.experiment{indPP}{indBlock}.RespT));
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------   Epoching data    -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters.epochingData;
% parameters.applyCSD; % NOTE already applied in the epoched EEG data
keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------   Behavioural data    -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we create equal sized Reaction Time (RT) bins per participant per conditions
% to check the alignment of the CPP with the actual median RT as the CPP
% peak should be highly linked with the reaction time.
parameters.binRTs(1, [1/3 2/3], [1 2]); % as there are sometimes to little trials per all conditions, we need to resort to a smaller number. 

% Include the context as such.. although we aren't really using it anymore.
parameters.conditions{4} = [1 2 3];

for indPP = 1:length(parameters.ppNames)
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 2, 4) = 3;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 1 & parameters.behaviour{indPP}.trialMatrix(:,2) == 70, 4) = 2;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 1 & parameters.behaviour{indPP}.trialMatrix(:,2) == 25, 4) = 1;
end


%% %%%%%%%%%%%%%%%%%%%% Pre-set figure parameters    %%%%%%%%%%%%%%%%%%%%%%

parameters.figLayOut.letterSize  = 9;
parameters.figLayOut.letterType  = 'Arial';
parameters.figLayOut.lineWidth   = 1.2;
parameters.figLayOut.lineType    = {'-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-'};
parameters.figLayOut.legNames{1} = {'Constant', 'Mixed'};
parameters.figLayOut.legNames{2} = {'Weak', 'Strong'}; % shift as they are 'order' by strength
parameters.figLayOut.legNames{3} = {'2', '4', '6', '8'};
parameters.figLayOut.legTitle{1} = 'Context';
parameters.figLayOut.legTitle{2} = 'Evidence strength';
parameters.figLayOut.legTitle{3} = 'ITI duration';

parameters.figLayOut.legNames{4}  = {'Weak', 'Strong', 'Mixed'};
parameters.figLayOut.legTitle{4} = 'Context';

parameters.figLayOut.plotCI      = 0.05; % get shaded area of not.
parameters.figLayOut.removeInter = 1;
parameters.figLayOut.saveDim     = [5 11];
parameters.figLayOut.colours     = ColoursCond12;

% these range are used for the initial plots to explore the full range.
% We later zoomed in to the new ranges to get a better view. 
parameters.figLayOut.targetLim   = [-200 0:500:1000];
parameters.figLayOut.responseLim = [-400 0 200];
parameters.figLayOut.plotRT      = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------     Waveform plotting       -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now transfer parameters to EEG object. These can be run without each
% other.
close all

parameters.eeg.baseline  = 0 + [-1 0]*1000/parameters.stim.freqSSVEP;
parameters.figLayOut.colours = ColoursCond12;
parameters.figLayOut.legends = {'Weak Context', 'Strong Context', 'Mixed (Weak)', 'Mixed (Strong)'};

%% ------------------- get CPP -----------------------------------------
% define topoplot information
trangeTopo  =  0 + [-2 -1]*1000/parameters.stim.freqSSVEP; % note if you add a second row, these will be used for FA plots. Only nesc. when you are using target-locked topos
TargetOrResponse = 2;       % as written 1) target-locked 2) response-locked
AverageOrSlope   = 1;       % as written 1) for grand average 2) slope

% define rest
baselineCorrect  = 1;
methodUsed       = 2;       % Method used can be 1) choose spec. electrodes, 2) choose cluster of electrodes and get the 3 best based on SNR, 3) apply lucalize
plotThis         = 2;       % 1) individual plots, 
                            % 2) average Hits and seperated average Misses and FA
                            % 3) include all average target-locked Misses
                            % 4) plot target and response seperately.
                            % 5) get first derivative
grouping         = 0;       % if you want to include seperated plots for a specific condition (NOTE THIS should be your first in the plotCondition)
plotCondition    = [1 2];
%
% Plot target and response locked CPP (target-locked including all misses)
% as a function of target context as well as target evidence strength.
%{
[~, fig]  = parameters.plotERP(methodUsed, 'CPP',  'A5 A19 A32 A18 A20 A31 A17 A21 A30' , baselineCorrect,...
   plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse, AverageOrSlope);
keyboard

figure(fig.Handle{1});
fig.Info{1}(1,2).select()
legend('off')
ylim([-3 35])
yticks([0:10:30])

getAxes = gca;
shaded = shadedErrorBar(trangeTopo, repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle  = 'none'; shaded.edge(2).LineStyle  = 'none';
shaded.patch.FaceAlpha    = 0.1;
fig.Info{1}(1,1).select()
yticks([0:10:30])

plotSave(gca, 'CPP_Shaded.png', fullfile(parameters.figFolder, 'groupAverage\HPF0_CSD1\CPP\'), parameters.figLayOut.saveDim)
%}

% plot to compare with the sim. CPP
%{
keyboard

saveHere = 'C:\Users\eng\Documents\2. Postdoc_Kelly\Project 5 - Neural correlates of static and dynamic decision-bound adjustments\Project 5.1 - Signatures of bound adjustment\Data\Figures\groupAverage\HPF0_CSD1\CPP\';

parameters.figLayOut.saveDim = [4 6];
parameters.figLayOut.targetLim   = [0 500:500:1000];
plotThis = 5;
parameters.figLayOut.plotCI = 0;

[~, fig]  = parameters.plotERP(methodUsed, 'CPP',  'A5 A19 A32 A18 A20 A31 A17 A21 A30' , baselineCorrect,...
    plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse, AverageOrSlope);

figure(fig.Handle{1})
ylim([-1 30])
yticks([0:10:30])
legend('off')
ylabel('')
plotSave(gca, 'CPPModel.png', saveHere, parameters.figLayOut.saveDim);

% plot to compare with the sim. CPPP
parameters.figLayOut.targetLim   = [0 500:500:1000];
plotThis = 5;
parameters.figLayOut.plotCI = 0;

[~, fig]  = parameters.plotERP(methodUsed, 'CPP',  'A5 A19 A32 A18 A20 A31 A17 A21 A30' , baselineCorrect,...
    plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse, AverageOrSlope);

figure(fig.Handle{1})
ylim([-0.2 0.3])
yticks([-0.2:0.2:.2])
legend('off')
plotSave(gca, 'derCPPModel.png', saveHere, [parameters.figLayOut.saveDim(1) 6.27]);

%}

% plot Difference topo plots
%{
% parameters.figLayOut.plotCI    = 0.2;
% parameters.figLayOut.targetLim  = -200:200:1000;
% parameters.figLayOut.saveDim    = [5 11];
% parameters.plotDifferenceTopo('CPP', 'A5 A19 A32 A18 A20 A31 A17 A21 A30', [1 2], trangeTopo, TargetOrResponse, AverageOrSlope);
%}

% get waveform parameters.
%{

TargetOrResponse = 2;
trangeTopo(1,:)  = 0 + [-2 -1]*1000/parameters.stim.freqSSVEP;
trangeTopo(2,:)  = [-250 -1*1000/parameters.stim.freqSSVEP];
PeakMeanOrMax    = 1;
plotCondition    = [1 2];
negOrPos         = 2;

parameters.getWaveformParameters(methodUsed, 'CPP', [1/3 2/3], baselineCorrect,...
    1, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos, PeakMeanOrMax);

tblCPP = parameters.getWaveformParameters(methodUsed, 'CPP', [ ], baselineCorrect,...
    1, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos, PeakMeanOrMax);


lmeCPP = fitglme(tblCPP,  'Peak ~  1 +  Context*Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');

parameters.plotDifferenceTopo('CPP', [], [1 2], trangeTopo(1,:), TargetOrResponse, AverageOrSlope, 0, 1, [15:30]);

%}

% ITI duration
%{
baselineCorrect = 0;
parameters.figLayOut.colours = ColoursCond12([1 1 1 1 2 2 2 2 3 3 3 3],:);
parameters.figLayOut.lineType = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.', ':'};

% ITI plotting.
plotCondition = [4 3]; grouping = 4;
parameters.plotITI(methodUsed, 'CPP',  'A5 A19 A32 A18 A20 A31 A17 A21 A30' , 1,...
    plotCondition, grouping);

parameters.figLayOut.colours = repmat(ColoursCond12(1:3,:),4,1);
parameters.figLayOut.lineType = repmat({'-'},12, 1)';


% ITI plotting.
plotCondition = [4]; grouping = 0;
parameters.plotITI( methodUsed, 'CPP',  'A5 A19 A32 A18 A20 A31 A17 A21 A30' , 1,...
    plotCondition, grouping);
%}


%% ------------------- get Motor preperation ---------------------------
%
trangeTopo  = 0 + [-2 -1]*1000/parameters.stim.freqSSVEP;
AverageOrSlope   = 1;
TargetOrResponse = 2;

baselineCorrect  = 2;
parameters.eeg.baseline =  0 + [-2 0]*1000/parameters.stim.freqSSVEP; % Motor Execution threshold
parameters.figLayOut.targetLim   = [-200 500:500:1000];
parameters.figLayOut.plotCI = 0.05;
parameters.figLayOut.saveDim     = [5 11];

methodUsed       = 1;
grouping         = 0;
plotCondition    = [1 2];
plotThis         = 2; % 1) plot individual, 2) plot only averages.

% use specific code as the people werent instructed to use the right hand. 
% handmess(parameters, 'Beta ERD', [], trangeTopo, TargetOrResponse, AverageOrSlope, 1, [13:30]);
possElec{1} = 'D20 D19 D18 D27 D28 D17'; % right handed response
possElec{2} = 'B21 B22 B23 B19 B18 B17'; % left handed response
Elec{1} = possElec{1};
Elec{2} = possElec{1};
Elec{3} = possElec{1};
Elec{4} = possElec{1};
Elec{5} = possElec{2};
Elec{6} = possElec{1};
Elec{7} = possElec{1};
Elec{8} = possElec{1};
Elec{9} = possElec{1};
Elec{10} = possElec{1};
Elec{11} = possElec{1};
Elec{12} = possElec{2};
Elec{13} = possElec{1};
Elec{14} = possElec{1};

%{
keyboard
%{
[Quality, fig] = parameters.plotERP(methodUsed, 'Beta ERD', Elec,...
    baselineCorrect, plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse, AverageOrSlope, 2, 15:30);

figure(fig.Handle{1});
fig.Info{1}(1,2).select()

legend('off')
ylim([-0.2 0.8])
yticks([-0.2:0.2:0.8])
line([0 0], [-0.1 0.8], 'Color', 'k', 'LineWidth', 2)
% set(gca, 'YDir','reverse')

fig.Info{1}(1,1).select()
yticks([-0.2:0.2:0.8])
line([0 0], [-0.1 0.8], 'Color', 'k', 'LineWidth', 2)
getAxes = gca;
shaded = shadedErrorBar(trangeTopo, repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle  = 'none'; shaded.edge(2).LineStyle  = 'none';
shaded.patch.FaceAlpha    = 0.1;
% set(gca, 'YDir','reverse')

plotSave(gca, 'Beta ERD.png', fullfile(parameters.figFolder, 'groupAverage\HPF0_CSD1\Beta ERD\'), parameters.figLayOut.saveDim)
%}

% First plot Beta without baseline correction
%
[Quality, fig] = parameters.plotERP(methodUsed, 'Beta ERD', Elec,...
    0, plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse, AverageOrSlope, 2, 15:30);

parameters.figLayOut.colours = [ColoursCond12([1 1 1 2 2 2],:); ColoursCond12([3  3 3  4 4 4],:)];
parameters.figLayOut.lineType = {'-' '-.' ':' '-' '-.' ':' '-' '-.'  ':' '-' '-.' ':'}


[Quality, fig] = parameters.plotERP(methodUsed, 'Beta ERD', Elec,...
    0, 3, [1 2 0], 1, trangeTopo, TargetOrResponse, AverageOrSlope, 2, 15:30);

figure(fig.Handle{1});
fig.Info{1}(1,2).select()

legend('off')
ylim([3.5 5])
yticks([3.5:0.5:5])
line([0 0], [3.5 5], 'Color', 'k', 'LineWidth', 2)


figure(fig.Handle{2});
fig.Info{2}(1,2).select()

legend('off')
ylim([3.5 5])
yticks([3.5:0.5:5])
line([0 0], [3.5 5], 'Color', 'k', 'LineWidth', 2)

% set(gca, 'YDir','reverse')
getAxes = gca;
shaded = shadedErrorBar(0 + [-2 -1]*1000/parameters.stim.freqSSVEP, repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle  = 'none'; shaded.edge(2).LineStyle  = 'none';
shaded.patch.FaceAlpha    = 0.1;

fig.Info{1}(1,1).select()
ylim([3.5 5])
yticks([3.5:0.5:5])
line([0 0], [3.5 5], 'Color', 'k', 'LineWidth', 2)getAxes = gca;
shaded = shadedErrorBar( 0 + [-2 0]*1000/parameters.stim.freqSSVEP, repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle  = 'none'; shaded.edge(2).LineStyle  = 'none';
shaded.patch.FaceAlpha    = 0.1;
% set(gca, 'YDir','reverse')

plotSave(gca, 'Beta ERD_noBaseline.png', fullfile(parameters.figFolder, 'groupAverage\HPF0_CSD1\Beta ERD\'), parameters.figLayOut.saveDim)
%}

% parameters.plotDifferenceTopo('Beta ERD', [], [1 2], trangeTopo, TargetOrResponse, AverageOrSlope,1, [15:30]);

% ITI period
%{
keyboard

baselineCorrect = 3; %ColoursCond12;%
parameters.figLayOut.colours  = ColoursCond12([1 2 4],:); %ColoursCond12([1 1 1 1 2 2 2 2 4 4 4 4],:); %
parameters.figLayOut.lineType = repmat({'--'},4,1) ;%
parameters.figLayOut.saveDim 
for indPP = 1:length(parameters.ppNames)
    parameters.behaviour{indPP}.trialMatrix(:,5) = nan;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,3) == 8,5) = 1;
end
parameters.conditions{5} = 1;
parameters.figLayOut.legNames{5} = '8 sec';
parameters.figLayOut.legTitle{5} = 'ITI';

% ITI plotting.
plotCondition = [4 5]; grouping =0 ;
parameters.figLayOut.plotCI = 0;
fig = parameters.plotITI(methodUsed,  'Beta ERD', Elec, baselineCorrect,...
    plotCondition, grouping, -67 + [-1 0]*1000/parameters.stim.freqSSVEP);

getAxis = gca;
%getAxis.YAxis.Visible = 'off'; 
xlim([-8000 500])
ylim([0 0.5])
yticks([0:0.25:0.5])

xlabel (' ')
xticks([-8000:2000:0])
xticklabels([0:2:8])
%line([-8000 -8000], [180 0], 'Color', 'black', 'LineWidth', 2)
line( [-8000 500],[0 0], 'Color', 'black', 'LineWidth', 2)
legend('off')

% legend(fig.legendThis, {'Weak', 'Strong',  'Mixed'}, 'Location', 'northwest')
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);

plotSave(gca,  ['Beta ERD_ITI.png'], hereFigFolder, [4 6]);% 

fig.legendThis.delete

ylim([50 135])
getAxis.YAxis.Visible = 'off'; 
getAxis.XAxis.Visible = 'off'; 
xticks([-8000:2000:0])
xlim([-8000 1000])
xlabel('off')

plotSave(gca,  ['Beta ERD_ITI_model.png'], hereFigFolder, [1.4 2.1].*1.5);% 

keyboard
%}


%% Plot BETA as marker of motor preperation
%
keyboard
hereFigFolder = fullfile(parameters.figFolder,  'groupAverage',...
    ['HPF' num2str(parameters.eeg.HPFcutoff) '_CSD' num2str(parameters.eeg.applyCSD)],...
    'Beta ERD', 'RTBINS2/');

baselineCorrect  = 2;
TargetOrResponse = 1;
negOrPos = 1;
PeakMeanOrMax = 1;
trangeTopo(2,:)  =  [-500 0]; 

% 1) excursion
%
trangeTopo(1,:)  =  0 + [-2 -1]*1000/parameters.stim.freqSSVEP; %

parameters.getWaveformParameters(methodUsed, 'Beta ERD', [1/3 2/3], baselineCorrect,...
    1, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);

xlim([300 1050])
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'amplitudeExcursion.png', hereFigFolder, [4.6 6.4]);

tbl.Excursion = parameters.getWaveformParameters(methodUsed, 'Beta ERD', [], baselineCorrect,...
   1, [plotCondition], grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);

parameters.plotDifferenceTopo('Beta ERD', [], [1 2], trangeTopo(1,:), TargetOrResponse,...
    AverageOrSlope, 0, 1, [15:30], 0, 1);


% plot for ITI effects
%{
tbl.Excursion = parameters.getWaveformParameters(methodUsed, 'Beta ERD', [1/3 2/3], baselineCorrect,...
  1, [4 3], grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);

tbl.Excursion.Context = reordercats(tbl.Excursion.Context, [2 1 3]);

lme.ITI = fitglme(tbl.Excursion,  'Peak ~  1 + ITIduration*Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');
%}

ylim([0 0.8])
yticks(0:0.2:1)
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'BetaERD_exursion.png', hereFigFolder, [4.17 2.7]);

% check for RT interaction effects:
lme.RTinteraction = fitglme(tbl.Excursion,  'Peak ~  1 + Context*Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');

lme.Excursion = fitglme(tbl.Excursion,  'Peak ~  1 + Context*Evidencestrength + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.FixedExcursion = tbl.Excursion(double(tbl.Excursion.Context) == 1,:);
lme.FixedExcursion  = fitglme(tbl.FixedExcursion ,  'Peak ~  1 + Evidencestrength + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.MixedExcursion = tbl.Excursion(double(tbl.Excursion.Context) == 2,:);
lme.MixedExcursion = fitglme(tbl.MixedExcursion,  'Peak ~  1 + Evidencestrength + (1|ppNames)',...
    'FitMethod','Laplace');   

tbl.WeakExcursion = tbl.Excursion(double(tbl.Excursion.Evidencestrength) == 1,:);
lme.WeakExcursion  = fitglme(tbl.WeakExcursion ,  'Peak ~  1 + Context + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.StrongExcursion = tbl.Excursion(double(tbl.Excursion.Evidencestrength) == 2,:);
lme.StrongExcursion = fitglme(tbl.StrongExcursion,  'Peak ~  1 + Context + (1|ppNames)',...
    'FitMethod','Laplace'); 

diary(fullfile(parameters.logFolder, 'WaveformsStats', 'lmeERD.txt'))
clc
fprintf('Waveform stats checking for excursion effects.\nInitial assumption that they are context and RT-invariant by baseline correction to motor execution threshold.\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
lme.RTinteraction
fprintf('While there is an interaction effect (e.g. smaller excursion relates to faster RTs,\nthere are no interaction effects\n');
fprintf('Remove RT effects\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
lme.Excursion
fprintf('Interaction effect of Context and Motion Coherence\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
fprintf('Posthocs of Context effects show that there only is a difference between Weak and Strong in fixed condition=\n');
Fixed = lme.FixedExcursion
fprintf('And not in the mixed condition\n');
Mixed = lme.MixedExcursion
fprintf('--------------------------------------------------------------------------------------------------------------\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
fprintf('Posthocs of Motion Coherence effects show that there only is a difference between Strong in fixed and mixed condition=\n');
Fixed = lme.StrongExcursion
fprintf('but not for Weak\n');
Mixed = lme.WeakExcursion
fprintf('--------------------------------------------------------------------------------------------------------------\n');
diary off
%}
%
keyboard
% 2) baseline effects
%
PeakMeanOrMax    = 1;
plotCondition    = [1 2];
negOrPos = 1;
baselineCorrect  = 0;
TargetOrResponse = 1;
trangeTopo(1,:)  =  0 + [-4 -1]*1000/parameters.stim.freqSSVEP; %

%
% parameters.getWaveformParameters(methodUsed, 'Beta ERD', [1/3 2/3], baselineCorrect,...
%     1, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);
% 
% ylim([3.4 5])
% xlim([300 1050])
% xticks(400:200:1100);
% xticklabels((400:200:1100)./1000)
% set(gca,'FontSize', parameters.figLayOut.letterSize);
% set(gca,'FontName', parameters.figLayOut.letterType);
% % set(gca, 'YDir','reverse')
% 
% plotSave(gca, ['amplitudeBaseline.png'], hereFigFolder, [3.7 4.5]);%[4.3 obj.figLayOut.saveDim(2)/3]);

for indPP = 1:length(parameters.ppNames)  
    parameters.behaviour{indPP}.trialMatrix(:,5) = parameters.behaviour{indPP}.binRT;
    parameters.behaviour{indPP}.trialMatrix((parameters.behaviour{indPP}.trialMatrix(:,5) == 3),5) = nan;
    parameters.behaviour{indPP}.trialMatrix(:,6) = NaN;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,2) == 70,6) = ...
      parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,2) == 70,1);
end

parameters.conditions{5} = 'speed'; parameters.figLayOut.legTitle(5) = {'RTspeed'}; parameters.figLayOut.legNames{5} = {'fast', 'slow'};
parameters.conditions{6} = 'strong'; parameters.figLayOut.legTitle(6) = {'Strong'}; parameters.figLayOut.legNames{6} = {'fixed', 'mixed'};

parameters.plotDifferenceTopo('Beta ERD_baseline', [], [2 1], trangeTopo(1,:),...
    TargetOrResponse, AverageOrSlope, 0, 1, [15:30], 0, 0);


parameters.plotDifferenceTopo('Beta ERD_baseline', [], [5 6], trangeTopo(1,:),...
    TargetOrResponse, AverageOrSlope, 0, 1, [15:30], 0,0);

tbl.BetaBaseline = parameters.getWaveformParameters(methodUsed, 'Beta ERD', [], baselineCorrect,...
    0, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);


lme.RTinteraction = fitglme(tbl.BetaBaseline,  'Peak ~  1 + Context*Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');

% check context
tbl.Weak = tbl.BetaBaseline(double(tbl.BetaBaseline.Evidencestrength) == 1,:);
lme.Weak = fitglme(tbl.Weak,  'Peak ~  1 + Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong = tbl.BetaBaseline(double(tbl.BetaBaseline.Evidencestrength) == 2,:);
lme.Strong = fitglme(tbl.Strong,  'Peak ~  1 + Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong = tbl.BetaBaseline(double(tbl.BetaBaseline.Evidencestrength) == 2 & double(tbl.BetaBaseline.Context) == 1,:);
lme.Strong1 = fitglme(tbl.Strong,  'Peak ~  1 + RT + (1|ppNames)',...
    'FitMethod','Laplace');
tbl.Strong= tbl.BetaBaseline(double(tbl.BetaBaseline.Evidencestrength) == 2 & double(tbl.BetaBaseline.Context) == 2,:);
lme.Strong2 = fitglme(tbl.Strong,  'Peak ~  1 + RT + (1|ppNames)',...
    'FitMethod','Laplace');

diary(fullfile(parameters.logFolder, 'WaveformsStats', 'lmeERD_baseline.txt'))
clc
fprintf('Waveform stats checking for baseline effects.\nBaseline activity shows a signifant interaction effect bewteen RT*context*evidencestrength\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
lme.RTinteraction
fprintf('While lawfully lower beta synchronization causes faster reaction times in all conditions but the Strong Context\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
fprintf('--------------------------------- WEAK -------------------------------------------\n');
lme.Weak
fprintf('--------------------------------- STRONG -------------------------------------------\n');
lme.Strong
fprintf('--------------------------------- interaction effect -------------------------------------------\n');
lme.Strong1
lme.Strong2
fprintf('--------------------------------------------------------------------------------------------------------------\n');
diary off
%}


% 2) Pre response effects
%
baselineCorrect  = 0;
TargetOrResponse = 2;
PeakMeanOrMax    = 2;

trangeTopo(1,:)  =  0 + [-2 0]*1000/parameters.stim.freqSSVEP; %
trangeTopo(2,:)  =  0 + [-500 0]*1000/parameters.stim.freqSSVEP; %

% parameters.getWaveformParameters(methodUsed, 'Beta ERD', [1/3 2/3], baselineCorrect,...
%     1, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);
% 
% ylim([3.4 5])
% xlim([300 1050])
% xticks(400:200:1100);
% xticklabels((400:200:1100)./1000)
% 
% set(gca,'FontSize', parameters.figLayOut.letterSize);
% set(gca,'FontName', parameters.figLayOut.letterType);
% % set(gca, 'YDir','reverse')
% 
% plotSave(gca, ['amplitudePreResponse.png'], hereFigFolder, [3.7 4.5]);%[4.3 obj.figLayOut.saveDim(2)/3]);

tbl.BetaPreresponse = parameters.getWaveformParameters(methodUsed, 'Beta ERD', [], baselineCorrect,...
    0, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);

parameters.plotDifferenceTopo('Beta ERD_Preresponse', [], [2 1],...
    0 + [-2 2]*1000/parameters.stim.freqSSVEP, TargetOrResponse,...
    AverageOrSlope, 0, 1, [15:30],PeakMeanOrMax,0);

parameters.plotDifferenceTopo('Beta ERD_Preresponse', [], [6 5],...
    0 + [-2 2]*1000/parameters.stim.freqSSVEP, TargetOrResponse,...
    AverageOrSlope, 0, 1, [15:30],PeakMeanOrMax,0);
% 

lme.RTinteraction = fitglme(tbl.BetaPreresponse,  'Peak ~  1 + Context*Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');

% check context
tbl.Weak = tbl.BetaPreresponse(double(tbl.BetaPreresponse.Evidencestrength) == 1,:);
lme.Weak = fitglme(tbl.Weak,  'Peak ~  1 + Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong = tbl.BetaPreresponse(double(tbl.BetaPreresponse.Evidencestrength) == 2,:);
lme.Strong = fitglme(tbl.Strong,  'Peak ~  1 + Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong  = tbl.BetaPreresponse(double(tbl.BetaPreresponse.Evidencestrength) == 2 & double(tbl.BetaPreresponse.Context) == 1,:);
lme.Strong1 = fitglme(tbl.Strong,  'Peak ~  1 + RT + (1|ppNames)',...
    'FitMethod','Laplace');
tbl.Strong  = tbl.BetaPreresponse(double(tbl.BetaPreresponse.Evidencestrength) == 2 & double(tbl.BetaPreresponse.Context) == 2,:);
lme.Strong2 = fitglme(tbl.Strong,  'Peak ~  1 + RT + (1|ppNames)',...
    'FitMethod','Laplace');
% 
% diary(fullfile(parameters.logFolder, 'WaveformsStats', 'lmeERD_preResponse.txt'))
% clc
% fprintf('Waveform stats checking for baseline effects.\nBaseline activity shows a signifant interaction effect bewteen RT*context*evidencestrength\n');
% fprintf('--------------------------------------------------------------------------------------------------------------\n');
% lme.RTinteraction
% fprintf('While lawfully lower beta synchronization causes faster reaction times in all conditions but the Strong Context\n');
% fprintf('--------------------------------------------------------------------------------------------------------------\n');
% fprintf('--------------------------------- WEAK -------------------------------------------\n');
% lme.Weak
% fprintf('--------------------------------- STRONG -------------------------------------------\n');
% lme.Strong
% fprintf('--------------------------------- interaction effect -------------------------------------------\n');
% lme.Strong1
% lme.Strong2
% fprintf('--------------------------------------------------------------------------------------------------------------\n');
% diary off
% %}
%{
all = tbl.BetaPreresponse;
all.Baseline = tbl.BetaBaseline.Peak;

tbl.Weak  = all(double(tbl.BetaPreresponse.Evidencestrength) == 1 & double(tbl.BetaPreresponse.Context) == 1,:);
lme.Weak1 = fitglme(tbl.Weak,  'Peak ~  1 + RT*Baseline + (Baseline-1|ppNames)',...
    'FitMethod','Laplace');

tbl.Weak  = all(double(tbl.BetaPreresponse.Evidencestrength) == 1 & double(tbl.BetaPreresponse.Context) == 2,:);
lme.Weak2 = fitglme(tbl.Weak,  'Peak ~  1 + RT*Baseline + (Baseline-1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong  = all(double(tbl.BetaPreresponse.Evidencestrength) == 2 & double(tbl.BetaPreresponse.Context) == 1,:);
lme.Strong1 = fitglme(tbl.Strong,  'Peak ~  1 + RT*Baseline + (Baseline-1|ppNames)',...
    'FitMethod','Laplace');
tbl.Strong  = all(double(tbl.BetaPreresponse.Evidencestrength) == 2 & double(tbl.BetaPreresponse.Context) == 2,:);
lme.Strong2 = fitglme(tbl.Strong,  'Peak ~  1 + RT*Baseline + (Baseline-1|ppNames)',...
    'FitMethod','Laplace');

diary(fullfile(parameters.logFolder, 'WaveformsStats', 'lmeERD_preResponseconvBaseline.txt'))
clc
fprintf('--------------------------------- WEAK fixed -------------------------------------------\n');
lme.Weak1 

fprintf('--------------------------------- WEAK fixed -------------------------------------------\n');
lme.Weak2

fprintf('--------------------------------- Strong fixed -------------------------------------------\n');
lme.Strong1 

fprintf('--------------------------------- Strong fixed -------------------------------------------\n');
lme.Strong2
fprintf('--------------------------------------------------------------------------------------------------------------\n');
diary off
%}


% 
% 
% 
% tbl.Strong = tblAll(double(tblAll.Evidencestrength) == 2 & double(tblAll.Context) == 1,:);
% lme.Strong1 = fitglme(tbl.Strong,  'RT ~  1 + Preperation + (1|ppNames)',...
%     'FitMethod','Laplace');
% tbl.Strong = tblAll(double(tblAll.Evidencestrength) == 2 & double(tblAll.Context) == 2,:);
% lme.Strong2 = fitglme(tbl.Strong,  'RT ~  1 + Preperation + (1|ppNames)',...
%     'FitMethod','Laplace');



%tblN2.Excursion = tbl.Excursion.Peak;

%% ------------------- get SSVEP -----------------------------------------
% just to check for bleading and ensure that the baseline difference is not
% a effect of SSVEP.
%{
keyboard
trangeTopo  = 0 + [-2 -1]*1000/parameters.stim.freqSSVEP;
AverageOrSlope   = 1;
TargetOrResponse = 1;
baselineCorrect  = 0;

parameters.eeg.baseline =  0 + [-2 0]*1000/parameters.stim.freqSSVEP; % Motor Execution threshold
parameters.figLayOut.targetLim   = [-200 0:500:1000];
parameters.figLayOut.plotCI = 0.05;
parameters.figLayOut.saveDim     = [5 11];

methodUsed       = 1;
grouping         = 0;
plotCondition    = [1 2];
plotThis         = 3; % 1) plot individual, 2) plot only averages.


[Quality, fig] = parameters.plotERP(methodUsed, 'SSVEP', [],...
    baselineCorrect, plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse, AverageOrSlope, 1, 15);
fig.Info{1}(1,2).select()
legend('off')
fig.Info{1}(1,1).select()
getAxes = gca;
shaded = shadedErrorBar( 0 + [-2 0]*1000/parameters.stim.freqSSVEP, repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle  = 'none'; shaded.edge(2).LineStyle  = 'none';
shaded.patch.FaceAlpha    = 0.1;
plotSave(gca, 'SSVEP.png', fullfile(parameters.figFolder, 'groupAverage\HPF0_CSD1\SSVEP\'), parameters.figLayOut.saveDim)


parameters.plotDifferenceTopo('SSVEP', [], [2 1], trangeTopo, TargetOrResponse, AverageOrSlope,1, [15:30]);
%}

%% ------------------- get N2 -----------------------------------------
% define topoplot information
keyboard
trangeTopo  = [240 340] ;% 250 + [-2 -1]*1000/parameters.stim.freqSSVEP; % note if you add a second row, these will be used for FA plots. Only nesc. when you are using target-locked topos
TargetOrResponse = 1;       % as written 1) target-locked 2) response-locked

% define rest
baselineCorrect  = 1;
methodUsed       = 1;       % Method used can be 1) choose spec. electrodes, 2) choose cluster of electrodes and get the 3 best based on SNR, 3) apply lucalize
plotThis         = 6;      % 1) individual plots, 
                            % 2) average Hits and seperated average Misses and FA
                            % 3) include all average target-locked Misses
                            % 4) plot target and response seperately.
                            % 5) get first derivative
grouping         = 0;       % if you want to include seperated plots for a specific condition (NOTE THIS should be your first in the plotCondition)
plotCondition    = [1 2];
% Plot target and response locked N2 (target-locked including all misses)
% as a function of target context as well as target evidence strength.
%
[~, fig]  = parameters.plotERP(methodUsed, 'N2',  ['D30 D31 B12 B11'], baselineCorrect,...
   plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse, AverageOrSlope);

plotThis         = 2;
[~, fig]  = parameters.plotERP(methodUsed, 'N2',  ['D30 D31 B12 B11'], baselineCorrect,...
   plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse, AverageOrSlope);


figure(fig.Handle{1});
legend('off')
ylim([-26 5])
getAxes = gca;
ylabel('')
%}
%{
[hypo, timeSeries] = compareEffects_context(parameters, methodUsed, 'N2', baselineCorrect, plotCondition, TargetOrResponse, 50, 10, 0.001, 0);

context = hypo(5,:);
context(:, hypo(end,:) ~= 0) = 0;
context = find(context ~= 0);
%
shaded = shadedErrorBar(parameters.eeg.targetEpoch([ timeSeries(context(1)) - 50/2 timeSeries(context(end))+ 50/2]),...
    repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle = 'none'; shaded.edge(2).LineStyle = 'none';
shaded.patch.FaceAlpha = 0.05;
legend('off')
      
% 
plotSave(gca, 'N2_shaded.png', fullfile(parameters.figFolder, 'groupAverage\HPF0_CSD1\N2\'), [4.1 6.1])
%}
% parameters.plotDifferenceTopo('N2',  ['D30 D31 B12 B11'], [2 1], trangeTopo, TargetOrResponse, AverageOrSlope);
%%
hereFigFolder = fullfile(parameters.figFolder,  'groupAverage',...
    ['HPF' num2str(parameters.eeg.HPFcutoff) '_CSD' num2str(parameters.eeg.applyCSD)],...
    'N2', 'RTBINS2/');

TargetOrResponse = 1;
trangeTopo(1,:)  = parameters.eeg.targetEpoch(timeSeries(context)) + [-1 1]*1000/parameters.stim.freqSSVEP;
trangeTopo(2,:)  = [0 300];
PeakMeanOrMax    = 2;
plotCondition    = [1 2];
negOrPos         = 1;

parameters.getWaveformParameters(methodUsed, 'N2', [1/3 2/3], 1,...
    1, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos, PeakMeanOrMax);

set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
yticks([-20:10:0])
xlim([300 1050])
xticks(400:200:1100);
xticklabels((400:200:1100)./1000)
xlabel(' ')
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
% set(gca, 'YDir','reverse')
plotSave(gca, ['peak.png'], hereFigFolder, [4.1 3.8]);

tblN2 = parameters.getWaveformParameters(methodUsed, 'N2', [], baselineCorrect,...
    0, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos, PeakMeanOrMax);


lme.N2 = fitglme(tblN2,  'Peak ~  1 + Context*Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');


%{
% check context
tbl.Weak = tblN2(double(tblN2.Evidencestrength) == 1,:);
lme.Weak = fitglme(tbl.Weak,  'Peak ~  1 + Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong= tblN2(double(tblN2.Evidencestrength) == 2,:);
lme.Strong = fitglme(tbl.Strong,  'Peak ~  1 + Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');

% check evidence strength
tbl.Fixed = tblN2(double(tblN2.Context) == 1,:);
lme.Fixed  = fitglme(tbl.Fixed ,  'Peak ~  1 + Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.Mixed = tblN2(double(tblN2.Context) == 2,:);
lme.Mixed = fitglme(tbl.Mixed,  'Peak ~  1 + Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');   
%}

keyboard

%}
%% Possible mdeiation methods
%{
CPPCovariant = getCovariant(parameters, 2, 'CPP', baselineCorrect, plotThis, plotCondition, TargetOrResponse, tblN2.PeakLatency);
%}
%
tblAll = tbl.BetaBaseline;
tblAll.Preperation = tbl.BetaBaseline.Peak;
tblAll.Execution   = tbl.BetaPreresponse.Peak;
tblAll.N2          = tblN2.Peak;
tblAll.buildUp     = tblCPP.Peak;


close all
addpath(genpath('C:\Users\eng\Documents\MATLAB\dataAnalysis\Mediation'))
savePath = fullfile(parameters.figFolder, 'Stats', 'N2-CPP-RT');
if ~isfolder(savePath), mkdir(savePath); end


conditionNames = {'Weak fixed', 'Strong Fixed', 'Weak Mixed', 'Strong Mixed'};

Names = {'N2' 'RT' 'buildUp'};
RTMediation(parameters, tblAll, 'N2', 'RT', 'buildUp',[], [], 0, Names,...
    savePath , {'Evidencestrength'}, conditionNames)



%% ------------------- get CPP anterior ----------------------------------------------
%{
    trangeTopo  =  0 + [1 2]*1000/parameters.stim.freqSSVEP;
    TargetOrResponse = 2;
    AverageOrSlope   = 1;
    baselineCorrect  = 1;
    methodUsed       = 1;
    plotThis         = 3;
    grouping         = 0;
    
    [Quality, figInfo, panelInfo] = parameters.plotERP(methodUsed, 'CPP_anterior', 'A1 A2 B1 D15 D1 C1', baselineCorrect, plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse, AverageOrSlope);
    
    figure(figInfo{1})
    figInfo{1}(1,2).select()
    getAxes = gca;
    shaded = shadedErrorBar(trangeTopo, repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
    shaded.mainLine.LineStyle = 'none';
    plotSave(gca, 'CPP_anterior_shaded.png', 'C:\Users\eng\Documents\2. Postdoc_Kelly\Project 5 - Neural correlates of static and dynamic decision-bound adjustments\Project 5.1 - Signatures of bound adjustment\Data\Figures\groupAverage\HPF0_CSD1\', parameters.figLayOut.saveDim)

    parameters.plotDifferenceTopo('CPP_anterior', 'A1 A2 B1 D15 D1 C1', [2 1], trangeTopo, TargetOrResponse, AverageOrSlope); %'A4 A5 A19 A32 A18 A20 A31 A17 A21 A30'

        
    % get waveform parameters.
    trangeTopo(1,:) =  0 + [1 2]*1000/parameters.stim.freqSSVEP;
    trangeTopo(2,:) =  [-300 -1*1000/parameters.stim.freqSSVEP];

    tblCPP = parameters.getWaveformParameters(methodUsed, 'CPP_anterior', 'A1 A2 B1 D15 D1 C1' , baselineCorrect,...
        plotThis, plotCondition, grouping, trangeTopo);
%}
