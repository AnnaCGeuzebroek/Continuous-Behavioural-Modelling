classdef dataAnalysis < handle
    
    properties
        % set up parameters
        dataFolders
        inputFolder
        outputFolder
        logFolder
        figFolder
        
        dataFiles
        
        ppNames
        conditions
        condNames
        numBlocks
        numTrials
        order
        
        analysisEEG
        analysisEOG
        analysisEMG
        analysisEyelink
        
        stim
        triggers
        eeg
        system
        event
        experiment
        behaviour
        modelBehaviour
        DetectOrDisc
        
        eyelink
        
        % plotting information
        figLayOut
    end
    
    methods(Access = public)
        
        function getInformation(obj, preprocess, inputFolder, conditions, condNames, blockRandomized)
            %% getInformation(obj, inputFolder, conditions, condNames)
            % set all information preset in the object to later be used
            % throughout the functions.
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -----------------   Conditions  ----------------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % now still limited to the conditions that are block
            % randomized. Important to reimplement the other randomized
            % parameters. i.e. here the inter-stimuli-times.
            
            obj.conditions = conditions;
            obj.condNames  = condNames;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -------------- Folders and Participants  -------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get input folder and created the other needed folders e.g.
            % Processed data, logs and figure folders.
            obj.inputFolder	 = fullfile(inputFolder, 'Raw');
            obj.outputFolder = fullfile(inputFolder, 'Processed');  if ~exist(obj.outputFolder, 'dir'), mkdir(obj.outputFolder); end
            obj.logFolder    = fullfile(inputFolder, 'Logs');       if ~exist(obj.logFolder, 'dir'),    mkdir(obj.logFolder);    end
            obj.figFolder    = fullfile(inputFolder, 'Figures');    if ~exist(obj.figFolder, 'dir'),    mkdir(obj.figFolder);    end
            
            obj.dataFolders{1} = 'Trial files';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -------------- Set-up for EEG    ---------------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % check in the file if EEG is used, if so set-up the parameters
            % in EEGlab
            if obj.analysisEEG
                obj.dataFolders{2} = 'EEG data';
                switch lower(obj.system)
                    case 'brainvision'
                        % initiated eeglab
                        if ~exist('pop_loadbv', 'file')
                            eeglab;
                            close all;
                        end
                        obj.eeg.extention = '.vhdr';
                        obj.eeg.function  = @(folder, fileName)pop_loadbv(folder, fileName);
                    case 'biosemi'
                        if ~exist('pop_biosig', 'file')
                            eeglab;
                            close all;
                        end
                        obj.eeg.extention = '.bdf';
                        obj.eeg.function  = @(folder, fileName)pop_biosig(folder);
                    otherwise
                        obj.eeg.extention = obj.system;
                        obj.eeg.function  = @(folder, fileName)load(folder);
                end
                
                % create output folder for EEG data specific.
                if ~exist(fullfile(obj.outputFolder, 'EEG data'), 'dir')
                    mkdir(fullfile(obj.outputFolder, 'EEG data'));
                    mkdir(fullfile(obj.outputFolder, 'EEG data', 'groupAverage'))
                end
            else
                % create output folder for EEG data specific.
                if ~exist(fullfile(obj.outputFolder, 'Behavioural data'), 'dir')
                    mkdir(fullfile(obj.outputFolder, 'Behavioural data'));
                end
            end
            
            
            if obj.analysisEyelink
                obj.dataFolders{3} = 'Eyelink data';
                % initiated Gigdet
                if ~exist('gigadet', 'file')
                    addpath(genpath('C:\Users\eng\Documents\MATLAB\Toolbox_Visual\GigaDet'))
                end
                if ~exist(fullfile(obj.outputFolder, 'Eyelink data'), 'dir'), mkdir(fullfile(obj.outputFolder, obj.dataFolders{3})); end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -------------- Experimental information  -------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % always start from initial trial information .mat files saved
            % in './Trial files/' folder
            
            if preprocess == 1
                % get all participants names.
                obj.ppNames = dir(obj.inputFolder);
                obj.ppNames = {obj.ppNames.name};
                obj.ppNames(strncmp(obj.ppNames, '.',1)) = [];
                
                deletePP    = zeros(length(obj.ppNames),1);
                
                % run through all the participants
                for indPP = 1:length(obj.ppNames)
                    % preset current-folder
                    currFolder = fullfile(obj.inputFolder, obj.ppNames{indPP}, obj.dataFolders{1});
                    
                    obj.dataFiles{indPP} = dir(fullfile(currFolder, '*.mat'));
                    obj.dataFiles{indPP} = {obj.dataFiles{indPP}.name};
                    obj.dataFiles{indPP}(strncmp(obj.dataFiles{indPP}, '.',1)) = [];
                    
                    fileID = fopen(fullfile(obj.logFolder, [obj.ppNames{indPP} '.log']),'a');
                    if isempty(obj.dataFiles{indPP})
                        % when there wasn't a datafile
                        warning(['No data file detected for participant: ' obj.ppNames{indPP}])
                        fprintf(fileID, '-------------- %s --------------\n', obj.ppNames{indPP});
                        fprintf(fileID, 'No data file detected for participant\n');
                        fprintf(fileID, 'This participant is excluded and will not be further analysed');
                        obj.numBlocks(indPP,:) = 0;
                        obj.experiment{indPP}  = [];
                        deletePP(indPP) = 1;
                    else
                        obj.numBlocks(indPP,:) = length(obj.dataFiles{indPP});
                        fprintf(fileID, '-------------- %s --------------\n', obj.ppNames{indPP});
                        fprintf(fileID, 'Data file detected for participant\n');
                        fprintf(fileID, 'This participant is included, under condition that all blocks have the following:\n');
                        fprintf(fileID, 'EEG data and if eye movements are measured Eyelink files\n');
                        
                        % get all blocks and check if the data files are
                        % complete. e.g check for EEG data or Eyelink data if
                        % those are requested for the analysis.
                        deleteBlock = zeros(1,length(obj.dataFiles{indPP}));
                        for indBlock = 1:length(obj.dataFiles{indPP})
                            % load trial file to get all information of the individual
                            % blocks
                            
                            %                         if strfind(obj.ppNames{indPP}, 'PP11') & indBlock == 2
                            %                             keyboard
                            %                         end
                            blockInfo = load(fullfile(currFolder, obj.dataFiles{indPP}{indBlock}));
                            saveInfo  = fieldnames(blockInfo);
                            for indSaveInfo = 1:length(saveInfo)
                                stringToEval = sprintf('obj.experiment{%s}{%s}.%s = blockInfo.%s;',  num2str(indPP), num2str(indBlock), saveInfo{indSaveInfo}, saveInfo{indSaveInfo});
                                eval(stringToEval)
                            end
                            
                            % check if EEG should be analysed
                            if obj.analysisEEG
                                if isempty(dir(fullfile(obj.inputFolder, obj.ppNames{indPP}, obj.dataFolders{2}, [obj.dataFiles{indPP}{indBlock}(1:end-4)  obj.eeg.extention]))) %obj.dataFiles{indPP}{indBlock}(1:end-4)
                                    deleteBlock(indBlock) = 1;
                                    fprintf(fileID,'* Block %s named %s misses EEG data and is excluded completely\n', num2str(indBlock), obj.dataFiles{indPP}{indBlock});
                                end
                            end
                            
                            % check if Eyelink should be analysed
                            if obj.analysisEyelink
                                if  isempty(dir(fullfile(obj.inputFolder, obj.ppNames{indPP}, 'Eyelink data', [obj.dataFiles{indPP}{indBlock}(1:end-4) '.edf'])))
                                    deleteBlock(indBlock) = 1;
                                    fprintf(fileID,'* Block %s named %s misses only Eyelink data and is excluded completely\n', num2str(indBlock), obj.dataFiles{indPP}{indBlock});
                                end
                            end
                            
                            if ~deleteBlock(indBlock)
                                fprintf(fileID,'* Block %s named %s is complete\n', num2str(indBlock), obj.dataFiles{indPP}{indBlock});
                            end
                        end
                        
                        obj.dataFiles{indPP}(logical(deleteBlock))  = [];
                        obj.experiment{indPP}(logical(deleteBlock)) = [];
                        obj.numBlocks(indPP) = obj.numBlocks(indPP) - sum(logical(deleteBlock));
                        
                        % move is participant if none of the blocks was
                        % complete
                        
                        if sum(deleteBlock) == length(deleteBlock)
                            % when there wasn't a datafile
                            warning(['No data file detected for participant: ' obj.ppNames{indPP}])
                            fprintf(fileID, '-------------- %s --------------\n', obj.ppNames{indPP});
                            fprintf(fileID, 'No data file detected for participant\n');
                            fprintf(fileID, 'This participant is excluded and will not be further analysed');
                            deletePP(indPP) = 1;
                        end
                        fclose(fileID);
                    end
                end
                
                deletePP = logical(deletePP);
                obj.ppNames(deletePP)    = [];
                obj.dataFiles(deletePP)  = [];
                obj.numBlocks(deletePP)  = [];
                obj.experiment(deletePP) = [];
                
                % find block order using the Eyelink timestamp.
                if blockRandomized && obj.analysisEyelink && ~exist(fullfile(obj.outputFolder, 'Trial files', 'BlockOrder.mat'), 'file')
                    for indPP = 1:length(obj.ppNames)
                        currentFile = fullfile(obj.inputFolder, obj.ppNames{indPP}, 'Eyelink data');
                        eyeFiles    = dir(fullfile(currentFile, '*.edf')); eyeFiles = {eyeFiles(:).name};
                        
                        time = NaT(length(eyeFiles),1);
                        for indBlock = 1:length(eyeFiles)
                            
                            checkRealBlock = find(strncmp(obj.dataFiles{indPP}, eyeFiles{indBlock}(1:end-4), length(eyeFiles{indBlock}(1:end-4))));
                            if length(checkRealBlock) > 1, checkRealBlock(1:end-1) = []; end
                            clear checkRealBlock
                            
                            % if Eyelink, then we can use it to see the order in
                            % which it is collected
                            [ ~, collectInfo]= edfImport(fullfile(currentFile, eyeFiles{indBlock}), [1 1 1]);
                            dateInd = strfind(collectInfo, 'DATE');
                            endInd = strfind(collectInfo, '**');
                            
                            time(indBlock) = datetime(collectInfo(dateInd+10:endInd(2)-2),'InputFormat','MM dd HH:mm:ss yy');
                            
                        end
                        timeTable = table(time,(1:length(eyeFiles))', 'VariableNames', {'time', 'blockNum'});
                        timeTable = sortrows(timeTable,'time');
                        obj.order{indPP} = reshape(repmat(timeTable.blockNum', obj.numTrials,1), [], 1);
                    end
                    
                    BlockOrder = obj.order;
                    save(fullfile(obj.outputFolder, 'Trial files', 'BlockOrder.mat'), 'BlockOrder');
                    
                elseif blockRandomized && exist(fullfile(obj.outputFolder, 'Trial files', 'BlockOrder.mat'), 'file')
                    load(fullfile(obj.outputFolder, 'Trial files', 'BlockOrder.mat'), 'BlockOrder');
                    obj.order = BlockOrder;
                else
                    for indPP = 1:length(obj.ppNames)
                        for indBlock = 1:length(obj.dataFiles{indPP})
                            A = regexp(obj.dataFiles{indPP}{indBlock}(length(obj.ppNames{indPP})+1:end),'\d*','Match');
                            order(str2double(A{:})) = indBlock;
                        end
                        
                        obj.order{indPP} = reshape(repmat(order,obj.numTrials,1),1,[])';
                    end
                end
            else
                % get all participants names.
                obj.ppNames = dir( fullfile(obj.outputFolder, 'EEG data'));
                obj.ppNames = {obj.ppNames.name};
                obj.ppNames(strncmp(obj.ppNames, '.',1)) = [];
                obj.ppNames(strcmp(obj.ppNames, 'groupAverage')) = [];
            end
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% -----------------   Eyelink processing   ---------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function extractEyelink(obj)
            %% extractEyelink
            % extractEyelink can be used to scroll through participants
            % eyelink data files (.edf). It checks if the .edf is
            % presented and if it has already been visual processed, e.g.
            % if a .mat file exist. This requires the researcher to
            % manually go through the file, and if nesceccary adapt the DET
            % file used in GIGADET.
            % This function is dependent on GIGADET as written by Jeroen
            % Goossens and adapted by Anna Geuzebroek. download the Eyelink
            % matlab functions written by Cornellisen from SR-Research.
            % Be sure to check if matlab is able to read .mex files
            % (check by mex -setup).
            
            % now load the file in gigadet and visually inspect data. Also
            % save .DET file in the .log file.
            
            for indPP = 1:length(obj.ppNames)
                currentFile = fullfile(obj.inputFolder, obj.ppNames{indPP}, obj.dataFolders{3});
                eyeFiles    = dir(fullfile(currentFile, '*.edf')); eyeFiles = {eyeFiles(:).name};
                fileID      = fopen(fullfile(obj.logFolder, [obj.ppNames{indPP} '.log']),'a');
                
                for indBlock = 1:length(eyeFiles)
                    
                    checkRealBlock = find(strncmp(obj.dataFiles{indPP}, eyeFiles{indBlock}(1:end-4), length(eyeFiles{indBlock}(1:end-4))));
                    if length(checkRealBlock) > 1, checkRealBlock(1:end-1) = []; end
                    
                    if ~exist(fullfile(obj.outputFolder, 'Eyelink data', obj.ppNames{indPP}, [obj.ppNames{indPP} '_' num2str(checkRealBlock) '_preprocessedEyelink.mat' ]), 'file')
                        mkdir(fullfile(obj.outputFolder, 'Eyelink data', obj.ppNames{indPP}));
                        fprintf(fileID,'----------------- EYELINK -------------------\n');
                        fprintf(fileID,'Eyelink data is analysis to extract blinks and saccades\n');
                        
                        Trials = edfImport(fullfile(currentFile, eyeFiles{indBlock}), [1 1 1]);
                        
                        %% timeseries and event times of the Eyelink and for EEG seperately
                        EL.timeSeries = double(Trials.Samples.time  - Trials.Samples.time(1));
                        eventTimes = Trials.Events.sttime - Trials.Samples.time(1);
                        
                        % Determine the samplingRate of the Eyelink:
                        % 1) EEG block as eeg tiggers are reliable.
                        try
                            BlockLength = obj.event{indPP}{checkRealBlock}.Times(end) - obj.event{indPP}{checkRealBlock}.Times(obj.event{indPP}{checkRealBlock}.Number == 1);
                            BlockLength = BlockLength/obj.eeg.SampleRate;
                        catch
                            keyboard
                        end
                        % first detect blockStart using TASK_START which is
                        % send out at the same time as (1) to the EEG system
                        % and is this the synchtime.
                        for indEvent = 1:length(Trials.Events.message)
                            if ~isempty(Trials.Events.message{indEvent}) && strcmp(Trials.Events.message{indEvent}, 'TASK_START')
                                blockStart = find(EL.timeSeries == eventTimes(indEvent));
                            end
                        end
                        
                        % rerun through the message to get trialStarts (later
                        % used for artefact rejection) as well as the blockend.
                        % And use block start to cut off the trials.
                        
                        trialStart = [];
                        for indEvent = 1:length(Trials.Events.message)
                            if ~isempty(Trials.Events.message{indEvent}) && strcmp(Trials.Events.message{indEvent}(end), num2str(obj.triggers.stimulusON))
                                trialStart(end+1) = find(EL.timeSeries == eventTimes(indEvent)) - blockStart;
                            elseif ~isempty(Trials.Events.message{indEvent})
                                blockEnd   = find(EL.timeSeries == eventTimes(indEvent)) - blockStart;
                            end
                        end
                        
                        % downsampling
                        EL.SampleRate  = BlockLength/blockEnd;
                        StepSize       = EL.SampleRate;
                        
                        EL.timeSeries(blockStart:end) = [0:StepSize:(length(EL.timeSeries(blockStart:end))*StepSize-StepSize)];
                        EL.timeSeries(1:blockStart-1) = NaN;
                        
                        % use blocklength to get actual blockEnd
                        [~, blockEnd] = min(abs(EL.timeSeries - (BlockLength + EL.SampleRate*2)));
                        EL.timeSeries = EL.timeSeries(blockStart:blockEnd);
                        
                        % cutting all important sequences.
                        cathesian  = [Trials.Samples.hx(2, blockStart:blockEnd); Trials.Samples.hy(2, blockStart:blockEnd)];
                        pupilSize  = Trials.Samples.pa(2, blockStart:blockEnd);
                        
                        %% timeseries and event times of the Eyelink and for EEG seperately
                        % get amplitude, this the rx/ry
                        amplitude = cart2sph(15000, cathesian(1,:), -cathesian(2,:));
                        
                        % get blinks
                        blinks      = detrend(pupilSize);
                        mBlinks     = mean(blinks); thresBlinks = mBlinks - (2*std(blinks));
                        indBlinks   = logical([0 blinks < thresBlinks 0]); blinks(indBlinks(2:end-1)) = nan;
                        blinkOnset  = find(diff(indBlinks) > 0);
                        blinkOffset = find(diff(indBlinks) < 0);
                        
                        amplitude(indBlinks(2:end-1)) = nan;
                        
                        for indBlink = 1:length(blinkOnset)
                            % figure, plot(abs(diff(blinks(blinkOnset(indBlink)-150:blinkOffset(indBlink)+150))).^2)
                            % hold on, plot(blinks(blinkOnset(indBlink)-150:blinkOffset(indBlink)+150))
                            try
                                turningPoint = find((abs(diff(blinks(blinkOnset(indBlink)-150:blinkOffset(indBlink)+150))) > 4));
                            catch
                                try
                                    turningPoint = find((abs(diff(blinks(blinkOnset(indBlink)-150:blinkOffset(indBlink)))) > 4));
                                catch
                                    try
                                        turningPoint = find((abs(diff(blinks(blinkOnset(indBlink):blinkOffset(indBlink)+150))) > 4));
                                    catch
                                        turningPoint = find((abs(diff(blinks(blinkOnset(indBlink):end))) > 4));
                                    end
                                end
                            end
                            
                            try
                                blinks((blinkOnset(indBlink)-150 + turningPoint(1)):(blinkOnset(indBlink)-150 + turningPoint(end))) = nan;
                                amplitude((blinkOnset(indBlink)-150 + turningPoint(1)):(blinkOnset(indBlink)-150 + turningPoint(end))) = nan;
                            end
                        end
                        EL.indBlinks = isnan(blinks);
                        
                        % calculated velocity and use it to get onset/offset
                        velocity     = diff(amplitude)./StepSize;
                        medianVel    = nanmedian(velocity);
                        stdVel       = nanstd(velocity);
                        thresholdVel = medianVel + stdVel.*6;
                        indSaccade   = [0 abs(velocity) > thresholdVel 0];
                        %{
                    figure, plot(EL.timeSeries(1:end-1), abs(EL.velocity))
                    hold on,
                    tmpLine = line([EL.timeSeries(1) EL.timeSeries(end-1)], [medianVel medianVel]);       tmpLine.Color = 'r';
                    tmpLine = line([EL.timeSeries(1) EL.timeSeries(end-1)], [thresholdVel thresholdVel]); tmpLine.Color = 'r';
                    plot(EL.timeSeries(1:end-1), abs(EL.velocity).*(abs(EL.velocity) > thresholdVel), '*' )
                        %}
                        
                        % index sacOnset and Offset
                        saccOnset  = find(diff(indSaccade) > 0);
                        saccOffset = find(diff(indSaccade) < 0);
                        timingEL = [EL.timeSeries(saccOnset); EL.timeSeries(saccOffset)];
                        
                        % Saccade length in ms.
                        saccLength = abs(timingEL(1,:) - timingEL(2,:));
                        
                        % delete saccades that are shorter than 8 ms. TODO not
                        % hardprogram.
                        saccOnset(saccLength  < 8) = [];
                        saccOffset(saccLength < 8) = [];
                        timingEL = [EL.timeSeries(saccOnset); EL.timeSeries(saccOffset)];
                        
                        % delete saccades that follow shortly after the
                        % saccades. The inter saccadic period needs to be
                        % longer than 50 ms.
                        saccOnset(logical( [0 (timingEL(1, 2:end) - timingEL(2,1:end-1)) < 50])) = [];
                        saccOffset(logical([0 (timingEL(1, 2:end) - timingEL(2,1:end-1)) < 50])) = [];
                        
                        fprintf(fileID,'* block %s: %s saccades detected\n', num2str(indBlock), num2str(length(saccOnset)));
                        
                        EL.indSacc = zeros(1, length(EL.timeSeries));
                        for indSac = 1:length(saccOnset)
                            
                            % figure, plot(amplitude(saccOnset(indSac)-150:saccOffset(indSac)+150))
                            try
                                EL.saccAmp(indSac) = abs(diff([nanmean(amplitude(saccOnset(indSac)-10:saccOnset(indSac)))  nanmean(amplitude(saccOffset(indSac):saccOffset(indSac)+10))]));
                                EL.saccVel(indSac) = max(abs(diff(amplitude(saccOnset(indSac)-10:saccOffset(indSac)+10))./StepSize));
                                EL.indSacc(saccOnset(indSac):saccOffset(indSac)) = 1;
                            catch
                                try
                                    EL.saccAmp(indSac) = abs(diff([nanmean(amplitude(saccOnset(indSac):saccOnset(indSac)))  nanmean(amplitude(saccOffset(indSac):saccOffset(indSac)+10))]));
                                    EL.saccVel(indSac) = max(abs(diff(amplitude(saccOnset(indSac):saccOffset(indSac)+10))./StepSize));
                                    EL.indSacc(saccOnset(indSac):saccOffset(indSac)) = 1;
                                catch
                                    EL.saccAmp(indSac) = abs(diff([nanmean(amplitude(saccOnset(indSac)-10:saccOnset(indSac)))  nanmean(amplitude(saccOffset(indSac):saccOffset(indSac)))]));
                                    EL.saccVel(indSac) = max(abs(diff(amplitude(saccOnset(indSac)-10:saccOffset(indSac)))./StepSize));
                                    EL.indSacc(saccOnset(indSac):saccOffset(indSac)) = 1;
                                end
                            end
                        end
                        
                        EL.indSacc = logical(EL.indSacc);
                        obj.eyelink{indPP}{checkRealBlock} = EL;
                        %{
                        figure, hold on,
                        plot(EL.timeSeries, amplitude)
                        plot(EL.timeSeries(saccOnset), amplitude(saccOnset), '*')
                        
                        saveas(gcf,fullfile(obj.figFolder, ['Eyelink' obj.ppNames{indPP} '_' num2str(indBlock) '.png']))
                        %}
                        save(fullfile(obj.outputFolder, 'Eyelink data',  obj.ppNames{indPP}, [obj.ppNames{indPP} '_' num2str(checkRealBlock) '_preprocessedEyelink' ]), 'EL');
                        close all
                        
                    else
                        load(fullfile(obj.outputFolder, 'Eyelink data',  obj.ppNames{indPP}, [obj.ppNames{indPP} '_' num2str(checkRealBlock) '_preprocessedEyelink' ]), 'EL');
                        obj.eyelink{indPP}{checkRealBlock} = EL;
                        % hold on, plot(EL.saccAmp, EL.saccVel, '*');
                    end
                end
                fclose(fileID);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% -----------------   EEG peprocessing   ---------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dependables for preprocessing, set-up parameters for the preprocessing including filters ect.
        
        function applyPreprocessing(obj)
            %% preprocessing(obj)
            % This function goes throught the inputfolder and get all participants
            % names and intially search .eeg data in this folder. If not it will look
            % for the EEG folder with a subfolder data and grabs the .eeg dat
            % there independent on the task.
            % To use this part of the code you need eeglab with the
            % appropiated extention: Biosig for the Biosemi system and
            % bva-io for the BrainVision system.
            if obj.analysisEEG
                
                % preset the filters.
                obj.setFilters;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% ----------------    Pre-processe eeg data  -------------
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Now loop through the subjects and extract the single trial
                % ERPs and the other important information about every trial:
                
                for indPP = 1:length(obj.ppNames)
                    fprintf('Get pre-processing participant %i out of %i\n', indPP, length(obj.ppNames))
                    close all
                    
                    currInput     = fullfile(obj.inputFolder, obj.ppNames{indPP}, obj.dataFolders{2});
                    currOutput    = fullfile(obj.outputFolder, obj.dataFolders{2}, obj.ppNames{indPP});
                    currFigOutput = fullfile(obj.figFolder, 'preProcessing data', ['HPF_' num2str(obj.eeg.HPFcutoff)], obj.ppNames{indPP});
                    currEOGOutput = fullfile(obj.outputFolder, 'EOG data', obj.ppNames{indPP});
                    currEMGOutput = fullfile(obj.outputFolder, 'EMG data', obj.ppNames{indPP});
                    
                    if ~exist(currOutput, 'dir'); mkdir(currOutput); end
                    if obj.analysisEOG && ~exist(currEOGOutput, 'dir');     mkdir(currEOGOutput); end
                    if obj.analysisEMG && ~exist(currEMGOutput, 'dir');     mkdir(currEMGOutput); end
                    if obj.eeg.badChannel && ~exist(currFigOutput, 'dir');  mkdir(currFigOutput); end
                    
                    
                    % already set for the possibility to work with Biosemi.
                    blockFiles = dir(fullfile(currInput, ['PP*' obj.eeg.extention])); blockFiles = {blockFiles.name};
                    fileID = fopen(fullfile(obj.logFolder, [obj.ppNames{indPP} '.log']),'a');
                    
                    blockCounter = 1;
                    while blockCounter <= length(blockFiles)
                        
                        % save as two options to get first the output and
                        % then bad channel output.
                        currFilter   = fullfile(currOutput, [obj.ppNames{indPP} '_' num2str(blockCounter) '_filteredEEG_HPF' num2str(obj.eeg.HPFcutoff) '.mat' ]);
                        currPrepross = fullfile(currOutput, [obj.ppNames{indPP} '_' num2str(blockCounter) '_preprocessedEEG_HPF' num2str(obj.eeg.HPFcutoff) '.mat' ]);
                        
                        % check for LPF data (e.g. pre-processed data)
                        currFilterLPF = fullfile(currOutput, [obj.ppNames{indPP} '_' num2str(blockCounter) '_filteredEEG_HPF0.mat' ]);
                        if ~exist(currPrepross, 'file')
                            if ~exist(currFilter, 'file')
                                
                                if blockCounter == 1
                                    fprintf(fileID,'----------------- PREPROCESSING -------------------\n');
                                    fprintf(fileID,'1) Cutting access EEG after last event\n');
                                    
                                    if obj.eeg.applyLPF
                                        fprintf(fileID, '2) Apply low-pass filter of %s with a L Hamming window of %s\n', num2str(obj.eeg.cuttoffLPF), num2str(obj.eeg.LHamming));
                                    else
                                        fprintf(fileID, '2) No low-pass filter applied\n');
                                    end
                                    if obj.eeg.applyDetrend
                                        fprintf(fileID,'3) Apply detrending taking possible DCC into account\n');
                                    else
                                        fprintf(fileID,'3) No detrending applied\n');
                                    end
                                    if obj.eeg.applyHPF
                                        fprintf(fileID, '4) Apply high-pass filter of %s\n', num2str(obj.eeg.HPFcutoff));
                                    else
                                        fprintf(fileID, '4) No high-pass filter applied\n');
                                    end
                                end
                                
                                if ~exist(currFilterLPF, 'file') || obj.eeg.applyLPF == 0
                                    switch lower(obj.system)
                                        case 'brainvision'
                                            keyboard
                                            EEG = obj.eeg.function(currInput, blockFiles{blockCounter}); % read in EEG - this is an EEGLAB function that outputs a structure 'EEG' with lots of fileds, the most important of which is EEG.data - the actual EEG data!
                                            findEvents.Number = [];
                                            findEvents.Times  = [];
                                            badchannels.DCC   = [];   % DC onset
                                            
                                            for indEvent = 1:length(EEG.event)
                                                codestr = EEG.event(indEvent).code;
                                                if contains(codestr, 'Stimulus')
                                                    findEvents.Number(end+1) = str2double(regexp(EEG.event(indEvent).type, '\d*', 'Match')); % trigger codes
                                                    findEvents.Times(end+1)  = round(EEG.event(indEvent).latency);	% sample points in continuous data at which the triggers were sent
                                                elseif contains(codestr, 'DC')
                                                    badchannels.DCC(end+1) = round(EEG.event(indEvent).latency);
                                                end
                                            end
                                            
                                            % Trim data before first stimuli presentation and
                                            % after last trigger.
                                            if findEvents.Times(end)+obj.eeg.SampleRate*2 < size(EEG.data,2)
                                                EEG.data(:,findEvents.Times(end)+obj.eeg.SampleRate*2:end) = []; % get rid of any extra EEG beyond 2s after last trigger (subject might move or jump or dance)
                                            end
                                            
                                            EEG.data(:,1:findEvents.Times(findEvents.Number == 1)) = [];
                                            
                                            % remove first nummbers until
                                            % task_trial (also to align eyelink data)
                                            findEvents.Times = (findEvents.Times - findEvents.Times(findEvents.Number == 1));
                                            
                                            % For the 2-s of data following any DCC, set to nan
                                            for d=1:length(badchannels.DCC)
                                                EEG.data(:, badchannels.DCC(d)+(1:obj.eeg.postDCC)) = nan;
                                            end
                                            
                                        case 'biosemi'
                                            
                                            EEG = obj.eeg.function(fullfile(currInput, blockFiles{blockCounter})); % read in EEG - this is an EEGLAB function that outputs a structure 'EEG' with lots of fileds, the most important of which is EEG.data - the actual EEG data!
                                            if isempty(EEG.data)
                                                blockCounter = blockCounter + 1;
                                                continue
                                            end
                                            findEvents.Number = [];
                                            findEvents.Times  = [];
                                            badchannels.DCC   = [];   % DC onset
                                            
                                            for indEvent = 1:length(EEG.event)
                                                codestr = EEG.event(indEvent).type;
                                                
                                                findEvents.Number(end+1) = codestr; % trigger codes
                                                findEvents.Times(end+1)  = round(EEG.event(indEvent).latency);	% sample points in continuous data at which the triggers were sent
                                            end
                                            
                                            % Trim data before first stimuli presentation and
                                            % after last trigger.
                                            try
                                                if findEvents.Times(end)+obj.eeg.SampleRate*2 < size(EEG.data,2)
                                                    EEG.data(:,findEvents.Times(end)+obj.eeg.SampleRate*2:end) = []; % get rid of any extra EEG beyond 2s after last trigger (subject might move or jump or dance)
                                                end
                                            catch
                                                keyboard
                                            end
                                            if findEvents.Number(1) == obj.triggers.start
                                                EEG.data(:, 1:findEvents.Times(1)) = [];
                                                findEvents.Times = findEvents.Times - findEvents.Times(1);
                                            elseif findEvents.Number(1) ~= obj.triggers.start
                                                indStart = find(findEvents.Number == obj.triggers.start,1);
                                                if ~isempty(indStart)
                                                    EEG.data(:, 1:findEvents.Times(indStart)) = [];
                                                    findEvents.Number = findEvents.Number(indStart:end);
                                                    findEvents.Times  = findEvents.Times(indStart:end) - findEvents.Times(indStart);
                                                else
                                                    findEvents.Times  = findEvents.Times - findEvents.Times(1);
                                                end
                                            end
                                            
                                    end
                                end
                                
                                if obj.analysisEOG && ~exist(fullfile(currEOGOutput, [obj.ppNames{indPP} '_' num2str(blockCounter) '_EOG' ]), 'file')
                                    % apply a high-pass filter on EOG data,
                                    % always.. e.g. 0.5 to remove
                                    % trends. Here the high-pass
                                    % distrotion is actually
                                    % beneficial.
                                    filters = obj.eeg;
                                    filters.HPFcutoff    = 0.5;
                                    filters.applyDetrend = 0;
                                    
                                    % get the EOG data
                                    EOG = EEG; EOG.data = [];
                                    
                                    EOG.data(1,:) = diff(EEG.data(obj.eeg.NumberOfChannels+obj.eeg.channelVEOG,:)); % get vertical data
                                    EOG.data(2,:) = diff(EEG.data(obj.eeg.NumberOfChannels+obj.eeg.channelHEOG,:)); % get horizontal data
                                    
                                    
                                    EOG = obj.applyFilters(EOG, badchannels, findEvents, 0, filters); % get figure as output as well.
                                    [~, EOG.rho] = cart2pol(EOG.data(2,:), EOG.data(1,:));
                                    
                                    %% TODO add addition of using frontal electrodes if EOG isnt good.
                                    % if you want to plot it.
                                    % try this
                                    %{
                                    f = figure('visible','off');
                                    
                                    % get vertical data
                                    subplot(1,3,1); hold on;
                                    plot(EOG.data(1,:))
                                    
                                    % get vertical data
                                    subplot(1,3,2); hold on;
                                    plot(EOG.data(2,:))
                                    
                                    subplot(1,3,3); plot(detrend(EOG.rho))
                                    linkaxes
                                    
                                    saveas(gcf,fullfile(currFigOutput,[obj.ppNames{indPP} '_' num2str(blockCounter) '_EOG.png']))
                                    %}
                                    save(fullfile(currEOGOutput, [obj.ppNames{indPP} '_' num2str(blockCounter) '_EOG' ]), 'EOG');
                                end
                                
                                if obj.analysisEMG && ~exist(fullfile(currEMGOutput, [obj.ppNames{indPP} '_' num2str(blockCounter) '_EMG' ]), 'file')
                                    EMG = EEG; EMG.data =[];
                                    % get the EMG data
                                    EMG.data(1,:) = double(diff(EEG.data(obj.eeg.NumberOfChannels+obj.eeg.channelLEMG,:))); % get left data
                                    EMG.data(2,:) = double(diff(EEG.data(obj.eeg.NumberOfChannels+obj.eeg.channelREMG,:))); % get right data
                                    
                                    Fs    = obj.eeg.SampleRate;
                                    Fnyq  = Fs/2;
                                    fco   = 20;
                                    [b,a] = butter(2,fco*1.25/Fnyq);
                                    
                                    for indChannel = 1:size(EMG.data)
                                        EMG.data(indChannel,:) = double(abs(EMG.data(indChannel,:) -mean(EMG.data(indChannel,:))));
                                        EMG.data(indChannel,:) = filtfilt(b,a,EMG.data(indChannel,:));
                                    end
                                    
                                    % if you wanna plot
                                    %{
                                    % try this
                                    f = figure('visible','off');
                                    
                                    % get left movement data
                                    subplot(1,2,1); hold on;
                                    plot(1:length(EMG.data(1,:)),EMG.data(1,:),'k');
                                    xlabel('Time (s)'); ylabel('Left EMG (V)');
                                    legend('Raw (offset)','Rectified','Linear envelope');
                                    plot(findEvents.Times(findEvents.Number == 101), repmat(500, 1, sum(findEvents.Number == 101)),'r*')
                                    
                                    % get right movement data
                                    subplot(1,2,2); hold on;
                                    plot(1:length(EMG.data(2,:)),EMG.data(2,:),'k');
                                    xlabel('Time (s)'); ylabel('Right EMG (V)');
                                    legend('Raw (offset)','Rectified','Linear envelope');
                                    plot(findEvents.Times(findEvents.Number == 103), repmat(500, 1, sum(findEvents.Number == 103)),'r*')
                                    
                                    linkaxes
                                    saveas(gcf,fullfile(currFigOutput, [obj.ppNames{indPP} '_' num2str(blockCounter) '_EMG.png']))
                                    %}
                                    save(fullfile(currEMGOutput, [obj.ppNames{indPP} '_' num2str(blockCounter) '_EMG' ]), 'EMG');
                                end
                                
                                if ~exist(currFilterLPF, 'file') || obj.eeg.applyLPF == 0
                                    EEG.data(obj.eeg.NumberOfChannels+1:end,:)=[];
                                    EEG.nbchan   = obj.eeg.NumberOfChannels;
                                    EEG.chanlocs = obj.eeg.chanlocs;
                                else
                                    load(currFilterLPF);
                                end
                                
                                %% %%%%%%%%%%%%%% Preprocessing     %%%%%%%%%%%%%%%%%
                                EEG = obj.applyFilters(EEG, badchannels, findEvents, 0); % get figure as output as well.
                                
                                %% %%%%%%%%%%%%%% save preprocessed data %%%%%%%%%%%%
                                save(currFilter, 'EEG', 'findEvents', 'badchannels');
                                obj.event{indPP}{blockCounter} = findEvents;
                                fprintf('.')
                            else
                                load(currFilter, 'EEG', 'findEvents', 'badchannels');
                                % obj.event{indPP}{blockCounter} = findEvents;
                            end
                            
                            
                            %% %%%%%%%%%%%%%% find Bad channels %%%%%%%%%%%%
                            if obj.eeg.badChannel
                                fprintf(fileID,'5) Apply bad channel selection with the adapted PREP algorithm\n');
                                fprintf(fileID,'Bad are identified as well as saturated channels and interpolated\n');
                                
                                fprintf(fileID,'* Block %s: Bad or saturated channels detection\n', num2str(blockCounter));
                                
                                % try this
                                output = findNoisyChannels(EEG);
                                saveas(gcf, fullfile(currFigOutput, [obj.ppNames{indPP} '_topoPlot' num2str(blockCounter) '.png']))
                                
                                output.noisyChannels.all(output.noisyChannels.all == obj.eeg.NumberOfChannels) = [];
                                
                                [badchannels.rejectedChannels, h] = interactiveFigure(obj, EEG, output.noisyChannels.all, 0);
                                
                                saveas(h,fullfile(currFigOutput, [obj.ppNames{indPP} '_channels' num2str(blockCounter) '.png']))
                                
                                badchannels.rejectedChannels(badchannels.rejectedChannels == obj.eeg.NumberOfChannels+1) = [];
                                if ~isempty(badchannels.rejectedChannels)
                                    for indChan = 1:length(badchannels.rejectedChannels)
                                        fprintf(fileID,'\t- Bad channel: %s\n', num2str(badchannels.rejectedChannels(indChan)));
                                    end
                                    EEG = eeg_interp(EEG, badchannels.rejectedChannels, 'spherical');
                                else
                                    fprintf(fileID, '\t- No bad channels were detected\n');
                                end
                                
                            else
                                fprintf(fileID,'5) No bad channel selection applied, adviced to apply it later on the epoched data\n');
                            end
                            
                            % re-references to average
                            fprintf(fileID,'6) Re-reference the data to average signal.\n');
                            EEG.data(end+1,:) = 0;  % add reference channel
                            EEG.data = EEG.data-repmat(nanmean(EEG.data),[obj.eeg.NumberOfChannels+1,1]);
                            
                            save(currPrepross, 'EEG', 'findEvents', 'badchannels');
                        else
                            load(currPrepross, 'findEvents');
                            obj.event{indPP}{blockCounter} = findEvents;
                        end
                        blockCounter = blockCounter + 1;
                    end
                    close all;
                    fclose(fileID);
                end
                fprintf('All participants have been pre-processed.\n')
                fprintf('This cost the most time, so you are almost there!\n')
                fprintf('ADVICE: Check the EEG and optinal EOG and EMG figure manual.\n')
            end
        end
        
        function setFilters(obj)
            %% setFilters(obj)
            % sets up all filters ect. to be used on the data.
            
            % --------------  Low-pass filter     ------------------------
            % Here is the best low-pass filter for ERP analysis ever
            % - check freqz - It has a 3dB freq a bit lower than fc, at
            % about 35 Hz, and it has specially strong attenuation at exactly
            % 50Hz(mains in Ireland)/60 Hz(mains in USA).
            %
            % Get FIR filter weights 'LPK' for Hamming-windowed sinc
            % filter (see Semmlow 2004 signal proc book):
            
            wc  = obj.eeg.cuttoffLPF/obj.eeg.SampleRate*2*pi;
            LPK = nan(obj.eeg.LHamming,1);
            for i=1:obj.eeg.LHamming
                n = i-ceil(obj.eeg.LHamming/2);
                if n == 0, LPK(i) = wc/pi;
                else, LPK(i) = sin(wc*n)/(pi*n); end
            end
            obj.eeg.lowPassKernel = LPK.*hamming(obj.eeg.LHamming);
            
            % figure;
            % freqz(obj.eeg.lowPassKernel,1,10000,obj.eeg.SampleRate);
            
            %% --------------  Baseline -----------------------------------
            % define the baseline window:
            % not hard program
            if ~isfield(obj.eeg, 'baseline')
                % this is a round random epoch choice to use as a baseline.
                % However when you wanna work with a SSVEP you need a more
                % specific number from exactly two cylces of the SSVEP. e.g.
                % [-2*1000/obj.stim.freqSSVEP 0]. Therefore, just simple
                % overwrite the parameter.
                obj.eeg.baseline = [-0.1 0];
            end
        end
        
        function [rejectChannels, figHandle] = interactiveFigure(obj, EEG, rejectChannels, ManualCheck)
            data = EEG.data;
            figHandle = figure('units','normalized','outerposition',[0 0 1 1],'visible','off'); hold on;
            for indChan = 1:obj.eeg.NumberOfChannels-1
                plot((indChan-1)*100+data(indChan,:), 'k', 'tag',sprintf('Channel %d', indChan));
            end
            
            for indChan = rejectChannels
                plot((indChan-1)*100+data(indChan,:), 'r', 'tag',sprintf('Channel %d', indChan));
            end
            
            title('Bad channel selection - in RED bad channel selected on raw data')
            ylim([0-300 obj.eeg.NumberOfChannels*100+300]);
            yticks((0:obj.eeg.NumberOfChannels)*100);
            yticklabels(obj.eeg.ChannelsName)
            ax = gca;
            ax.FontSize = 8;
            xlim([0 size(data,2)])
            drawnow
            pause(1)
            
            % manual check if more than 10 channels are bad.
            if ManualCheck
                set(0, 'DefaultFigureVisible', 'on')
                fprintf('Manually check electrodes')
                global channel
                datacursormode on
                dcm = datacursormode(gcf);
                set(dcm,'UpdateFcn',@myupdatefcn)
                returnKey = 0;
                tmpChan   = 0;
                
                while returnKey == 0
                    pause(0.001)
                    if ~isempty(channel) && tmpChan == 0
                        tmpChan = channel;
                        [~, ~, kbCode] = KbCheck;
                        if find(kbCode) == 32 % add to rejectChannels backspace
                            rejectChannels = union(rejectChannels, channel);
                            plot((channel-1)*100+data(channel,:), 'r', 'tag',sprintf('Channel %d', channel));
                            kbCode = [];
                        elseif find(kbCode) == 78 % remove to rejectChannels backspace
                            rejectChannels(rejectChannels==channel) = [];
                            plot((channel-1)*100+data(channel,:), 'k', 'tag',sprintf('Channel %d', channel));
                            kbCode = [];
                        elseif find(kbCode) == 8 % remove to rejectChannels backspace
                            for indChanHer = 1:length(rejectChannels)
                                plot((rejectChannels(indChanHer)-1)*100+data(rejectChannels(indChanHer),:), 'k', 'tag',sprintf('Channel %d', rejectChannels(indChanHer)));
                            end
                            kbCode = [];
                            rejectChannels = [];
                        end
                    elseif tmpChan ~= channel
                        [~, ~, kbCode] = KbCheck;
                        if find(kbCode) == 32 % add to rejectChannels backspace
                            rejectChannels = union(rejectChannels, channel);
                            plot((channel-1)*100+data(channel,:), 'r', 'tag',sprintf('Channel %d', channel));
                            kbCode = [];
                        elseif find(kbCode) == 78 % remove to rejectChannels backspace
                            rejectChannels(rejectChannels==channel) = [];
                            plot((channel-1)*100+data(channel,:), 'k', 'tag',sprintf('Channel %d', channel));
                            kbCode = [];
                        elseif find(kbCode) == 8 % remove to rejectChannels backspace
                            for indChanHer = 1:length(rejectChannels)
                                plot((rejectChannels(indChanHer)-1)*100+data(rejectChannels(indChanHer),:), 'k', 'tag',sprintf('Channel %d', rejectChannels(indChanHer)));
                            end
                            kbCode = [];
                            rejectChannels = [];
                        end
                    else
                        [~, ~, kbCode] = KbCheck;
                        if find(kbCode) == 32 % add to rejectChannels backspace
                            rejectChannels = union(rejectChannels, channel);
                            plot((channel-1)*100+data(channel,:), 'r', 'tag',sprintf('Channel %d', channel));
                            kbCode = [];
                        elseif find(kbCode) == 78 % remove to rejectChannels backspace
                            rejectChannels(rejectChannels==channel) = [];
                            plot((channel-1)*100+data(channel,:), 'k', 'tag',sprintf('Channel %d', channel));
                            kbCode = [];
                        elseif find(kbCode) == 8 % remove to rejectChannels backspace
                            for indChanHer = 1:length(rejectChannels)
                                plot((rejectChannels(indChanHer)-1)*100+data(rejectChannels(indChanHer),:), 'k', 'tag',sprintf('Channel %d', rejectChannels(indChanHer)));
                            end
                            kbCode = [];
                            rejectChannels = [];
                        end
                    end
                    
                    if find(kbCode) == 13 % if enter renew the topoplots
                        rejectChannels(rejectChannels == 97) = [];
                        if ~isempty(rejectChannels)
                            %                             for indChan = 1:length(rejectChannels)
                            %                                 fprintf(fileID,'\t- Bad channel: %s\n', num2str(rejectChannels(indChan)));
                            %                             end
                            tmpEEG = eeg_interp(EEG, rejectChannels, 'spherical');
                            findNoisyChannels(tmpEEG);
                            
                            %                         else
                            %                             fprintf(fileID, '\t- No bad channels were detected\n');
                        end
                    end
                    
                    if find(kbCode) == 27
                        returnKey = 1;
                    end
                end
                
                clearvars -global channel
            end
        end
        
        function [EEG, figHandle] = applyFilters(obj, EEG, badChannels, events, plotFigure, filters)
            %% applyFilters
            % code will apply the pre-processing as set-up by the user.
            % This includes detrending the data, high-pass and low-pass
            % filtering. The preprocessed data can possibility be plotted.
            if ~exist('filters', 'var')
                filters = obj.eeg;
            end
            
            if isfield(badChannels, 'satChannels.Samples')
                Samples  = badChannels.satChannels.Samples;
                Channels = badChannels.satChannels.Channels;
            else
                Channels = [];
                Samples  = [];
            end
            
            %% %%%%%%%%%%%%%%%% detrend (DCCs)  %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.eeg.applyDetrend
                
                startT = max(events.Times(1)-obj.eeg.postDCC-EEG.srate, -obj.eeg.postDCC);   % go one second back from first stimulus trigger
                DetrendBoundaries = [startT badChannels.DCC(logical(badChannels.DCC>startT)) size(EEG.data,2)];
                % make a detrended version
                detrdata = EEG.data;
                for n=1:length(DetrendBoundaries)-1
                    interval = DetrendBoundaries(n)+obj.eeg.postDCC+1:DetrendBoundaries(n+1);
                    if ~isempty(interval)
                        detrdata(:,interval) = detrend(detrdata(:,interval)')';
                    end
                end
                
                % if there are channels that saturated, detrend the bits in between the saturated bits
                if ~isempty(Samples)
                    for q=1:length(Samples)
                        breaks = find(diff(Samples{q})>1);
                        sat_int_end = [Samples{q}(breaks) Samples{q}(end)];
                        sat_int_start = [Samples{q}(1) Samples{q}(breaks+1)];
                        interval = 1:sat_int_start(1)-1;
                        detrdata(Channels(q),interval) = detrend(EEG.data(Channels(q),interval)')';
                        for n=2:length(sat_int_start)
                            interval = sat_int_end(n-1)+1:sat_int_start(n)-1;
                            detrdata(Channels(q),interval) = detrend(EEG.data(Channels(q),interval)')';
                        end
                        if sat_int_end(end)~=size(EEG.data,2)
                            interval = sat_int_end(end)+1:size(EEG.data,2);
                            detrdata(Channels(q),interval) = detrend(EEG.data(Channels(q),interval)')';
                        end
                    end
                end
                
                if obj.eeg.applyDetrend==1
                    chans2detrend = (1:size(EEG.data,1))';
                elseif obj.eeg.applyDetrend==2
                    % only detrend the channels that benefit from it, in terms of a
                    % reduction of stdev by at least 25%:
                    chans2detrend = find(std(detrdata,[],2)<0.75*std(EEG.data,[],2));
                end
                EEG.data(chans2detrend,:) = detrdata(chans2detrend,:);
            end
            
            % High pass filter
            if filters.applyHPF
                if filters.simonHPF
                    [B,A] = butter(3,filters.HPFcutoff/(filters.SampleRate/2),'high');
                    startT = max(-filters.postDCC-EEG.srate,-filters.postDCC);   % go one second back from first stimulus trigger
                    HPFBoundaries = [startT size(EEG.data,2)];
                    
                    % first HPF the channels that didn't saturate:
                    for n=1:length(HPFBoundaries)-1
                        interval = HPFBoundaries(n)+filters.postDCC+1:HPFBoundaries(n+1);
                        if ~isempty(interval)
                            EEG.data(:,interval) = filtfilt(B,A,double(EEG.data(:,interval))')';
                        end
                    end
                else
                    EEG.data = eegfilt(EEG.data,obj.eeg.SampleRate,obj.eeg.HPFcutoff,0);
                end
            end
            
            
            % First LP Filter
            if obj.eeg.applyLPF
                %                 EEG.data = eegfilt(EEG.data,obj.eeg.SampleRate,0,obj.eeg.cuttoffLPF);
                
                for indChan = 1:size(EEG.data,1)
                    EEG.data(indChan,:) = conv(EEG.data(indChan,:),  filters.lowPassKernel, 'same');
                end
            end
            
            if plotFigure
                figHandle = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
                for indChan = 1:size(EEG.data,1)
                    plot((indChan-1)*100+EEG.data(indChan,:), 'k', 'tag',sprintf('Channel %d', indChan));
                end
                
                title('Preprocessed data')
                ylim([0-300 size(EEG.data,1)*100+300]);
                yticks([1 10:10:size(EEG.data,1) size(EEG.data,1)]*100);
                yticklabels([1 10:10:size(EEG.data,1) size(EEG.data,1)])
                xlim([0 size(EEG.data,2)])
                drawnow
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% -----------------   Epoching    --------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % collection of function that are used for the epoching.
        
        function epochingData(obj,redoEpoch)
            %% epochingData(obj,redoEpoch)
            % Function is used to extract the epoch data around the target
            % as well as around the reponse. For the target the behaviour
            % is as well extracted.
            
            if ~exist('redoEpoch', 'var'), redoEpoch = 0; end
            
            % Include part of the 'waiting' period before each trial. Initially we
            % extract a larger epoch to apply the ERP analysis on.
           
            if strcmpi(obj.stim.timing, 'past')
                obj.eeg.epoch = round((obj.stim.epoch(1) - max(obj.stim.lengthITI)) * obj.eeg.SampleRate):...
                    round(obj.stim.epoch(2)*obj.eeg.SampleRate);
            elseif strcmpi(obj.stim.timing, 'future')
                obj.eeg.epoch = round(obj.stim.epoch(1) * obj.eeg.SampleRate):...
                    round((obj.stim.epoch(2) + min(obj.stim.lengthITI)) * obj.eeg.SampleRate);
            elseif isempty(obj.stim.timing)
                obj.eeg.epoch = round(obj.stim.epoch(1)*obj.eeg.SampleRate):...
                    round(obj.stim.epoch(2)*obj.eeg.SampleRate);
            end
            
            obj.eeg.epochPlot = obj.eeg.epoch*1/obj.eeg.SampleRate;
            maxERPLength = length(obj.eeg.epoch);
            
            
            %% %%%%%%%%%%%%   Loop through participants to get ERP   %%%%%%
            for indPP = 1:length(obj.ppNames)
                
                currInput  = fullfile(obj.outputFolder, 'EEG data', obj.ppNames{indPP});
                currOutput = fullfile(obj.outputFolder, 'EEG data', obj.ppNames{indPP}, [obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.stim.timing '.mat']);
                
                if ~exist(currOutput, 'file') || redoEpoch
                    fprintf(['Epoching participant ' obj.ppNames{indPP} ])
                    
                    %% %%%%%%%%%%%%  Loop through block files   %%%%%%%
                    % Allocated information
                    
                    if obj.analysisEEG, ERP = nan(obj.eeg.NumberOfChannels, maxERPLength,  max(obj.numTrials*obj.numBlocks)); end % channels x timeseries (longest one 10000+1800+300) x epochs.
                    if obj.analysisEyelink, EYE = nan(2, maxERPLength,  max(obj.numTrials*obj.numBlocks)); end % only when eyelink data is measured will this fill up.
                    if obj.analysisEOG, EOG = nan(1, maxERPLength,  max(obj.numTrials*obj.numBlocks)); end% only when eyelink data is measured will this fill up.
                    if obj.analysisEMG, EMG = nan(2, maxERPLength,  max(obj.numTrials*obj.numBlocks)); end % only when eyelink data is measured will this fill up.
                    
                    tmpBehaviour = [];
                    tmpBehaviour.RT            = [];
                    tmpBehaviour.indReaction   = [];
                    tmpBehaviour.Reaction	   = [];
                    tmpBehaviour.Misses        = [];
                    tmpBehaviour.FalseAlarm    = [];
                    tmpBehaviour.indFalseAlarm = [];
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% --------------- Loop through block files   ---------
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    targCounter  = 1;	% accumulating across blocks
                    for indBlock = 1:obj.numBlocks(indPP)
                        fprintf('.')
                        if  obj.analysisEEG
                            EEG = [];
                            try
                                load(fullfile(currInput, [obj.ppNames{indPP} '_' num2str(indBlock) '_preprocessedEEG_HPF' num2str(obj.eeg.HPFcutoff) '.mat']), 'EEG');
                            catch
                                continue
                            end
                        end
                        
                        % load eyelink data if exist.
                        if obj.analysisEyelink && exist(fullfile(obj.outputFolder,'Eyelink data', obj.ppNames{indPP}, [obj.ppNames{indPP} '_' num2str(indBlock) '_preprocessedEyelink.mat']), 'file')
                            EL = [];
                            load(fullfile(obj.outputFolder,'Eyelink data', obj.ppNames{indPP}, [obj.ppNames{indPP} '_' num2str(indBlock) '_preprocessedEyelink' ]), 'EL');
                        end
                        
                        % load EMG if exist
                        if obj.analysisEOG && exist(fullfile(obj.outputFolder, 'EOG data', obj.ppNames{indPP}, [obj.ppNames{indPP} '_' num2str(indBlock) '_EOG.mat']), 'file')
                            tmpEOG = [];
                            tmpEOG = load(fullfile(obj.outputFolder, 'EOG data', obj.ppNames{indPP}, [obj.ppNames{indPP} '_' num2str(indBlock) '_EOG.mat']));
                            
                            tmpEOG = tmpEOG.EOG.rho;
                            %TODO get also frontal electrodes
                        end
                        
                        % load EMG if exist
                        if obj.analysisEMG && exist(fullfile(obj.outputFolder, 'EMG data', obj.ppNames{indPP}, [obj.ppNames{indPP} '_' num2str(indBlock) '_EMG.mat']), 'file')
                            tmpEMG = [];
                            tmpEMG = load(fullfile(obj.outputFolder, 'EMG data', obj.ppNames{indPP}, [obj.ppNames{indPP} '_' num2str(indBlock) '_EMG.mat']), 'EMG');
                            tmpEMG = tmpEMG.EMG.data;
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% --------------- Loop through events  ---------
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % index of event triggers to the events we
                        % want to time-lock to (targets).
                        
                        % Start creating the trialMatrix, a matrix that
                        % relates all epochs to the respective conditions.
                        
                        currCond = NaN(obj.numTrials, length(obj.condNames));
                        for indCond = 1:length(obj.conditions)
                            try
                                tmpCond  = eval(sprintf('obj.experiment{%i}{%i}.%s', indPP, indBlock, obj.condNames{indCond}));
                            catch
                                keyboard
                            end
                            if length(tmpCond) == 1
                                currCond(:,indCond) = repmat(tmpCond,1, obj.numTrials);
                            else
                                currCond(:,indCond) = tmpCond;
                            end
                        end
                        
                        trialITI = eval(sprintf('obj.experiment{%i}{%i}.%s', indPP, indBlock, obj.stim.namesITI));
                        
                        if obj.analysisEEG
                            
                            targets = find(ismember(obj.event{indPP}{indBlock}.Number, obj.eeg.epochLock));
                            currTimes = obj.event{indPP}{indBlock}.Times(targets);
                            if length(targets) ~= obj.numTrials
                                keyboard
                                %                             if length(targets) == obj.numTrials +1
                                %                                 targets(1) = [];
                                %                             elseif length(targets) > obj.numTrials+1 && length(find(obj.event{indPP}{indBlock}.Number == obj.triggers.start)) > 1
                                %                                 index = find(obj.event{indPP}{indBlock}.Number == obj.triggers.start);
                                %                                 targets = targets(targets > index(end));
                                %                                 currTimes = obj.event{indPP}{indBlock}.Times(targets);
                            elseif isempty(targets)
                                keyboard
                                targets   = 1:length(obj.experiment{indPP}{indBlock}.RespLR) - 1;
                                currTimes = round((obj.experiment{indPP}{indBlock}.RespT - obj.experiment{indPP}{indBlock}.RespT(1)) .* obj.eeg.SampleRate); % in ms --> convert to samples
                                currTimes(1) = [];
                            end
                            
                            if strcmp(obj.stim.timing, 'future')
                                targets(end) = [];
                            end
                        else
                            targets = obj.experiment{indPP}{indBlock}.TargOnT;
                        end
                        
                        for indTarget = 1:length(targets)
                            %% %%%%%%%%%%%  Epoch definition     %%%%%%%%%%%%%%%%%%%%%%%%%%
                            % note the use of integer number of cycles of SSVEP
                            % hence get the list of sample points relative to a given
                            % event that we need to extract from the continuous EEG
                            
                            % Keep track of the task condition for each trial
                            tmpBehaviour.trialMatrix(targCounter,:) = currCond(indTarget,:);
                            
                            if obj.analysisEEG
                                if isa(obj.stim.epoch ,'function_handle')
                                    currEpoch = obj.stim.epoch(indPP, indBlock, indTarget);
                                else
                                    currEpoch = obj.stim.epoch;
                                end
                                
                                if strcmpi(obj.stim.timing, 'past')
                                    currITI = trialITI(indTarget);
                                    indCurrEpoch = round((0-(currITI-0.5))*obj.eeg.SampleRate):...
                                        round(currEpoch(2)*obj.eeg.SampleRate);
                                elseif strcmpi(obj.stim.timing, 'future')
                                    currITI = trialITI(indTarget+1);
                                    indCurrEpoch = round((currEpoch(1))*obj.eeg.SampleRate):...
                                        round((currEpoch(2)+currITI)*obj.eeg.SampleRate);
                                elseif isempty(obj.stim.timing)
                                    currITI = trialITI(indTarget);
                                    indCurrEpoch = round((currEpoch(1))*obj.eeg.SampleRate):...
                                        round(currEpoch(2)*obj.eeg.SampleRate);
                                end
                                
                                
                                % 1) check whether the epoch relative to
                                % the current event is within the bounds of
                                % the EEG data
                                if (currTimes(indTarget) + indCurrEpoch(1) <= 0) || (currTimes(indTarget) + indCurrEpoch(end) > size(EEG.data,2))
                                    tmpBehaviour.RT(targCounter,:)            = NaN;
                                    tmpBehaviour.indReaction(targCounter,1)   = NaN;
                                    tmpBehaviour.Reaction(targCounter,:)	  = NaN;
                                    tmpBehaviour.Misses(targCounter,1)        = 0;
                                    tmpBehaviour.FalseAlarm(targCounter,:)    = 0;
                                    tmpBehaviour.indFalseAlarm(targCounter,:) = NaN;
                                    targCounter = targCounter + 1;
                                    continue
                                end
                                
                                % 2) Now extract the epoch
                                epoch = EEG.data(1:obj.eeg.NumberOfChannels, currTimes(indTarget) + indCurrEpoch);
                                
                                % now add the current epoch onto the growing
                                %'erp' matrix by concatenation along the
                                % 3rd dimension.
                                if strcmpi(obj.stim.timing, 'past')
                                    ERP(:, end-length(indCurrEpoch)+1:end, targCounter) = epoch;
                                elseif strcmpi(obj.stim.timing, 'future')
                                    ERP(:, 1:length(indCurrEpoch), targCounter) = epoch;
                                elseif isempty(obj.stim.timing)
                                    ERP(:,1:length(indCurrEpoch),targCounter) = epoch;
                                end
                            end
                            
                            if obj.analysisEyelink && ~isempty(EL) && round(EL.timeSeries(end)) > EEG.times(end)
                                
                                timePoints = EEG.times(:,currTimes(indTarget)) + indCurrEpoch([1 end]);
                                
                                [~, indEvent(1)] = min(abs(EL.timeSeries - timePoints(1)));
                                [~, indEvent(2)] = min(abs(EL.timeSeries - timePoints(2)));
                                
                                if diff(indEvent) + 1 >  length(epoch)
                                    indEvent(2) = indEvent(2) - ((diff(indEvent) + 1) - length(epoch));
                                elseif diff(indEvent) + 1 <  length(epoch)
                                    indEvent(2) = indEvent(2) + (length(epoch) - (diff(indEvent) + 1));
                                end
                                
                                if indEvent(2) < length(EL.indSacc)
                                    epochEL = [EL.indBlinks(indEvent(1):indEvent(2)); EL.indSacc(indEvent(1):indEvent(2))];
                                else
                                    epochEL = nan;
                                end
                                
                                if strcmpi(obj.stim.timing, 'past')
                                    EYE(:, end-length(indCurrEpoch)+1:end, targCounter) = epoch;
                                elseif strcmpi(obj.stim.timing, 'future')
                                    EYE(:, 1:length(indCurrEpoch), targCounter) = epoch;
                                elseif isempty(obj.stim.timing)
                                    EYE(:,1:length(indCurrEpoch),targCounter) = epoch;
                                end
                            end
                            
                            if obj.analysisEOG
                                % 2) Now extract the epoch
                                epoch = tmpEOG(:, currTimes(indTarget) + indCurrEpoch);
                                
                                if strcmpi(obj.stim.timing, 'past')
                                    EOG(:, end-length(indCurrEpoch)+1:end, targCounter) = epoch;
                                elseif strcmpi(obj.stim.timing, 'future')
                                    EOG(:, 1:length(indCurrEpoch), targCounter) = epoch;
                                elseif isempty(obj.stim.timing)
                                    EOG(:,1:length(indCurrEpoch),targCounter) = epoch;
                                end
                                
                            end
                            
                            if obj.analysisEMG
                                % 2) Now extract the epoch
                                epoch = tmpEMG(:, currTimes(indTarget) + indCurrEpoch);
                                if strcmpi(obj.stim.timing, 'past')
                                    EMG(:, end-length(indCurrEpoch)+1:end, targCounter) = epoch;
                                elseif strcmpi(obj.stim.timing, 'future')
                                    EMG(:, 1:length(indCurrEpoch), targCounter) = epoch;
                                elseif isempty(obj.stim.timing)
                                    EMG(:,1:length(indCurrEpoch),targCounter) = epoch;
                                end
                            end
                            
                            % get the response - first one after the cue, may include jumpguns
                            nextresp = find(obj.experiment{indPP}{indBlock}.RespT   > obj.experiment{indPP}{indBlock}.TargOnT(indTarget)...
                                & obj.experiment{indPP}{indBlock}.RespT < obj.experiment{indPP}{indBlock}.TargOnT(indTarget) + obj.stim.RTdeadLine(2)); % 19/4/18 (had been +1.6) pick up any button click within 200 ms of offset of motion
                            
                            numResp = max([length(nextresp), 1, size(tmpBehaviour.Reaction,2)]);
                            tmpBehaviour.RT(targCounter,1:numResp)          = NaN;
                            tmpBehaviour.indReaction(targCounter,1:numResp) = NaN;
                            tmpBehaviour.Reaction(targCounter,1:numResp)	= 0;
                            tmpBehaviour.Misses(targCounter,1) = 0;
                            tmpBehaviour.ITI(targCounter,1)  = currITI;
                            
                            if ~isempty(nextresp)
                                for indResp = 1:length(nextresp)
                                    
                                    tmpBehaviour.Reaction(targCounter,indResp) = obj.experiment{indPP}{indBlock}.RespLR(nextresp(indResp));
                                    
                                    tmpBehaviour.RT(targCounter,indResp) = obj.experiment{indPP}{indBlock}.RespT(nextresp(indResp)) - obj.experiment{indPP}{indBlock}.TargOnT(indTarget); % in ms!
                                    if obj.analysisEEG
                                        [~, tmpBehaviour.indReaction(targCounter,indResp)] = min(abs(obj.eeg.epochPlot - tmpBehaviour.RT(targCounter,indResp)));
                                    end
                                end
                            else
                                tmpBehaviour.Misses(targCounter,1)  = 1;
                            end
                            
                            % find False Alarms.
                            if strcmpi(obj.stim.timing, 'future')
                                FAresp = find((obj.experiment{indPP}{indBlock}.RespT < obj.experiment{indPP}{indBlock}.TargOnT(indTarget)...
                                    & obj.experiment{indPP}{indBlock}.RespT > obj.experiment{indPP}{indBlock}.TargOnT(indTarget) -  obj.eeg.epochPlot(1)) | ...
                                    (obj.experiment{indPP}{indBlock}.RespT > obj.experiment{indPP}{indBlock}.TargOnT(indTarget) + obj.stim.RTdeadLine(2)...
                                    & obj.experiment{indPP}{indBlock}.RespT < obj.experiment{indPP}{indBlock}.TargOnT(indTarget) + obj.stim.RTdeadLine(2) + currITI));
                            else
                                if exist('currITI', 'var')
                                    FAresp = find(obj.experiment{indPP}{indBlock}.RespT <= obj.experiment{indPP}{indBlock}.TargOnT(indTarget)...
                                        & obj.experiment{indPP}{indBlock}.RespT > obj.experiment{indPP}{indBlock}.TargOnT(indTarget) - (currITI - (obj.stim.RTdeadLine(2) - obj.stim.duration)));
                                else
                                    FAresp = [];
                                end
                            end
                            
                            numFA = max([length(FAresp), 1, size(tmpBehaviour.FalseAlarm,2)]);
                            tmpBehaviour.FalseAlarm(targCounter,1:numFA)    = 0;
                            tmpBehaviour.indFalseAlarm(targCounter,1:numFA) = NaN;
                            
                            if ~isempty(FAresp)
                                for indFA = 1:length(FAresp)
                                    tmpBehaviour.FalseAlarm(targCounter,indFA) = 1;
                                    tmpBehaviour.indFalseAlarm(targCounter,indFA) = obj.experiment{indPP}{indBlock}.RespT(FAresp(indFA)) - obj.experiment{indPP}{indBlock}.TargOnT(indTarget);     % in ms!
                                end
                            end
                            
                            targCounter = targCounter + 1;
                        end
                    end
                    
                    saveFiles = {'tmpBehaviour'};
                    if obj.analysisEEG
                        ERP = ERP(:,:,1:size(tmpBehaviour.indReaction,1));
                        
                        % artifact rejection
                        [tmpBehaviour.artifacts, tmpBehaviour.blinks, ERP] = artifactReject(obj, ERP, EOG, tmpBehaviour.RT);
                        saveFiles = {saveFiles{:}, 'ERP'};
                        %{
                        if  obj.eeg.applyCSD
                            goodTrials = tmpBehaviour.artifacts & tmpBehaviour.blinks;%
                            tic
                            % load matrices that CSD function needs to run for the specific biosemi 128 cap:
                            normalcsdERP(1:obj.eeg.NumberOfChannels,:,goodTrials) = CSD(ERP(1:obj.eeg.NumberOfChannels,:,goodTrials),...
                                obj.eeg.transChanlocs.G,obj.eeg.transChanlocs.H);
                            timing.normalCSD(indPP) = toc
                            
                             % load matrices that CSD function needs to run for the specific biosemi 128 cap:
                             tic
                            csdERP(1:obj.eeg.NumberOfChannels,:,goodTrials) = FastCSD(ERP(1:obj.eeg.NumberOfChannels,:,goodTrials),...
                                obj.eeg.transChanlocs.G,obj.eeg.transChanlocs.H);
                            timing.fastCSD(indPP) = toc
                            saveFiles = {saveFiles{:}, 'csdERP'};
                        end
                        %}
                    end
                    % save trials
                    
                    if obj.analysisEyelink && ~isempty(EL)
                        saveFiles = {saveFiles{:}, 'EYE'};
                    end
                    
                    if obj.analysisEMG
                        saveFiles = {saveFiles{:}, 'EMG'};
                    end
                    
                    if obj.analysisEOG
                        saveFiles = {saveFiles{:}, 'EOG'};
                    end
                    
                    save(currOutput, saveFiles{:}, '-v7.3');
                    obj.behaviour{indPP} = tmpBehaviour;
                    
                    fprintf('\n')
                else
                    load(currOutput, 'tmpBehaviour');
                    
                    obj.behaviour{indPP} = tmpBehaviour;
                    
                end
                
            end
        end
        
        function [artifact, blinks, ERP] = artifactReject(obj, ERP, EOG, RT)
            %% %%%%%%%%%%%%%% Artifact detection    %%%%%%%%%%%%%%%%%
            % quick function to determine epochs with artifact to be
            % removed later onwards for CSD or for general data used in
            % sortERPs.
            fprintf('Artifact rejection\n');
            
            % 1) Check quality of the epoch using both possible
            % Eyelink data, as well as a artifact EEG data.
            % Preallocated the good trials.
            
            baseline = obj.eeg.epochPlot >= obj.eeg.baseline(1) & obj.eeg.epochPlot < obj.eeg.baseline(2);
            tmpERP = ERP - repmat(nanmean(ERP(:,baseline,:),2),1,size(ERP,2),1);
            EOG = EOG - repmat(nanmean(EOG(:,baseline,:),2),1,size(EOG,2),1);
            
            trangeTarget = obj.eeg.epochPlot >= obj.eeg.targetEpoch(1) & obj.eeg.epochPlot <= obj.eeg.targetEpoch(end);
            artifact = squeeze(max(abs(tmpERP(:,trangeTarget,:)),[],2) >= obj.eeg.artifactThres);
            if obj.eeg.intArtifact
                load('blankEEG.mat')
            end
            
            for indEpoch = 1:size(ERP,3)
                if obj.eeg.intArtifact && any(artifact(:,indEpoch)) && nanmean(artifact(:,indEpoch)) < 0.1
                    EEG.data = ERP(:,:,indEpoch);
                    EEG.chanlocs = obj.eeg.chanlocs;
                    EEG.nbchan   = obj.eeg.NumberOfChannels;
                    tmp = eeg_interp(EEG, find(artifact(:,indEpoch)), 'spherical');
                    ERP(:,:,indEpoch) = tmp.data;
                end
                
                trangeTarget = obj.eeg.epochPlot >= -100 & obj.eeg.epochPlot <= RT(indEpoch,1)+0.2;
                
                if obj.analysisEyelink
                    blinks(indEpoch) = squeeze(~any(EOG(1, trangeTarget, indEpoch)));
                elseif obj.analysisEOG
                    blinks(indEpoch) = squeeze(~any(EOG(:,trangeTarget,indEpoch) > obj.eeg.artifactEOG));
                end
            end
            if obj.eeg.intArtifact
                artifact = nanmean(artifact) <= 0.1;
            else
                artifact = ~any(artifact);
            end
            
        end
        
        function applyCSD(obj)
            
            if obj.eeg.applyCSD
                fprintf('Now apply CSD')
                for indPP = 1:length(obj.ppNames)
                    fprintf('.')
                    goodTrials = obj.behaviour{indPP}.artifacts & obj.behaviour{indPP}.blinks;%
                    
                    currOutput = fullfile(obj.outputFolder, 'EEG data', obj.ppNames{indPP}, [obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.stim.timing '.mat']);
                    
                    m = matfile(currOutput,'Writable',true);
                    varlist = whos(m); names = {varlist.name};
                    
                    if ~any(contains(names, 'csdERP'))
                        ERP = m.ERP; csdERP = NaN(size( m.ERP));
                        
                        % load matrices that CSD function needs to run for the specific biosemi 128 cap:
                        csdERP(1:obj.eeg.NumberOfChannels,:,goodTrials) = CSD(ERP(1:obj.eeg.NumberOfChannels,:,goodTrials),...
                            obj.eeg.transChanlocs.G,obj.eeg.transChanlocs.H);
                        
                        m.csdERP = csdERP;
                    end
                    clear m ERP csdERP
                end
                
                fprintf('\n');
            end
        end
        
        function binRTs(obj, firstRT, qps, theseCond)
            %% function binRTs
            % We sorted trials according to RT and divided them into two
            % equal-sized bins. Importantly, RT binning was done
            % within each trial condition (i.e., intertrial interval
            % and target contrast direction), thus eliminating confounding
            % factors known to have an influence on RT.
            
            if ~exist('firstRT', 'var'); firstRT = 0; end
            if ~exist('qps', 'var');  qps = 0.5; end
            if ~exist('theseCond', 'var'); theseCond = 1:length(obj.conditions); end
            
            for indPP = 1:length(obj.ppNames)
                obj.behaviour{indPP}.RT(obj.behaviour{indPP}.RT == 0) = NaN;
                % remove all the 'early' response
                indEarly = obj.behaviour{indPP}.RT <= obj.stim.RTCutOff;
                %obj.behaviour{indPP}.Early = zeros(size(indEarly,1),1);
                
                if any(indEarly(:))
                     
                    % move early response to a miss as well as adding them
                    % to false alarms.
                    [r,c] = find(indEarly);
                    early = obj.behaviour{indPP}.RT(find(indEarly));
                    %
                    numFA = size(obj.behaviour{indPP}.FalseAlarm,2);
                    obj.behaviour{indPP}.FalseAlarm(:,end+1:end+max(c))= 0;
                    obj.behaviour{indPP}.indFalseAlarm(:,end+1:end+max(c)) = NaN;
                    
                    for currEarly = 1:length(r)
                        obj.behaviour{indPP}.FalseAlarm(r(currEarly),numFA+c(currEarly))= 1;
                        obj.behaviour{indPP}.indFalseAlarm(r(currEarly),numFA+c(currEarly)) = early(currEarly);
                    end
                    
                    % move early response to a miss as well as adding them
                    % to false alarms.
                    obj.behaviour{indPP}.RT(find(indEarly)) = NaN;
                    obj.behaviour{indPP}.indReaction(find(indEarly)) = NaN;
                    obj.behaviour{indPP}.Reaction(find(indEarly)) = NaN;
                    obj.behaviour{indPP}.Misses(any(indEarly,2)) = 1; 
                 
                end
                
                if firstRT
                    % only get the first response that is larger then
                    % the RT cutoff
                    for indTrial = find(~obj.behaviour{indPP}.Misses(:,1))'
                        try
                            if  ~isempty(find(~isnan(obj.behaviour{indPP}.RT(indTrial,:)),1))
                                obj.behaviour{indPP}.RT(indTrial,1) = obj.behaviour{indPP}.RT(indTrial, find(~isnan(obj.behaviour{indPP}.RT(indTrial,:)),1));
                            end
                        catch
                            keyboard
                        end
                    end
                    obj.behaviour{indPP}.RT(:,2:end) = [];
                    obj.behaviour{indPP}.indReaction(:,2:end) = [];
                    obj.behaviour{indPP}.Reaction(:,2:end) = [];
                end
                
                obj.behaviour{indPP}.binRT = nan(size(obj.behaviour{indPP}.RT));
                
                posCombination = unique(obj.behaviour{indPP}.trialMatrix(:,theseCond), 'rows');
                
                for indCond = 1:length(posCombination)
                    currCond = find(sum(obj.behaviour{indPP}.trialMatrix(:,theseCond) == posCombination(indCond,:),2) == size(posCombination, 2));
                    
                    % remove all conditions that didn't have a
                    % response, e.g. the nans.
                    currRT = obj.behaviour{indPP}.RT(currCond,:);
                    Order = discretize(currRT, [0 quantile(currRT,qps)-0.5/60 obj.stim.RTdeadLine(end)], 'IncludedEdge','right');
                    obj.behaviour{indPP}.binRT(currCond,:) = Order; % bin data in order to get average RT quantiles
                    
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% -----------------	Statistics     ----------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [rm, results, tbl] = doRANOVA(obj, plotThis, plotComb, TimeOnTask, TimeOnExperiment)
            %% plotRANOVA(obj)
            % plot behaviour wise, e.g. general experiment and RT
            % distributions.
            %   plotThis:   plots the results or the depending variables on the
            %               measurements (e.g. Hits, FA and reaction time)
            %               when set to 1.
            %   plotComb:   gives which variable should be plotted on the
            %               x-axis, as grouping or as seperated plots.
            %
            % Additional it allows you to check for certain 'time' effects:
            %   TimeOnTask: How long is a participant in the block (e.g.
            %               comparing beginning of the block to the end.
            %   TimeOnExperiment:   How long into the experiment is the
            %                       participant, e.g. practice
            %                       effect/fatigue ect.
            
            close all
            
            % Do we want to include any 'partice effect'
            if ~exist('TimeOnTask', 'var'), TimeOnTask = 0; end
            if ~exist('TimeOnExperiment', 'var'), TimeOnExperiment = 0; end
            
            %% set up model parameters
            % as there are different parameters that you might wanna add,
            % like time dependency. For now it just simple time on task and
            % time within the experiment.
            
            [within, between, modelFun, allTrialMatrix] = getConditions(obj, plotComb, TimeOnTask, TimeOnExperiment);
            % we want to formally test if there are effects of Evidence stength and
            % context while accounting for participants differences as well as
            % Reaction time effects
            
            for indCond = 1:length(plotComb)
                eval(sprintf('tbl.%s =  reshape((allTrialMatrix(:, indCond, :)), [],1);',  within.Name{indCond}));
            end
            
            % extracted the means per participant per conditions
            reactionTimes = nan(length(obj.ppNames), within.numCond);
            hitRate       = nan(length(obj.ppNames), within.numCond);
            faRate        = nan(length(obj.ppNames), within.numCond);
            Misses        = nan(length(obj.ppNames), within.numCond);
            Accuracy      = nan(length(obj.ppNames), within.numCond);
            Errors        = nan(length(obj.ppNames), within.numCond);
            
            tbl.ppNames   = [];
            tbl.FAcount   = [];
            tbl.Accuracy  = [];
            
            tbl.Misses    = [];
            tbl.RT        = [];
            tbl.early     = [];
            
            for indPP = 1:length(obj.ppNames)
                obj.behaviour{indPP}.indFalseAlarm(obj.behaviour{indPP}.indFalseAlarm == 0) = NaN;
                posCombination = unique(allTrialMatrix(:,:, indPP), 'rows');
                posCombination(any(isnan(posCombination),2),:) = [];
                tmpAccuracy = nan(size(allTrialMatrix,1),1);
                tmpFAcount  = nan(size(allTrialMatrix,1),1);
                tmpMisses = nan(size(allTrialMatrix,1),1);
                tmpRT = nan(size(allTrialMatrix,1),1);
                tmpEarly = nan(size(allTrialMatrix,1),1);
                
                for indCond  = 1:size(posCombination,1)
                    currCond = sum(allTrialMatrix(:,:,indPP) == posCombination(indCond,:),2) == size(posCombination, 2);
                    
                    % remove all conditions that didn't have a
                    % response, e.g. the nans.
                    currRT      = obj.behaviour{indPP}.RT(currCond,:);
                    currMisses  = obj.behaviour{indPP}.Misses(currCond);
                    
                    currindFA   = isnan(obj.behaviour{indPP}.indFalseAlarm(currCond,:));
                    currFA      = obj.behaviour{indPP}.FalseAlarm(currCond,:);
                    currFA(currindFA) = 0;
                    
                    currEarly   = obj.behaviour{indPP}.Early(currCond,:);
                    
                    if ~isempty(obj.stim.FACutOff)
                        currFA = double(obj.behaviour{indPP}.indFalseAlarm(currCond,:) >= obj.stim.FACutOff);
                    end
                    
                    if obj.DetectOrDisc
                        currAccuracy = obj.behaviour{indPP}.Reaction(currCond,:);
                        tmpHit       = currAccuracy ~= 0 & ~isnan(currAccuracy);
                        if any(sum(tmpHit,2) > 1)
                            tmpHit(sum(tmpHit,2) > 1,2:end) = 0;
                        end
                        currErrors   = double(~(currAccuracy == obj.behaviour{indPP}.trialMatrix(currCond,1)));
                        currAccuracy = double(currAccuracy   == obj.behaviour{indPP}.trialMatrix(currCond,1));
                        
                        currAccuracy(isnan(currRT)) = nan;
                        currErrors(isnan(currRT)) = nan;
                        
                        Accuracy(indPP,indCond) = nanmean(currAccuracy(tmpHit))*100;
                        Errors(indPP,indCond)   = nanmean(currErrors(tmpHit))*100;
                    else
                        currAccuracy = obj.behaviour{indPP}.Reaction(currCond,1);
                        currAccuracy(isnan(currAccuracy)) = 0;
                        
                        % get signal detection information.
                        hitRate(indPP,indCond) = nanmean(currAccuracy)*100;
                    end
                    
                    % get mean reaction times
                    reactionTimes(indPP,indCond) = nanmedian(currRT);
                    
                    % Misses are in common.
                    Misses(indPP,indCond) = nanmean(currMisses)*100;
                    
                    tmpMisses(currCond)   = currMisses;
                    tmpEarly(currCond)    = currEarly;
                    
                    tmpAccuracy(currCond) = currAccuracy;
                    
                    tmpFAcount(currCond)  = nansum(currFA,2);
                    tmpRT(currCond)       = currRT;
                    
                    faRate(indPP,indCond) = sum(currFA(:));%/length(currCond).*100;
                end
                tbl.early    = [tbl.early; tmpEarly];
                tbl.RT       = [tbl.RT; nanzscore(tmpRT)];
                tbl.Misses   = [tbl.Misses; tmpMisses];
                tbl.Accuracy = [tbl.Accuracy; tmpAccuracy];
                tbl.FAcount  = [tbl.FAcount; tmpFAcount];%./(obj.behaviour{indPP}.ITI/min(obj.behaviour{indPP}.ITI))];
                
                tbl.ppNames  = [tbl.ppNames; repmat(indPP, size(tmpFAcount))];
            end
            
            tbl.ppNames = nominal(tbl.ppNames);
            tbl = struct2table(tbl);
            
            
            strCond = 'Conditions: ';
            for indCond = 1:size(within.Name,2)
                strCond = strcat(strCond, [within.Name{indCond} ', ']);
                eval(sprintf('within.Table.%s = categorical(within.Table.%s);', within.Name{indCond}, within.Name{indCond}));
            end
            
            if sum(~within.Par)
                for indCond = 1:size(between.Name,2)
                    strCond = strcat(strCond, ['Groups:' between.Name{indCond} ', ']);
                end
            end
            
            %% run the ranovas on the experiment.
            % for log files
            fileID = fopen(fullfile(obj.logFolder, 'statistics.log'),'w');
            fprintf(fileID, '-------------- Statistical analysis --------------\n');
            fprintf(fileID, 'comparing reaction time, and accuracy variables between the different conditions\n');
            fprintf(fileID, 'An repeated-measurement ANOVA was applied to see if \n');
            fprintf(fileID, 'the reaction times for the different conditions come from the same distribution.\n');
            
            % ------------- Reaction time ---------------------------------
            [rm.RT, results.RT] = preformRANOVA(obj, fileID, 'REACTION TIME', strCond, reactionTimes, between, within, modelFun);
            if obj.DetectOrDisc
                % ------------- Correct response   ----------------------------------
                [rm.Accuracy, results.Accuracy] = preformRANOVA(obj, fileID, 'ACCURACY', strCond, Accuracy, between, within, modelFun);
                
                % ------------- Mistakes   ----------------------------------
                [rm.Error, results.Error] = preformRANOVA(obj, fileID, 'ERRORS', strCond, Errors, between, within, modelFun);
            else
                % ------------- Hit rate     ----------------------------------
                [rm.Hit, results.Hit] = preformRANOVA(obj, fileID, 'HIT RATE', strCond, hitRate, between, within, modelFun);
            end
            
            % ------------- False alarms ----------------------------------
            [rm.FA, results.FA] = preformRANOVA(obj, fileID, 'FALSE ALARMS', strCond, faRate, between, within, modelFun);
            
            % ------------- Misses ----------------------------------
            [rm.Misses, results.Misses] = preformRANOVA(obj, fileID, 'MISSES', strCond, Misses, between, within, modelFun);
            
            fclose(fileID);
            
            if plotThis
                %% plotting of the different variables.
                % here we plot the depending parameters as a function of
                % the independent variables (y). How this is plotted depends on
                % the plotComb. e.g. which are the x-axis, which are
                % plotted as group parameters and which are plotted as
                % seperated parameters.
                
                % ------------------- reaction time -----------------------
                plotAverageBehaviour(obj, 'RT (s)', rm.RT, plotComb);
                
                if obj.DetectOrDisc
                    % ------------------- Accuracy -----------------------
                    plotAverageBehaviour(obj, 'Response accuracy (%)', rm.Accuracy, plotComb);
                    
                    % ------------------- Errors -----------------------
                    plotAverageBehaviour(obj, 'Errors (%)', rm.Error, plotComb)
                else
                    % ------------------- Hit rate -----------------------
                    plotAverageBehaviour(obj, 'Hit rate (%)', rm.Hit, plotComb)
                end
                
                % ------------------- FA rate -----------------------
                plotAverageBehaviour(obj, '# False alarm rate', rm.FA, plotComb)
                
                % ------------------- Misses -----------------------
                plotAverageBehaviour(obj, 'Misses (%)', rm.Misses, plotComb)
            end
        end
        
        function [rm, results] = preformRANOVA(obj, fileID, titleName, strCond, input, between, within, modelFun)
            % repeating function to perform the RANOVA for difference
            % behavioural parameters.
            
            % write away the descriptives statistics
            fprintf(fileID, '----- Descriptive Statistics %s:\n', titleName);
            fprintf(fileID, [strCond(1:end-1) '\n']);
            
            for indGroup = unique(between.Table)
                for indVar = 1:within.numCond
                    currMean = round(nanmean(input(between.Table == indGroup,indVar)));
                    currStd  = round(nanstd(input(between.Table == indGroup,indVar)));
                    fprintf(fileID, '%s: Condition %s, group %s: mu: %i, sigma: %i\n',...
                        titleName, num2str(within.Condition(indVar,:)), num2str(indGroup),...
                        currMean, currStd);
                end
            end
            
            fprintf(fileID, '-----  Statistics testing %s:\n', titleName);
            inputTable = array2table(input);
            
            if sum(~within.Par)
                eval(sprintf('inputTable.%s = categorical(betweenTable'');', between.Name{1}));
                rm = fitrm(inputTable, sprintf('input1-input%i ~ %s', size(inputTable,2), between.Name{1}), 'WithinDesign', within.Table);
            else
                rm = fitrm(inputTable, sprintf('input1-input%i ~ 1', size(inputTable,2)), 'WithinDesign', within.Table);
            end
            
            spherTest = mauchly(rm);
            results   = ranova(rm, 'WithinModel', modelFun);
            
            % check for test for sphericity
            results.pValueHF = [];
            results.pValueLB = [];
            if spherTest.pValue < 0.05
                fprintf(fileID, 'Mauchlys Test of Sphericity indicated that the assumption of sphericity had been violated (X^2(%d) = %3.2f, p = %0.4f).\n', ...
                    spherTest.DF, spherTest.ChiStat, spherTest.pValue);
                fprintf(fileID, 'From here onwards the Greenhouse-Geisser corrected p-values are reported to compare %s.\n', titleName);
                pvalues = results.pValueGG(:);
                results.pValue = pvalues;
            else
                fprintf(fileID, 'Mauchlys Test of Sphericity indicated that the assumption of sphericity had not been violated (X^2(%d) = %3.2f, p = %0.4f).\n', ...
                    spherTest.DF, spherTest.ChiStat, spherTest.pValue);
                % pvalues = results.pValue(:);
            end
            results.pValueGG = [];
            
            % report sigficance
            for indPVal = find(~contains(results.Properties.RowNames, 'Error'))'
                if results.pValue(indPVal) > 0.05
                    fprintf(fileID, 'There is no significant effect of %s on %s (F(%d, %d) = %3.2f, p = %0.4f).\n', ...
                        results.Properties.RowNames{indPVal}, titleName, results.DF(indPVal), results.DF(indPVal+1), results.F(indPVal), results.pValue(indPVal));
                else
                    fprintf(fileID, 'There is a significant effect of %s on %s (F(%d, %d) = %3.2f, p = %0.4f).\n', ...
                        results.Properties.RowNames{indPVal}, titleName, results.DF(indPVal), results.DF(indPVal+1), results.F(indPVal), results.pValue(indPVal));
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% -----------------    plot ERPS    -------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract and plot ERPs.
        function [Quality, fig] = plotERP(obj, methodUsed, plotName, Elec, BaselineCorrect, plotThis, plotComb,...
                grouping, trangeTopo, TargetOrResponse, MorletOrsFFT, freqRange)
            %% [Quality, figInfo, tbl] = plotERP(obj, plotName, Elec, plotThis, plotComb, grouping, trangeTopo,
            %       TargetOrResponse, MorletOrsFFT, freqRange)
            % function will extract topoplot to get best electrodes, if
            % Elec is not defined. Topoplot uses trangeTopo to identify the
            % range and you can either choose TargetOrResponse.
            % It further plots target and response-locked
            % average plot using the selected electrodes. Not this
            % standardly baseline-correct the data!
            % Parameters:
            %   MethodUsed          =  gives the possibility to just use
            %                          1) the pre-set electrodes
            %                          2) get the 3 best electrodes with highest SNR
            %                          3) Lucalize a cluster of electrodes
            %   plotName            =  given name to the plots (e.g. CPP)
            %   Elec                =  selected electrodes (either number
            %                          or names), when empty the electrodes
            %                          will be asked to be selected during
            %                          the topoplot phase.
            %   BaselineCorrect     =   0 - no, 1 - yes
            %   plotThis            =  	0) only topoplots
            %                           1) individual plots,
            %                           2) average Hits and seperated average Misses and FA
            %                           3) ,, include all average target-locked Misses
            %                           4) plot target and response seperately.
            %                           5) ,, include all average target-locked
            %                               Misses  (NOTE this doesn't change the
            %                               response-locked ones).
            %                           6) get first derivative
            %   plotComb            =  Give vector, e.g. [1 2], to
            %                          selected which condition and in
            %                          what order.
            %   grouping            =  Seperately plot the topo's depending
            %                          on this condition
            %   trangeTopo          =  Range in ms.
            %   TargetOrResponse    =  1 - Target, 2 - Response
            %   MorletOrsFFT        =  0 - just ERP, 1 - Morlet, 2- sFFT
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -------------- PRE-SET PARAMETERS  -------------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % check for parameters.
            if ~exist('TargetOrResponse', 'var');   TargetOrResponse = 2;  end  % standard on response-locked
            if ~exist('MorletOrsFFT', 'var');       MorletOrsFFT = 0;      end  % standard on average
            
            % preset output folders
            currOutput = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName);
            if ~exist(currOutput, 'dir'); mkdir(currOutput); end
            
            currTopo = fullfile(currOutput, [plotName, 'Topo_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD) '.mat']);
            
            % preset figure folder
            averageFolder = fullfile(obj.figFolder,  'groupAverage', ['HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)], [plotName '/']);
            if ~exist(averageFolder, 'dir'); mkdir(averageFolder); end
            
            indivFolder = fullfile(obj.figFolder,  'individualPlots', ['HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)], [plotName '/']);
            if ~exist(indivFolder, 'dir'); mkdir(indivFolder); end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -------------- CREATE AVERAGE TOPOPLOT  --------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for electrode selection, we first create an average topoplot.
            % This allows us to visualize and properly chooise the cluster
            % of electrodes the check.
            % pre-set filename
            
            trangeBaseline = obj.eeg.epochPlot > obj.eeg.baseline(1) & obj.eeg.epochPlot < obj.eeg.baseline(2);
            
            if TargetOrResponse(1) == 1
                trangeTopoFA = obj.eeg.responseEpoch > trangeTopo(end,1) & obj.eeg.responseEpoch < trangeTopo(end,2);
                trangeTopo   = obj.eeg.targetEpoch > trangeTopo(1,1) & obj.eeg.targetEpoch < trangeTopo(1,2);
                warning('Will need a second input for False alarms topoplot');
            elseif TargetOrResponse(1)  == 2
                trangeTopo   = obj.eeg.responseEpoch > trangeTopo(1,1) & obj.eeg.responseEpoch < trangeTopo(1,2);
                trangeTopoFA = trangeTopo;
            end
            
            if ~exist(currTopo, 'file')
                fprintf('Extract data for avarage topoplot to determine best %s channels.\n', plotName)
                
                % 1) extract all topoplot per participant
                ERPTopo   = nan(1,obj.eeg.NumberOfChannels, max(obj.numBlocks)*obj.numTrials , length(obj.ppNames));
                ERPTopoFA = nan(1,obj.eeg.NumberOfChannels, max(obj.numBlocks)*obj.numTrials , length(obj.ppNames));
                
                for indPP = 1:length(obj.ppNames)
                    clear ERP tmp*
                    fprintf(['Now processing participant ' obj.ppNames{indPP} ' to get topoplot\n'])
                    currInput  = fullfile(obj.outputFolder,'EEG data', obj.ppNames{indPP});
                    
                    % load the EEG data.
                    if obj.eeg.applyCSD
                        load(fullfile(currInput,[obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.stim.timing '.mat']), 'csdERP');
                        ERP = csdERP;
                    else
                        load(fullfile(currInput,[obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.stim.timing '.mat']), 'ERP');
                    end
                    
                    % Baseline-correct the data for the target and
                    % response. Again not for thERPe ERPwhole
                    ERP = ERP - repmat(nanmean(ERP(:, trangeBaseline, :),2), 1, size(ERP,2), 1);
                    
                    if MorletOrsFFT ~= 0
                        clear down*
                        parfor indChan = 1:obj.eeg.NumberOfChannels-1
                            downERD(indChan,:,:,:) = obj.shortfft(ERP(indChan,:,:), freqRange);
                        end
                        [downERD(obj.eeg.NumberOfChannels,:,:,:), timeSeries] = obj.shortfft(ERP(obj.eeg.NumberOfChannels,:,:), freqRange);
                        if size(downERD,4) > 1
                            downERD = squeeze(nanmean(downERD,2));
                        end
                        
                        ERP = nan(size(downERD,1), size(obj.eeg.epochPlot,2), size(downERD,3));
                        
                        for indChan = 1:size(downERD,1)
                            for indEpoch = 1:size(downERD,3)
                                ERP(indChan,timeSeries(1):timeSeries(end)+ diff(timeSeries(1:2))-1, indEpoch) = interp(downERD(indChan,:, indEpoch), diff(timeSeries(1:2)));
                            end
                        end
                        if BaselineCorrect
                            ERP = ERP - repmat(nanmean(ERP(:, trangeBaseline, :),2), 1, size(ERP,2), 1);
                        end
                    end
                    
                    % calculated ERP topography
                    if TargetOrResponse == 1
                        [~, tmpPlot,~, tmpFAPlot] = sortERPs(obj, ERP, indPP, 0, 1);
                    elseif TargetOrResponse == 2
                        [~, ~, tmpPlot, tmpFAPlot]= sortERPs(obj, ERP, indPP, 0, 1);
                    end
                    
                    ERPTopo(1,:,:,indPP) = squeeze(nanmean(tmpPlot(:,trangeTopo,:,1),2));
                    
                    % calculated FA ERP topography
                    tmpFAPlot = nanmean(tmpFAPlot,4);
                    ERPTopoFA(1,:,:,indPP) = squeeze(nanmean(tmpFAPlot(:,trangeTopoFA,:),2));
                end
                
                save(currTopo, 'ERPTopo', 'ERPTopoFA');
                
                % 2) plot average
                obj.plotAverageTopo(ERPTopo, {plotName}, [], grouping)
                
                % temporary dock this figure as to select electrodes
                currFig = gcf;
                set(currFig,'WindowStyle','docked');
                
                % 3) hand pick a cluster of electrodes to convine the
                % electrode selection later on.
                
                if ~exist('Elec', 'var') || isempty(Elec)
                    tmpChan = input(sprintf('Please selected %s channels from topoplot', plotName), 's');
                    tmpChan = split(tmpChan);
                    for indChan = 1:length(tmpChan)
                        Elec(indChan) = find(strcmp(tmpChan{indChan}, obj.eeg.ChannelsName));
                    end
                elseif ischar(Elec)
                    tmpChan = split(Elec); clear Elec;
                    for indChan = 1:length(tmpChan)
                        Elec(indChan) = find(strcmp(tmpChan{indChan}, obj.eeg.ChannelsName));
                    end
                elseif iscell(Elec)
                    for indPP = 1:length(Elec)
                        tmpChan = split(Elec{indPP});
                        for indChan = 1:length(tmpChan)
                            indElec(indPP, indChan) = find(strcmp(tmpChan{indChan}, obj.eeg.ChannelsName));
                        end
                    end
                    clear tmp*
                    
                    % 3) plot average to check location of electrodes.
                    obj.plotAverageTopo(ERPTopo, {plotName}, unique(indElec)', grouping)
                else
                    obj.plotAverageTopo(ERPTopo, {plotName}, Elec, grouping)
                end
                
                save(currTopo, 'ERPTopo', 'ERPTopoFA','Elec');
                
            else
                load(currTopo,  'ERPTopo', 'ERPTopoFA');
                
                if ischar(Elec)
                    %load(currTopo,  'Elec');
                    tmpChan = split(Elec); clear Elec;
                    for indChan = 1:length(tmpChan)
                        Elec(indChan) = find(strcmp(tmpChan{indChan}, obj.eeg.ChannelsName));
                    end
                    tmpElecForPlot = Elec;
                elseif iscell(Elec)
                    for indPP = 1:length(Elec)
                        tmpChan = split(Elec{indPP});
                        for indChan = 1:length(tmpChan)
                            indElec(indPP, indChan) = find(strcmp(tmpChan{indChan}, obj.eeg.ChannelsName));
                        end
                    end
                    tmpElecForPlot = unique(indElec)';
                else
                    tmpElecForPlot = Elec;
                end
                
                obj.plotAverageTopo(ERPTopo, {plotName}, tmpElecForPlot, grouping)
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% ------------ CREATE AVERAGE TIMESERIES  --------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for electrode selection, we first create an average topoplot.
            % This allows us to visualize and properly chooise the cluster
            % of electrodes the check.
            if plotThis ~= 0
                clear tmp*
                if methodUsed == 1
                    fprintf('Processing raw data to extract %s using a pre-set number of channels\n', plotName)
                    currOutput   = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName,...
                        [plotName 'ChoiceElec_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)  '.mat']);
                    
                elseif methodUsed == 2
                    fprintf('Processing  raw data to extract %s  by selected the electrodes with highest signal-to-noise at response\n', plotName)
                    currOutput = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName,...
                        [plotName 'BestElec_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD) '.mat']);
                    indElec = nan(length(obj.ppNames), length(Elec));
                    
                elseif methodUsed == 3
                    fprintf('Processing  raw data to extract %s by lucalizing the cluster of best electrodes\n', plotName)
                    currOutput   = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName,...
                        [plotName 'Luc_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD) '.mat']);
                    
                end
                
                % first get number of responses (for now just first RT)
                numResp = zeros(length(obj.ppNames),1); numFA = zeros(length(obj.ppNames),1);
                for indPP = 1:length(obj.ppNames)
                    numResp(indPP) =  size(obj.behaviour{indPP}.binRT,2);
                    numFA(indPP)   =  size(obj.behaviour{indPP}.FalseAlarm,2);
                end
                
                if ~exist(currOutput, 'file')
                    ERPWhole    = nan(1, length(obj.eeg.epochPlot),     max(obj.numTrials*obj.numBlocks), length(obj.ppNames));
                    ERPTarget   = nan(1, length(obj.eeg.targetEpoch),   max(obj.numTrials*obj.numBlocks), length(obj.ppNames));
                    ERPResponse = nan(1, length(obj.eeg.responseEpoch), max(obj.numTrials*obj.numBlocks),      max(numResp), length(obj.ppNames));
                    ERPFA       = nan(1, length(obj.eeg.responseEpoch)+512, max(obj.numTrials*obj.numBlocks),  max(numFA), length(obj.ppNames));
                    
                    % We add some quality control stuff including:
                    %   - number of good trials (e.g. after artifact reject including blinks and exceed threshold)
                    %   - Standard mean error (suggested by Luck as a measurement of 'clean data', basically variance around your region of interest divided by number of trials)
                    %   - Average signal/SME
                    Quality.numGood = nan(1, length(obj.ppNames));
                    Quality.SME     = nan(1, length(obj.ppNames));
                    Quality.SNR     = nan(1, length(obj.ppNames));
                    
                    for indPP = 1:length(obj.ppNames)
                        clear ERP
                        fprintf(['Now processing extracting %s for participant ' obj.ppNames{indPP} '\n'], plotName)
                        currInput  = fullfile(obj.outputFolder,'EEG data', obj.ppNames{indPP});
                        
                        
                        if exist('indElec', 'var') && size(indElec,1) == length(obj.ppNames)
                            if ~all(isnan(indElec(indPP,:)))
                                Elec = indElec(indPP,:);
                            end
                        elseif exist('indElec', 'var') && obj.DetectOrDisc
                            Elec = indElec;
                        end
                        
                        % load the EEG data.
                        if obj.eeg.applyCSD
                            load(fullfile(currInput,[obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.stim.timing '.mat']), 'csdERP');
                            ERP = csdERP;
                        else
                            load(fullfile(currInput,[obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.stim.timing '.mat']), 'ERP');
                        end
                        
                        % Baseline-correct the data for the target and
                        % response. Again not for the ERPwhole
                        if BaselineCorrect
                            ERP  = ERP - repmat(nanmean(ERP(:, trangeBaseline, :),2), 1, size(ERP,2), 1);
                        end
                        
                        if methodUsed == 1
                            for indGroup = 1:size(Elec,1)
                                
                                currERP  = ERP(Elec(indGroup,:),:,:,:);
                                if MorletOrsFFT == 1
                                    currFreq = freqRange;
                                    if ~isnan(obj.stim.freqSSVEP) & currFreq ~= obj.stim.freqSSVEP
                                        currFreq(currFreq >= obj.stim.freqSSVEP-1 & currFreq <= obj.stim.freqSSVEP+1) = [];
                                    end
                                    
                                    ERD = obj.morletWavelength(currERP, currFreq);
                                    if size(ERD, 4) > 1
                                        currERP = squeeze(nanmean(ERD,1));
                                    else
                                        currERP = ERD;
                                    end
                                    
                                elseif MorletOrsFFT == 2
                                    [downERP, timeSeries] = obj.shortfft(currERP, freqRange);
                                    downERP = squeeze(nanmean(downERP,2));
                                    currERP = nan(size(downERP,1), size(obj.eeg.epochPlot,2), size(downERP,3));
                                    
                                    for indChan = 1:size(currERP,1)
                                        for indEpoch = 1:size(currERP,3)
                                            currERP(indChan,timeSeries(1):timeSeries(end)+ diff(timeSeries(1:2))-1, indEpoch) = interp(downERP(indChan,:, indEpoch), diff(timeSeries(1:2)));
                                        end
                                    end
                                    
                                end
                                
                                % set away the whole epoch in ERP whole
                                ERPWhole(indGroup,:,1:size(currERP,3),indPP) = nanmean(currERP,1);
                                [Quality.numGood(indGroup,indPP), ERPTarget(indGroup,:,:,indPP),...
                                    ERPResponse(indGroup,:,:,1:numResp(indPP),indPP),...
                                    ERPFA(indGroup,:,:,1:numFA(indPP),indPP)] = ...
                                    sortERPs(obj, ERPWhole(indGroup,:,:,indPP), indPP, 0, 1);
                            end
                        elseif methodUsed == 2
                            %% %%%%%%%%% Select electrodes %%%%%%%%%%%%%%%%%%
                            % selected electrodes based on topoplot.
                            % remove 0 line of virtual channel 97;
                            % 1) extract time range from response locked.
                            % calculated ERP topography
                            if TargetOrResponse == 1
                                [~, tmpERP] = sortERPs(obj, ERP, indPP, 1, 1);
                            elseif TargetOrResponse == 2
                                [~, ~, tmpERP] = sortERPs(obj, ERP, indPP, 2, 1);
                            end
                            
                            % 2) extract training samples.
                            Signal = squeeze(nanmean(tmpERP(Elec,trangeTopo,:,1),2))';
                            
                            % 3) calculated SNR
                            evalSME = nanstd(Signal)./sqrt(sum(~isnan(nanmean(Signal,2))));
                            evalSNR = nanmean(Signal)./evalSME;
                            
                            [~, indBestElec] = sort(evalSNR, 'descend');
                            indElec(indPP, 1:3) = Elec(indBestElec(1:3));
                            
                            % set away the whole epoch in ERP whole
                            ERPWhole(1,:,1:size(ERP,3),indPP) = nanmean(ERP(indElec(indPP,1:3),:,:),1);
                            
                            [Quality.numGood(1,indPP), ERPTarget(1,:,:,indPP),...
                                ERPResponse(1,:,:,1:numResp(indPP),indPP),...
                                ERPFA(1,:,:,1:numFA(indPP),indPP)] = ...
                                sortERPs(obj, ERPWhole(1,:,:,indPP), indPP, 0, 1);
                            
                        elseif methodUsed == 3
                            %% %%%%%%%%% Lucasization method %%%%%%%%%%%%%%%%%%
                            % get two small epochs. One right before the
                            % response and one of the same length in the ITI.
                            % 1) extract time range from response locked.
                            % calculated ERP topography
                            if TargetOrResponse == 1
                                [~, tmpERP] = sortERPs(obj, ERP, indPP, 1, 1);
                            elseif TargetOrResponse == 2
                                [~, ~, tmpERP] = sortERPs(obj, ERP, indPP, 2, 1);
                            end
                            
                            ERP1 = reshape(tmpERP,size(tmpERP,1), size(tmpERP,2),[]);
                            
                            % 2) same sized area from ITI
                            if strcmp(obj.stim.timing, 'future')
                                keyboard
                                ERP2 = ERP(:,obj.eeg.epochPlot  >= 1 &  obj.eeg.epochPlot  < 1 + size(ERP1,2)*(1/obj.eeg.SampleRate),:);
                            else
                                tmpWindow = obj.eeg.epochPlot  >= -(min(obj.stim.lengthITI)-1) &  obj.eeg.epochPlot  < -(min(obj.stim.lengthITI)-1)+size(ERP1,2)*(1/obj.eeg.SampleRate);
                                ERP2 = ERP(:, tmpWindow,:);
                            end
                            
                            % 4) baseline-correct both with the first 100
                            % samples.
                            ERP1 = ERP1 - nanmean(ERP1(:,1:100,:),2);
                            ERP2 = ERP2 - nanmean(ERP2(:,1:100,:),2);
                            
                            % 5) extract training samples.
                            ERP1 = ERP1(Elec, obj.eeg.responseEpoch >= trangeTopo(1) & obj.eeg.responseEpoch <= trangeTopo(2), :);
                            ERP2 = ERP2(Elec, obj.eeg.responseEpoch >= trangeTopo(1) & obj.eeg.responseEpoch <= trangeTopo(2), :);
                            
                            % 6) Concatenated in a 2d array
                            ERP1 = ERP1(:,:); ERP1(:,any(isnan(ERP1),1))  = [];
                            ERP2 = ERP2(:,:); ERP2(:,any(isnan(ERP2),1))  = [];
                            
                            % 7) get the means and covariance:
                            mu1 = mean(ERP1,2);
                            mu2 = mean(ERP2,2);
                            
                            cov1 = cov(ERP1');
                            cov2 = cov(ERP2');
                            
                            Rp = (cov1+cov2)/2; % average the covariance matrix across conditions
                            % Rp = cov1;        % Lucas thought it was also worth trying just the covariance of one condition
                            
                            % get w, the direction that will separate conditions to provide maximal linear classification accuracy:
                            w = Rp\(mu1-mu2)*100;
                            % set away the whole epoch in ERP whole
                            ERPWhole(1,:,1:size(ERP,3),indPP) =reshape(w'*ERP(Elec,:),[size(ERP,2) size(ERP,3)]);
                            
                            [Quality.numGood(1,indPP), ERPTarget(1,:,:,indPP),...
                                ERPResponse(1,:,:,1:numResp(indPP),indPP),...
                                ERPFA(1,:,:,1:numFA(indPP),indPP)] = ...
                                sortERPs(obj, ERPWhole(1,:,:,indPP), indPP, 0, 1);
                        end
                        
                        %% %%%%%%%%%%%%%%% calculated SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Do you know whether it enhances signal to noise ratio
                        % by doing this? e.g. you could use a metric like average
                        % CPP amplitude divided by the across-trial standard deviation
                        % of amplitude (in the same response window or in the baseline I guess)
                        % to show that the procedure is really providing a benefit.
                        for acc = 1:size(ERPTarget,1)
                            if TargetOrResponse == 1
                                Signal = nanmean(squeeze(ERPTarget(acc,trangeTopo,:,indPP)));
                            elseif TargetOrResponse == 2
                                Signal = nanmean(squeeze(ERPResponse(acc,trangeTopo,:,:,indPP)));
                            end
                            
                            Quality.SME(acc,indPP) = nanstd(Signal)./sqrt(sum(~isnan(Signal)));
                            Quality.SNR(acc,indPP) = nanmean(Signal)./Quality.SME(acc,indPP);
                        end
                    end
                    
                    if methodUsed == 2
                        save(currOutput,  'ERPWhole', 'ERPTarget', 'ERPResponse', 'ERPFA', 'Quality', 'indElec', '-v7.3');
                    else
                        save(currOutput,  'ERPWhole', 'ERPTarget', 'ERPResponse', 'ERPFA', 'Quality', '-v7.3');
                    end
                else
                    if methodUsed == 2
                        load(currOutput,  'ERPWhole', 'ERPTarget', 'ERPResponse', 'ERPFA', 'Quality', 'indElec');
                    else
                        load(currOutput,  'ERPWhole', 'ERPTarget', 'ERPResponse', 'ERPFA', 'Quality');
                    end
                end
                
                %% actually doing the plotting
                % create low-pass filter for plotting purposes only!
                Fs    = obj.eeg.SampleRate;
                fco   = 6;
                [b,a] = butter(2,fco*2/Fs);
                
                plotTarget   = nan(size(ERPTarget));
                plotResponse = nan(size(ERPResponse));
                plotFA       = nan(size(ERPFA));
                
                for indPP = 1:length(obj.ppNames)
                    fprintf('.')
                    % filter ERPs as to remove SSVEP for example, e.g.
                    % smoothing the signal.
                    clear tmp*
                    for indEpoch = 1:size(obj.behaviour{indPP}.trialMatrix,1)
                        for indChanComb = 1:size(ERPTarget,1)
                            isNotNaN = ~isnan(ERPTarget(indChanComb,:,indEpoch,indPP));
                            plotTarget(indChanComb,isNotNaN,indEpoch,indPP) = (filtfilt(b,a,ERPTarget(indChanComb,isNotNaN,indEpoch,indPP))); % for plotting purposes only
                            
                            for indResponse = 1:size(ERPResponse,4)
                                isNotNaN = ~isnan(ERPResponse(indChanComb,:,indEpoch,indResponse,indPP));
                                plotResponse(indChanComb,isNotNaN,indEpoch,indResponse,indPP) = (filtfilt(b,a, ERPResponse(indChanComb,isNotNaN,indEpoch,indResponse,indPP)));
                            end
                            
                            % remove baseline
                            if BaselineCorrect == 1
                                tmpBaseline = nanmean(plotTarget(indChanComb,obj.eeg.targetEpoch >= obj.eeg.baseline(1) & obj.eeg.targetEpoch <= obj.eeg.baseline(2),indEpoch,indPP),2);
                            elseif BaselineCorrect == 2
                                
                                % forLaterBaseline(indChanComb,indEpoch) = nanmean(plotResponse(indChanComb,obj.eeg.responseEpoch >= obj.eeg.baseline(1) & obj.eeg.responseEpoch <= obj.eeg.baseline(2),indEpoch,1, indPP),2);
                                tmpBaseline = 0;

                            else
                                tmpBaseline = 0;
                            end
                            
                            plotTarget(indChanComb,:,indEpoch,indPP) = plotTarget(indChanComb,:,indEpoch, indPP) - tmpBaseline;
                            for indResponse = 1:size(ERPResponse,4)
                                plotResponse(indChanComb,:,indEpoch,indResponse, indPP) = plotResponse(indChanComb,:,indEpoch,indResponse, indPP)  - tmpBaseline;
                            end
                            
                            for indFA = 1:size(obj.behaviour{indPP}.FalseAlarm,2)
                                if obj.behaviour{indPP}.FalseAlarm(indEpoch,indFA)
                                    try
                                    isNotNaN = ~isnan(ERPFA(indChanComb,:,indEpoch,indFA,indPP));
                                    plotFA(indChanComb,isNotNaN,indEpoch,indFA,indPP) = filtfilt(b,a, ERPFA(indChanComb,isNotNaN,indEpoch,indFA,indPP));
                                    catch
                                        keyboard
                                    end
                                end
                            end
                        end
                    end
                     
                    if  obj.DetectOrDisc && size(plotTarget, 1) > 1
                        
%                         plotResponse (1,:,obj.behaviour{indPP}.Reaction == 1 & obj.behaviour{indPP}.Correct, :,indPP) = nan;
%                         plotResponse (2,:,obj.behaviour{indPP}.Reaction == 2 & obj.behaviour{indPP}.Correct, :,indPP) = nan;
%                         plotResponse (:,:,obj.behaviour{indPP}.Reaction == 0, :,indPP) = nan;
%                         
                        plotTarget(1,:,:,indPP)     = nanmean(plotTarget(:,:,:,indPP),1);
                        plotResponse(1,:,:,:,indPP) = nanmean(plotResponse(:,:,:,:,indPP),1);
                    end
                    
                    % save current plot
                    if plotThis == 1 || plotThis == 0
                        fprintf('Plotting individual plots.\n')
                        
                        try
                            figHandle = obj.plotIndividualTimeSeries(indPP, plotTarget(1,:,:,indPP), plotResponse(1,:,:,:,indPP), plotName, plotComb, grouping);
                            if exist('indElec', 'var')
                                Elec = indElec(indPP,:);
                            end
                            for indHandle = 1:length(figHandle)
                                figure(figHandle{indHandle})
                                obj.plotIndividualTopo(ERPTopo(:,:,:, indPP), Elec)
                                plotSave(gca,  [plotName '_' obj.ppNames{indPP} num2str(indHandle) '.png'], indivFolder, obj.figLayOut.saveDim);
                            end
                        end
                        figInfo = [];
                    end
                end
                
                % plot grand-averages (1 is always plot individual)!
                if  plotThis ~= 1
                    fig = obj.plotAverageTimeseries(plotTarget(1,:,:,:), plotResponse(1,:,:,:,:), plotFA, BaselineCorrect, plotName, plotComb, grouping, plotThis);  % plot these average.
                    for indHandle = 1:length(fig.Handle)
                        
                        figure(fig.Handle{indHandle})
                        plotSave(gca,  [plotName num2str(indHandle)  num2str(plotThis) '.png'], averageFolder, obj.figLayOut.saveDim);
                        
                        legend(obj.figLayOut.legends, 'Location', 'best');
                        plotSave(gca,  [plotName num2str(indHandle) num2str(plotThis) 'Legend.png'], averageFolder, obj.figLayOut.saveDim);
                        
                        if plotThis >= 4
                            figure(fig.ResponseHandle{indHandle})
                            plotSave(gca,  [plotName num2str(indHandle) 'Response.png'], averageFolder, [obj.figLayOut.saveDim(1) obj.figLayOut.saveDim(2).*0.35] );
                        end
                        
                        if obj.stim.FA && plotThis == 2
                            figure(fig.FAHandle{indHandle})
                            plotSave(gca,  [plotName num2str(indHandle) 'FA.png'], averageFolder, obj.figLayOut.saveDim);
                        end
                        
                        if plotThis == 2
                            figure(fig.MissesHandle{indHandle})
                            plotSave(gca,  [plotName num2str(indHandle) 'Misses.png'], averageFolder, obj.figLayOut.saveDim);
                        end
                    end
                else
                    fig = [];
                end
            else
                fig = [];
            end
            
            fprintf('.\n')
            fprintf('Done, go get yourself a coffee :) \n')
        end
        
        
        function [donetbl, waveform, plotData] = getWaveformParameters(obj, methodUsed, plotName, qps, BaselineCorrect, plotThis, plotComb,...
                grouping, trangeTopo, TargetOrResponse, negOrPos, PeakMeanOrMax, amplitudeAxis)
            % check for parameters.
            if ~exist('TargetOrResponse', 'var');   TargetOrResponse = 2;  end  % standard on response-locked
            if ~exist('plotComb', 'var'); plotComb = []; end
            if ~exist('qps', 'var'); qps   = [0.33 0.66]; end
            if ~exist('amplitudeAxis', 'var'); amplitudeAxis   = [ ]; end
            
            % preset folders
            saveFolder = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName, ['RTBINS\']);
            if ~exist(saveFolder, 'dir'); mkdir(saveFolder); end
            
            % preset figure save folder
            hereFigFolder = fullfile(obj.figFolder,  'groupAverage', ['HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)], plotName, ['RTBINS' num2str(length(qps)) '/']);
            if ~exist(hereFigFolder, 'dir'); mkdir(hereFigFolder); end
            
            % for Individual plots
            indivFolder = fullfile(obj.figFolder,  'individualPlots', ['HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)], plotName, ['RTBINS' num2str(length(qps)) '/']);
            if ~exist(indivFolder, 'dir'); mkdir(indivFolder); end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -------------- Set-up conditions   -------------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% set up model parameters
            % as there are different parameters that you might wanna add,
            % like time dependency. For now it just simple time on task and
            % time within the experiment.
            
            [within, between, modelFun, allTrialMatrix] = getConditions(obj, plotComb);
            
            % as we want to include Reaction time effect we here are
            % setting up the quantiles and RT as an extra within parameter.
            
            if ~isempty(qps)
                within.Name = {within.Name{:} 'RTs' };
                within.Par  = [within.Par 1];
                within.Var  = {within.Var{:} 1:size(qps)+1};
                
                tmpCondition = [];
                
                for indCond = 1:within.numCond
                    tmpCondition = [tmpCondition; repmat(within.Condition(indCond,:),length(qps)+1,1), (1:length(qps)+1)'];
                end
                
                within.Condition = tmpCondition;
                within.Table     = array2table(within.Condition, 'VariableNames', within.Name);
                
                % model set-up for ranova and plotting.
                switch length(within.Var)
                    case 1
                        modelFun = sprintf('%s',...
                            within.Name{1});  % main effects
                    case 2
                        modelFun = sprintf('%s + %s + %s * %s',...
                            within.Name{1}, within.Name{2},... % main effects
                            within.Name{1}, within.Name{2});
                    case 3
                        modelFun = sprintf('%s + %s + %s + %s * %s + %s * %s + %s * %s + %s * %s * %s',...
                            within.Name{1}, within.Name{2}, within.Name{3},... % main effects
                            within.Name{1}, within.Name{2},  within.Name{1}, within.Name{3},... % two-way interactions
                            within.Name{2}, within.Name{3},...
                            within.Name{1}, within.Name{2}, within.Name{3}); % three-way interactions
                    case 4
                        modelFun = sprintf('%s + %s + %s + %s + %s * %s + %s * %s + %s * %s + %s * %s + %s * %s + %s * %s + %s * %s * %s  + %s * %s * %s  + %s * %s * %s  + %s * %s * %s * %s',...
                            within.Name{1}, within.Name{2}, within.Name{3}, within.Name{4},... % main effects
                            within.Name{1}, within.Name{2},  within.Name{1}, within.Name{3},  within.Name{1}, within.Name{4},... % two-way interactions
                            within.Name{2}, within.Name{3}, within.Name{2}, within.Name{4},...
                            within.Name{3}, within.Name{4},...
                            within.Name{1}, within.Name{2}, within.Name{3}, within.Name{1}, within.Name{2}, within.Name{4},...% three-way interactions
                            within.Name{2}, within.Name{3}, within.Name{4},...
                            within.Name{1}, within.Name{2}, within.Name{3},  within.Name{4}); % four-way interactions
                    otherwise
                        keyboard
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% ---- create signal and extract possible DM parameters ------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Re-load the waveforms here and sort on z-score transform RT
            % per participant, signals will than be sorted, binned based on
            % conditions after which the DM accumuatlion paramters are
            % extracted. This includes 1) Onset, 2) Slope and 3) peak
            % average.
            
            if ~exist(fullfile(saveFolder, ['Parameters ' num2str(plotComb) num2str(TargetOrResponse) '_' num2str(BaselineCorrect) '.mat']), 'file')
                clear tmp*
                if methodUsed == 1
                    currOutput	= fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName, [plotName 'ChoiceElec_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)  '.mat']);
                elseif methodUsed == 2
                    currOutput  = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName, [plotName 'BestElec_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD) '.mat']);
                elseif methodUsed == 3
                    currOutput  = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName, [plotName 'Luc_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD) '.mat']);
                end
                
                load(currOutput, 'ERPTarget', 'ERPResponse');
                clc
                fprintf('Automatic detecting the onset/peak/slopes: 0 %%\n')
                
                % Prep signals for further parameters extraction.
                for indPP = 1:length(obj.ppNames)
                    
                    % 2) baseline-correct the trials if indicated
                    % remove baseline 
                    if size(ERPResponse,1) > 1
%                         if size(ERPResponse,1) > 1
%                             ERPResponse (1,:,obj.behaviour{indPP}.Reaction == 1 & obj.behaviour{indPP}.Correct, :,indPP) = nan;
%                             ERPResponse (2,:,obj.behaviour{indPP}.Reaction == 2 & obj.behaviour{indPP}.Correct, :,indPP) = nan;
%                             ERPResponse (:,:,obj.behaviour{indPP}.Reaction == 0, :,indPP) = nan;
%                         end
                        ERPResponse(1,:,:,:,indPP) = nanmean(ERPResponse(:,:,:,:,indPP),1);
                        ERPTarget(1,:,:,indPP)   = nanmean(ERPTarget(:,:,:,indPP),1);
                    end
                    
                    if BaselineCorrect == 1
                        tmpBaseline = nanmean(ERPTarget(1, obj.eeg.targetEpoch >= obj.eeg.baseline(1) & obj.eeg.targetEpoch <= obj.eeg.baseline(2), :, indPP),2);
                    elseif BaselineCorrect == 2
                        
%                         tmpBaseline = nanmean(ERPResponse(:, obj.eeg.responseEpoch >= obj.eeg.baseline(1) & obj.eeg.responseEpoch <= obj.eeg.baseline(2), :, 1, indPP),2);
                        tmpBaseline = 0;
                    else
                        tmpBaseline = 0;
                    end
                   
                                            
                    ERPTarget(:,:,:,indPP)     = ERPTarget(:,:,:, indPP)    - repmat(tmpBaseline,1,size(ERPTarget,2),1);
                    ERPResponse(:,:,:,:,indPP) = ERPResponse(:,:,:,:,indPP) - repmat(tmpBaseline,1,size(ERPResponse,2),1, size(ERPResponse,4));
                end
                
                
                % 1) get max amplitude according at response
                if TargetOrResponse == 1
                    PeakArea      = obj.eeg.targetEpoch > trangeTopo(1,1) & obj.eeg.targetEpoch < trangeTopo(1,2);
                    slopeArea     = obj.eeg.targetEpoch > trangeTopo(2,1) & obj.eeg.targetEpoch < trangeTopo(2,2);
                    timeSeries{1} = obj.eeg.targetEpoch(PeakArea);
                    timeSeries{2} = obj.eeg.targetEpoch(slopeArea);
                    ERPTarget(1,:,:,:) = nanmean(ERPTarget, 1);
                    Signal = ERPTarget;
                elseif TargetOrResponse == 2
                    PeakArea  = obj.eeg.responseEpoch > trangeTopo(1,1) & obj.eeg.responseEpoch < trangeTopo(1,2);
                    slopeArea = obj.eeg.responseEpoch > trangeTopo(2,1) & obj.eeg.responseEpoch < trangeTopo(2,2);
                    
                    timeSeries{1} = obj.eeg.responseEpoch(PeakArea);
                    timeSeries{2} = obj.eeg.responseEpoch(slopeArea);
                    Signal(1,:,:,:)  = squeeze(ERPResponse(:,:,:,1,:));
                end
                if size(Signal, 1) > 1
                    Signal(2,:,:,:) = [];
                end
                
                % TODO create proper table for LLM
                tbl.ppNames      = [];
                tbl.RT           = [];
                tbl.zScoreRT     = [];
                tbl.Slopes       = [];
                tbl.Peak         = [];
                tbl.PeakLatency	 = [];
                
                
                for indCond = 1:length(plotComb)
                    eval(sprintf('tbl.%s =  [];',  within.Name{indCond}));
                end
                
                for indPP = 1:length(obj.ppNames)
                    clc
                    fprintf('Automatic detecting the onset/peak/slopes: %i %% \n', round((indPP/length(obj.ppNames))*100))
                    
                    if nanmean(obj.behaviour{indPP}.artifacts & obj.behaviour{indPP}.blinks) >= 0.3
                        % 4) go through given conditions and sort, bin and
                        % average signals
                        goodTrials      = obj.behaviour{indPP}.artifacts' & obj.behaviour{indPP}.blinks';
                        VartrialMatrix  = obj.behaviour{indPP}.trialMatrix(:, plotComb);
                        posCombination  = unique(VartrialMatrix(~isnan(obj.behaviour{indPP}.RT(:,1)),:), 'rows');
                        posCombination(any(isnan(posCombination),2),:) = [];
                        
                        
                        currtblRT           = nan(size(goodTrials));
                        currtblzScoreRT     = nan(size(goodTrials));
                        currtblSlopes       = nan(size(goodTrials));
                        currtblPeak         = nan(size(goodTrials));
                        
                        currtblPeakLatency  = nan(size(goodTrials));
                        
                        acc = 1;
                        for indCond = 1:size(posCombination,1)
                            clear tmp*
                            tmpCond = find(sum(VartrialMatrix == posCombination(indCond,:),2) == size(posCombination, 2) & goodTrials);
                            
                            if ~isempty(tmpCond)
                                
                                % get the reaction times to create the
                                % timing there
                                tmpRT = obj.behaviour{indPP}.RT(tmpCond);
                                
                                % get equal size bins
                                if isempty(qps)
                                    tmpRTqps = discretize(tmpRT, [0 max(tmpRT)], 'IncludedEdge','right');
                                    tmpRTqps(isnan(tmpRTqps)) =  length(qps)+2; % add the misses
                                else
                                    tmpRTqps = discretize(tmpRT, [0 quantile(tmpRT,qps)-0.5/60 max(tmpRT)], 'IncludedEdge','right');
                                    tmpRTqps(isnan(tmpRTqps))   = length(qps)+2; % add the misses
                                end
                                
                                currtblRT(tmpCond)       = tmpRT;
                                currtblzScoreRT(tmpCond) = nanzscore(tmpRT);
                                % tmptblRTQuar(tmpCond)  = tmpRTqps;
                                
                                % extract the aligned ERPs for Amplitude
                                % and slope extraction
                                signalPeak  = squeeze(Signal(:,:,tmpCond,indPP)); 
                                if BaselineCorrect == 2
                                    tmpBaseline = nanmean(nanmean(ERPResponse(1, obj.eeg.responseEpoch >= obj.eeg.baseline(1) &...
                                        obj.eeg.responseEpoch <= obj.eeg.baseline(2), tmpCond, 1, indPP),2),3);
                                    signalPeak = signalPeak - tmpBaseline;
                                end
                                
                                % now extract all mean peak amplitudes
                                if PeakMeanOrMax(1) == 1
                                    tmpAmplitude = nanmean(signalPeak(PeakArea,:),1);
                                    tmpLatency   = nanmean(timeSeries{1});
                                    currtblPeakLatency(tmpCond) = tmpLatency;
                                    tmpLatency = repmat(tmpLatency, size(tmpCond));
                                elseif PeakMeanOrMax(1) == 2
                                    if negOrPos == 1
                                        [tmpAmplitude,tmpLatency] = min(signalPeak(PeakArea,:));
                                    elseif negOrPos == 2
                                        [tmpAmplitude,tmpLatency] = max(signalPeak(PeakArea,:));
                                    end
                                    tmpLatency = timeSeries{1}(tmpLatency)';
                                    currtblPeakLatency(tmpCond) = tmpLatency;
                                end
                                
                                % first also get quick slope determining if
                                % it is negative to exclude trials.
                                
                                signalSlope = squeeze(Signal(:,:,tmpCond,indPP)); signalSlope = sgolayfilt(signalSlope,1,51,[],1);
                                 

                                clear tmpSlope
                                for indEpoch = 1:length(tmpRT)
                                    if ~isnan(tmpRT(indEpoch))
                                        indSlope = obj.eeg.responseEpoch> tmpLatency(indEpoch)-(tmpRT(indEpoch)) &  obj.eeg.responseEpoch < tmpLatency(indEpoch); %
                                        x = obj.eeg.responseEpoch(indSlope);
                                        y = signalSlope(indSlope, indEpoch)';
                                        
                                        % Fit line to data using polyfit
                                        c = polyfit(x,y,1);
                                        % Evaluate fit equation using polyval
                                        tmpSlope(indEpoch,:) = c(1);
                                        
                                        % figure, plot(obj.eeg.responseEpoch, signalSlope(:, indEpoch))
                                        % hold on, plot(x,y)
                                    else
                                        tmpSlope(indEpoch,:) = NaN;
                                    end
                                end
                                
                                currtblSlopes(tmpCond) = tmpSlope;
                                currtblPeak(tmpCond)   = tmpAmplitude';
                                
                                % when not interested in rt effects (just
                                % collapsing over for plotting) we only add
                                % the hits when it is response locked. When
                                % we look at target-locked we could
                                % eventually add misses as well. 
                                
                                if isempty(qps) 
                                    numBins = 1;
                                else
                                    numBins = length(qps)+1;
                                end
                                
                                for indQps = 1:numBins
                                    
                                    currQps  = tmpRTqps == indQps;
                                    waveform.RT(indPP,acc) = nanmedian(tmpRT(currQps));
                                    
                                    %%%%%%%%%%%%%%%% extract peak %%%%%%%%%%%%
                                    if sum(currQps) == 0 || isnan(waveform.RT(indPP,acc))
                                        waveform.peak(indPP, acc)        = NaN;
                                        waveform.peakLatency(indPP, acc) = NaN;
                                        
                                    elseif PeakMeanOrMax(1) == 1 && ~isempty(currQps) && ~isnan(waveform.RT(indPP,acc))
                                    
                                        waveform.peak(indPP, acc)        = nanmean(nanmean(signalPeak(PeakArea,currQps),1),2);
                                        waveform.peakLatency(indPP, acc) = nanmean(timeSeries{1});
                                        
                                    elseif PeakMeanOrMax(1) == 2
                                        % 2.2) reshape rt into vector and sort and allRT to get line.
                                        tmpOnset = nanmean(signalPeak(PeakArea, currQps),2);
                                        
                                        if negOrPos == 1
                                            [peak,indMax] = min(tmpOnset);
                                        elseif negOrPos == 2
                                            [peak,indMax] = max(tmpOnset);
                                        end
                                        
                                        waveform.peak(indPP, acc) = peak;
                                        
                                        if ~isnan(peak)
                                            waveform.peakLatency(indPP, acc) = timeSeries{1}(indMax);
                                        else
                                            waveform.peakLatency(indPP, acc) = NaN;
                                        end
                                    end
                                    
                                    %%%%%%%%%%%%%%%% extract slope %%%%%%%%%%%%%%%%
                                    % first calculated onset
                                    currSignal = signalSlope(slopeArea, currQps)';
                                    
                                    % create the timeseries and calculated
                                    % slopes
                                    if ~isempty(currSignal) | ~all(isnan(currSignal))
                                        waveform.slopes(indPP, acc)     = nanmean(tmpSlope(currQps),1);
                                        waveform.numTrials(indPP, acc)	= sum(currQps);
                                    else
                                        waveform.slopes(indPP, acc)   = nan;
                                        waveform.interval(indPP, acc) = nan;
                                    end
                                    
                                    waveform.RT(indPP, acc) = nanmedian(tmpRT(currQps));
                                    
                                    acc = acc + 1;
                                end
                            end
                        end
                        
                        tbl.ppNames      = [tbl.ppNames;     repmat(indPP, size(currtblRT))];
                        tbl.RT           = [tbl.RT;          nanzscore(currtblRT)];
                        tbl.zScoreRT     = [tbl.zScoreRT;    currtblzScoreRT];
                        tbl.Slopes       = [tbl.Slopes;      currtblSlopes];
                        tbl.Peak         = [tbl.Peak;        currtblPeak];
                        tbl.PeakLatency	 = [tbl.PeakLatency; currtblPeakLatency];
                        
                    end
                    
                    % we want to formally test if there are effects of Evidence stength and
                    % context while accounting for participants differences as well as
                    % Reaction time effects
                    for indCond = 1:length(plotComb)
                        eval(sprintf('tbl.%s =  [tbl.%s; allTrialMatrix(1:length(goodTrials),indCond,indPP)];',  within.Name{indCond}, within.Name{indCond}));
                    end
                end
                
                save(fullfile(saveFolder, ['RTBins_' num2str(length(qps))  'Parameters ' num2str(plotComb) num2str(TargetOrResponse(1)) '_' num2str(BaselineCorrect) '.mat']), 'Signal', 'timeSeries', 'slopeArea', 'tbl', 'waveform')
            else
                load(fullfile(saveFolder, ['RTBins_' num2str(length(qps))  'Parameters ' num2str(plotComb) num2str(TargetOrResponse(1)) '_' num2str(BaselineCorrect) '.mat']), 'Signal', 'timeSeries', 'slopeArea', 'tbl', 'waveform')
            end
           
            % we want to formally test if there are effects of Evidence stength and
            % context while accounting for participants differences as well as
            % Reaction time effects
            tbl.ppNames = nominal(tbl.ppNames);
            
            donetbl = struct2table(tbl);
            donetbl(isnan(donetbl.RT),:) = [];
            
            fileID = fopen(fullfile(obj.logFolder, ['statisticsWaveform' plotName '.log']),'w');
            fprintf(fileID, '-------------- Statistical analysis --------------\n');
            fprintf(fileID, 'comparing waveform parameters between the different conditions\n');
            fprintf(fileID, 'An repeated-measurement ANOVA was applied to see if \n');
            fprintf(fileID, 'the reaction times for the different conditions come from the same distribution.\n');
            
            strCond = 'Conditions: ';
            for indCond = 1:size(within.Name,2)
                strCond = strcat(strCond, [within.Name{indCond} ', ']);
                eval(sprintf('within.Table.%s = categorical(within.Table.%s);', within.Name{indCond}, within.Name{indCond}));
            end
            
            if sum(~within.Par)
                for indCond = 1:size(between.Name,2)
                    strCond = strcat(strCond, ['Groups:' between.Name{indCond} ', ']);
                end
            end
            
            %% run the ranovas on the experiment.
            % ------------- Peak ---------------------------------
            % write away the descriptives statistics
            % remove inter-subject variabiltiy
            
            if plotThis ~= 0
                noEnoughPP = (sum(isnan(waveform.RT)) >= length(obj.ppNames)*0.5);
                
                within.Table(noEnoughPP,:) = [];
                waveform.RT(:,noEnoughPP)  = [];
                waveform.peak(:,noEnoughPP)  = [];
                waveform.slopes(:,noEnoughPP)  = [];

                % remove within-participant variability.
                indMean   = repmat(nanmean(waveform.peak,2), 1, size(waveform.peak,2));
                grandMean = nanmean( waveform.peak(:));
                waveform.peak = waveform.peak - indMean + grandMean;
                
                
                indMean   = repmat(nanmean(waveform.slopes,2), 1, size(waveform.slopes,2));
                grandMean = nanmean( waveform.slopes(:));
                waveform.slopes = waveform.slopes - indMean + grandMean;
                
                if any(strcmp(within.Name, 'RTs'))
                    posCombination = unique(within.Condition(:,1:end-1), 'row');
                else
                    posCombination = within.Condition;
                end
                
                [rm.Peak, results.Peak] = preformRANOVA(obj, fileID, 'Peak Amplitude', strCond,  waveform.peak, between, within, modelFun);
                [rm.Slopes, results.Slopes] = preformRANOVA(obj, fileID, 'Slopes Amplitude', strCond,  waveform.slopes, between, within, modelFun);

                rm.RT = preformRANOVA(obj, fileID, 'RT Amplitude', strCond,  waveform.RT, between, within, modelFun);
                
                if sum(~within.Par)
                    plotData   = margmean(rm.Peak, {within.Name between.Name{1}});
                    plotSlopes   = margmean(rm.Slopes, {within.Name between.Name{1}});

                    plotRT     = margmean(rm.RT, {within.Name between.Name{1}});
                else
                    plotData   = margmean(rm.Peak, within.Name);
                    plotSlopes   = margmean(rm.Slopes, within.Name);

                    plotRT     = margmean(rm.RT, within.Name);
                end
                
                fclose(fileID)
                
                % plot results
                fig1 =  figure('units','normalized','outerposition',[0 0 1 1]);
                hold on
                if isempty(qps)
                    clear averageSignal stdSignal averageRT
                    % NOTE THIS HAS BEEN LAZY CODING! AND REALLY NEEDS TO
                    % BE CHANGED :'(
                    if size(within.Condition,1) == 4
                        errorbar([1:2]-0.1, plotData.Mean([1 2]), plotData.StdErr([1 2]), 'LineWidth', 1, 'LineStyle', '-.',...
                            'Color', obj.figLayOut.colours(2,:), 'Marker', '+')
                        errorbar([1:2]+0.1, plotData.Mean([3 4]), plotData.StdErr([3 4]), 'LineWidth', 1, 'LineStyle', '-.',...
                            'Color', obj.figLayOut.colours(4,:), 'Marker', '+')
                        
                        xlim([1-0.2 2+0.2])
                        xticks([1 2])
                        %}
                    elseif size(within.Condition,1) == 8
                        errorbar([1:4]-0.1, plotData.Mean([1:4]), plotData.StdErr([1:4]), 'LineWidth', 1, 'LineStyle', '-.',...
                            'Color', obj.figLayOut.colours(1,:), 'Marker', '+')
                        errorbar([1:4]+0.1, plotData.Mean([5:8]), plotData.StdErr([5:8]), 'LineWidth', 1, 'LineStyle', '-.',...
                            'Color', obj.figLayOut.colours(2,:), 'Marker', '+')
                        xlim([1-0.2 4+0.2])
                        
                        xticks([1:4])
                    elseif size(within.Condition,1) == 16
                        keyboard
                         errorbar([1:4]-0.2, plotData.Mean([1:4]), plotData.StdErr([1:4]), 'LineWidth', 1, 'LineStyle', '-.',...
                            'Color', obj.figLayOut.colours(1,:), 'Marker', '+')
                        errorbar([1:4]-0.1, plotData.Mean([5:8]), plotData.StdErr([5:8]), 'LineWidth', 1, 'LineStyle', '-.',...
                            'Color', obj.figLayOut.colours(2,:), 'Marker', '+')
                           errorbar([1:4]+0.1, plotData.Mean([9:12]), plotData.StdErr([9:12]), 'LineWidth', 1, 'LineStyle', '-.',...
                            'Color', obj.figLayOut.colours(3,:), 'Marker', '+')
                       errorbar([1:4]+0.2, plotData.Mean([13:16]), plotData.StdErr([13:16]), 'LineWidth', 1, 'LineStyle', '-.',...
                            'Color', obj.figLayOut.colours(4,:), 'Marker', '+')
                    
                        xlim([1-0.2 4+0.2])
                        
                        xticks([1:4])
                    end
                    xticklabels({obj.figLayOut.legNames{plotComb(2)}{:}});
                    xlabel(' ');
                    
                    if ~isempty(amplitudeAxis)
                        ylim([amplitudeAxis]);
                    end
                    
                    set(gca,'FontSize', obj.figLayOut.letterSize);
                    set(gca,'FontName', obj.figLayOut.letterType);
                    
                    plotSave(gca, ['amplitude' plotName '.png'], hereFigFolder, obj.figLayOut.saveDim);
                else
                    figure(fig1)

                    posCombination = categorical(posCombination);
                    for indCond = 1: within.numCond
                        currPlot   = all( table2array ( plotData(:,1:size(posCombination,2))) == categorical(posCombination(indCond,:)), 2);
                        averageRT(indCond, :)     = plotRT.Mean(currPlot);
                        averageSignal(indCond, :) = plotData.Mean(currPlot);
                        stdSignal(indCond, :)     = plotData.StdErr(currPlot);
                        
                        legendThis(indCond) = errorbar(averageRT(indCond,:), averageSignal(indCond,:), stdSignal(indCond,:));
                        legendThis(indCond).Color = obj.figLayOut.colours(indCond,:);
                        legendThis(indCond).LineWidth = 1;
                        legendThis(indCond).LineStyle = obj.figLayOut.lineType{indCond};
                        legendThis(indCond).Marker = '.';
                    end
                    
                    if ~isempty(amplitudeAxis)
                        ylim([amplitudeAxis]);
                    end
                    plotSave(gca, ['amplitude' plotName '.png'], hereFigFolder, obj.figLayOut.saveDim);
                    
                    figure; hold on

                    posCombination = categorical(posCombination);
                    for indCond = 1: within.numCond
                        currPlot   = all( table2array ( plotData(:,1:size(posCombination,2))) == categorical(posCombination(indCond,:)), 2);
                        averageRT(indCond, :)     = plotRT.Mean(currPlot);
                        averageSignal(indCond, :) = plotSlopes.Mean(currPlot);
                        stdSignal(indCond, :)     = plotSlopes.StdErr(currPlot);
                        
                        legendThis(indCond) = errorbar(averageRT(indCond,:), averageSignal(indCond,:), stdSignal(indCond,:));
                        legendThis(indCond).Color = obj.figLayOut.colours(indCond,:);
                        legendThis(indCond).LineWidth = 1;
                        legendThis(indCond).LineStyle = obj.figLayOut.lineType{indCond};
                        legendThis(indCond).Marker = '.';
                    end
                    
                    if ~isempty(amplitudeAxis)
                        ylim([amplitudeAxis]);
                    end
                    plotSave(gca, ['slopes' plotName '.png'], hereFigFolder, obj.figLayOut.saveDim);
                end
            end
        end
        
        
        % extract and plot ITI
        function fig = plotITI(obj, methodUsed, plotName, Elec, BaselineCorrect, plotComb, grouping, thresArea)
            %% fig = plotITI(obj, methodUsed, plotName, Elec, BaselineCorrect, plotComb, grouping, thresArea)
            % TODO
            % Parameters:
            %   MethodUsed          =  gives the possibility to just use
            %                          1) the pre-set electrodes
            %                          2) get the 3 best electrodes with highest SNR
            %                          3) Lucalize a cluster of electrodes
            %   plotName            =  given name to the plots (e.g. CPP)
            %   Elec                =  selected electrodes (either number
            %                          or names), when empty the electrodes
            %                          will be asked to be selected during
            %                          the topoplot phase.
            %   BaselineCorrect     =  0 - no, 1 - yes
            %   plotThis            =  0 - nothing, 1 - individual or 2 - average
            %   plotComb            =  Give vector, e.g. [1 2], to
            %                          seleceted which condition and in
            %                          what order.
            %   grouping            =  Seperately plot the topo's depending
            %                          on this condition
            %   trangeTopo          =  Range in ms.
            %   TargetOrResponse    =  1 - Target, 2 - Response
            %   MorletOrsFFT        =  0 - just ERP, 1 - Morlet, 2- sFFT
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -------------- PRE-SET PARAMETERS  -------------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % create low-pass filter for plotting purposes only!
            Fs    = obj.eeg.SampleRate;
            fco   = 6;
            [b,a] = butter(2,fco*2/Fs);
            
            % preset figure save folder
            averageFolder = fullfile(obj.figFolder,  'groupAverage', ['HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)], [plotName '/ITI/']);
            if ~exist(averageFolder, 'dir'); mkdir(averageFolder); end
            
            indivFolder = fullfile(obj.figFolder,  'individualPlots', ['HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)], [plotName '/ITI/']);
            if ~exist(indivFolder, 'dir'); mkdir(indivFolder); end
            
            % preset folders
            clear tmp*
            if methodUsed == 1
                fprintf('Processing raw data to extract %s using a pre-set number of channels\n', plotName)
                currOutput   = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName, [plotName 'ChoiceElec_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)  '.mat']);
            elseif methodUsed == 2
                fprintf('Processing  raw data to extract %s  by selected the electrodes with highest signal-to-noise at response\n', plotName)
                currOutput = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName, [plotName 'BestElec_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD) '.mat']);
                indElec = nan(length(obj.ppNames), length(Elec));
            elseif methodUsed == 3
                fprintf('Processing  raw data to extract %s by lucalizing the cluster of best electrodes\n', plotName)
                currOutput   = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName, [plotName 'Luc_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD) '.mat']);
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% ------------ CREATE AVERAGE TIMESERIES  --------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for electrode selection, we first create an average topoplot.
            % This allows us to visualize and properly chooise the cluster
            % of electrodes the check.
            
            % pre-allocated the set-away matricis to be saved for plotting
            % (TODO save)
            load(currOutput,  'ERPWhole', 'ERPResponse');
            
            plotWhole   = nan(size(ERPWhole));
            plotResponse = nan(size(ERPResponse));
            for indPP = 1:length(obj.ppNames)
                fprintf('.')
                
                % filter ERPs as to remove SSVEP for example, e.g.
                % smoothing the signal.
                clear tmp*
                for indEpoch = 1:size(obj.behaviour{indPP}.trialMatrix,1)
                    for indChanComb = 1:size(ERPWhole,1)
                        
                        isNotNaN = ~isnan(ERPWhole(indChanComb,:,indEpoch,indPP));
                        
                        if ~all(isNotNaN == 0) && ~(obj.behaviour{indPP}.Misses(indEpoch))
                            
                            plotWhole(indChanComb,isNotNaN,indEpoch,indPP) = filtfilt(b,a,ERPWhole(indChanComb,isNotNaN,indEpoch,indPP));
                            
                            % remove baseline
                            if BaselineCorrect == 1
                                tmpBaseline     = nanmean(plotWhole(indChanComb,obj.eeg.epochPlot >= obj.eeg.baseline(1) & obj.eeg.epochPlot <= obj.eeg.baseline(2),indEpoch,indPP),2);
                            elseif BaselineCorrect == 2
                                currStart = find(isNotNaN);
                                if ~isempty(currStart)
                                    tmpBaseline = nanmean(plotWhole(indChanComb,currStart(1):currStart(100), indEpoch,indPP),2);
                                else
                                    tmpBaseline = NaN;
                                end
                            elseif BaselineCorrect == 3
                                tmpBaseline = 0;
                                plotResponse(indChanComb,:, indEpoch,1, indPP) = filtfilt(b,a,ERPResponse(indChanComb, :, indEpoch,1, indPP)); % & obj.eeg.responseEpoch >= Threshold(2)
                            else
                                tmpBaseline = 0;
                            end
                            
                            plotWhole(indChanComb,:,indEpoch,indPP) = plotWhole(indChanComb,:,indEpoch, indPP) - tmpBaseline;
                            
                        end
                    end
                end
                
                if obj.DetectOrDisc && size(plotWhole, 1) > 1
                    
                    %                     tmpplotTarget    = nan(1, size(plotTarget,2), size(plotTarget,3));
                    %                     tmpplotResponse  = nan(1, size(plotResponse,2), size(plotResponse,3), size(plotResponse,4));
                    %                     for indEpoch = 1:length(obj.behaviour{indPP}.Reaction)
                    %                         if obj.behaviour{indPP}.Reaction(indEpoch,1) == 1 && obj.behaviour{indPP}.correct(indEpoch)
                    %                             indChannels = 2;
                    %                             tmpplotTarget(1,:,indEpoch) = plotTarget(indChannels,:,indEpoch, indPP);
                    %                             tmpplotResponse(1,:,indEpoch,1) = plotResponse(indChannels,:,indEpoch, 1, indPP);
                    %                         elseif (obj.behaviour{indPP}.Reaction(indEpoch,1) == 3 || obj.behaviour{indPP}.Reaction(indEpoch,1) == 2) && obj.behaviour{indPP}.correct(indEpoch)
                    %                             indChannels = 1;
                    %                             tmpplotTarget(1,:,indEpoch) = plotTarget(indChannels,:,indEpoch, indPP);
                    %                             tmpplotResponse(1,:,indEpoch,1) = plotResponse(indChannels,:,indEpoch, 1, indPP);
                    %                         end
                    %                     end
                    %
                    plotWhole(1,:,:, indPP)     = nanmean(plotWhole(:,:,:, indPP),1);
                end
                
            end
            if obj.DetectOrDisc && size(plotWhole, 1) > 1
                plotResponse(1,:,:,:,:)     = nanmean(plotResponse,1);
                plotWhole(2,:,:,:) = [];
                plotResponse(2,:,:,:,:)  = [];
            end
            
            % 1) first get 'remove between participant variablity'.
            % get the individual means and the grandAverage
            indivMeanWhole   = squeeze(nanmean(plotWhole,3))';
            grandMeanWhole   = squeeze(nanmean(nanmean(plotWhole,3),4));
            
            indivMeanResp    = squeeze(nanmean(plotResponse,3))';
            grandMeanResp   = squeeze(nanmean(nanmean(plotResponse(:,:,:,1,:),3),5));
            currColor = obj.figLayOut.colours;
            
            [within, between, modelFun, allTrialMatrix] = getConditions(obj, plotComb);
            
            for indPP = 1:length(obj.ppNames)
                clear remove*
                for indEpoch = 1:size(plotWhole,3)
                    removeWhole(1,:,indEpoch)     = plotWhole(:,:,indEpoch,indPP);%- indivMeanWhole(indPP,:)   + grandMeanWhole;
                    
                    if BaselineCorrect == 3
                        removeResponse(1,:,indEpoch)  = squeeze(plotResponse(:,:,indEpoch,1,indPP));% - indivMeanResp(indPP,:)   + grandMeanResp;
                    end
                end
                
                VartrialMatrix   = allTrialMatrix(:,:, indPP);
                for indCond = 1:within.numCond
                    currCond = sum(VartrialMatrix == within.Condition(indCond,:),2) == size(within.Condition, 2);
                    
                    % remove all conditions that didn't have a
                    % response, e.g. the nans.
                    
                    currWhole(:,:, indCond, indPP) = nanmean(removeWhole(:,:, currCond),3);
                    
                    
                    if BaselineCorrect == 3
                        currCond = VartrialMatrix(:,1) == within.Condition(indCond,1);
                        
                        currResponse(:,:, indCond, indPP) = nanmean(removeResponse(:,:, currCond),3);
                        Threshold(indCond, indPP) =   nanmean(nanmean(removeResponse(:, ...
                            obj.eeg.responseEpoch>= thresArea(1) & obj.eeg.responseEpoch <= thresArea(2),...
                            currCond),3),2);
                    else
                        Threshold(indCond, indPP) = 0;
                    end
                end
            end
            
            %% %%%%%%%%%%%%%%%%%%%%%    timeSeries    %%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot timeSeries seperately for fast/slow and no response for
            % both condition both Epoch locked at stimulus onset and
            % response.
            if grouping ~= 0
                groups = unique(within.Condition(:,1))';
            else
                groups = 1;
                
            end
            groups(isnan(groups)) = [];
            
            acc = 1;
            for indGroup = 1:length(groups)
                fig.Handle{indGroup} = figure('units','normalized','outerposition',[0 0 1 1]);
                % create panel
                
                fig.Info{indGroup} = panel();
                fig.Info{indGroup}.pack(1, 1);
                fig.Info{indGroup}(1,1).marginright = 4;
                
                fig.Info{indGroup}.marginbottom = 10;
                fig.Info{indGroup}.marginright  = 10;
                fig.Info{indGroup}.fontsize = obj.figLayOut.letterSize;
                fig.Info{indGroup}.fontname = obj.figLayOut.letterType;
                for indBetween = unique(between.Table)
                    
                    if grouping ~= 0
                        numCond = sum(within.Condition(:,1) == groups(indGroup));
                    else
                        numCond = size(within.Condition,1);
                    end
                    for indCond = 1:numCond
                        
                        cutOut = obj.eeg.epochPlot > max(obj.stim.lengthITI) & obj.eeg.epochPlot < obj.stim.RTdeadLine(end);
                        
                        plotTimeSeries = obj.eeg.epochPlot(cutOut);
                        meanWhole = (nanmean(currWhole(:,cutOut,acc, between.Table == indBetween),4)) - nanmean(Threshold(acc,:),2);
                        CIWhole   = (1.96.*(nanstd(currWhole(:,cutOut,acc, between.Table == indBetween), [], 4)/sqrt(sum(between.Table == indBetween))));
                        
                        %forSave(:,indCond) = squeeze(meanWhole);
                        YStuff(acc,:) = [min(meanWhole) max(meanWhole)];
                        
                        %% plot Target-locked epochs
                        figure(fig.Handle{indGroup}); hold on
                        
                        % plot evidence onset
                        h = shadedErrorBar(plotTimeSeries, squeeze(meanWhole), squeeze(CIWhole));
                        h.mainLine.Color = currColor(acc,:);  h.patch.FaceColor = currColor(acc,:);
                        h.edge(1).Color  = currColor(acc,:);  h.edge(2).Color   = currColor(acc,:);
                        h.mainLine.LineWidth = 1; h.edge(1).LineStyle = 'none'; h.edge(2).LineStyle = 'none';
                        h.patch.FaceAlpha = 0;
                        
                        h.mainLine.LineStyle = '--';
                        fig.legendThis(acc) = h.mainLine;
                        
                        startITI = -1 * (max(obj.stim.lengthITI)  - (obj.stim.RTdeadLine(end) - obj.stim.duration));
                        del = find(isnan(meanWhole) | plotTimeSeries > 0 | plotTimeSeries < startITI); %
                        tmpPlotTimeSeries =  plotTimeSeries;
                        tmpPlotTimeSeries(del)=[];
                        meanWhole(del)=[];
                        
                        b = polyfit(tmpPlotTimeSeries,meanWhole,2);
                        urg = b(:,1)*plotTimeSeries.^2+b(:,2)*plotTimeSeries+b(:,3);
                        urg(plotTimeSeries< startITI) = nan;%
                        urg(plotTimeSeries>0) = urg(find(plotTimeSeries>=0,1));
                        
                        plot(plotTimeSeries,urg,...
                            'Color', currColor(acc,:), 'LineWidth', obj.figLayOut.lineWidth)
                        
                        acc = acc + 1;
                        
                    end
                end
                
                if strfind(plotName, 'ERD') ~= 0
                    set ( gca, 'ydir', 'reverse' )
                end
            end
            
            %% layouting
            YStuff = [min(YStuff(:,1)) max(YStuff(:,2))];
            if obj.eeg.applyCSD
                limits = [floor(YStuff(:,1)/5)*5 ceil(YStuff(:,2)/5)*5];
                if limits(1) == 0; limits(1) = -5; end
                if limits(2) == 0; limits(2) = 5; end
            else
                limits = [floor(YStuff(:,1)/1)*5 ceil(YStuff(:,1)/5)*5];
                if limits(1) == 0; limits(1) = -1; end
                if limits(2) == 0; limits(2) = 1; end
            end
            
            
            for indFig = 1:length(fig.Handle)
                figure(fig.Handle{indFig})
                %fig.Info{indFig}(1, 1).select(); hold on;
                
                ylim([limits]);
                line([0 0], limits, 'Color', 'k', 'LineWidth', 1.5)
            end
            
            for indHandle = 1:length(fig.Handle)
                figure( fig.Handle{indHandle})
                
                if grouping ~= 0
                    NameForPlot = obj.figLayOut.legNames{grouping}{indHandle};
                else
                    NameForPlot = obj.figLayOut.legNames{plotComb(1)}(indHandle,:);
                end
                
                %                 plotSave(gca,  [plotName NameForPlot '.png'], averageFolder, obj.figLayOut.saveDim);
                %                 leg =  legend(fig.legendThis(1:length(obj.figLayOut.legNames{plotComb(1)})),...
                %                     obj.figLayOut.legNames{plotComb(1)}, 'Location', 'northwest');
                %                 leg.Title.String = obj.figLayOut.legTitle{plotComb(1)};
                %                 plotSave(gca,  [plotName NameForPlot '_Legend.png'], averageFolder, obj.figLayOut.saveDim);
            end
            
            fprintf('.\n')
            fprintf('Done, go get yourself a coffee :) \n')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %%%%%%%%%%%%%    supporting function  %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [perGoodTrials, ERPTarget, ERPResponse, ERPFalseAlarm] = sortERPs(obj, ERP, indPP, indWhich, artifactReject)
            %% processERPs(obj)
            % function is used to extract the epoch data looking for the
            % stimulus-locked as well as the reponse-locked epochs. This
            % requires to a ERP and indPP. And gives the segmented
            % ERPTarget and ERPResponse back. Add indWhich to indicate what
            % to plot.
            
            % Some preset parameters
            if ~exist('artifactReject', 'var'), artifactReject = 1; end
            if ~exist('indWhich', 'var'),       indWhich = 0; end
            
            % identify good trials to only be used!
            if artifactReject
                goodTrials = obj.behaviour{indPP}.artifacts & obj.behaviour{indPP}.blinks;
            else
                goodTrials = ones(length(obj.behaviour{indPP}.Reaction ),1);
            end
            perGoodTrials = sum(goodTrials)./length(goodTrials);
            
            % preset everything and allocated the space.
            trangeTarget    = obj.eeg.epochPlot >= obj.eeg.targetEpoch(1) & obj.eeg.epochPlot <= obj.eeg.targetEpoch(end);
            
            ERPTarget       = nan(size(ERP,1), sum(trangeTarget),  max(obj.numTrials*obj.numBlocks));
            
            ERPFalseAlarm   = nan(size(ERP,1),  length(obj.eeg.responseEpoch)+512,  max(obj.numTrials*obj.numBlocks), size(obj.behaviour{indPP}.indFalseAlarm,2));
            ERPResponse     = nan(size(ERP,1),  length(obj.eeg.responseEpoch),  max(obj.numTrials*obj.numBlocks));
            
            if perGoodTrials > 0.3 % conditional only when 30% of the trials are left.
                
                % check how many trials are still left and if this person
                % should be include
                hitsTrial = goodTrials';
                
                %% ------------ Extracting Target-locked ------------------
                if indWhich == 0 || indWhich == 1
                    ERPTarget(:, :,  any(hitsTrial,2)) = ERP(:, trangeTarget,  any(hitsTrial,2));
                end
                
                %% ------------ Extracting Response-locked ----------------
                
                if indWhich == 0  || indWhich == 2
                    for indEpoch = 1:length(hitsTrial)
                        if any(hitsTrial(indEpoch,:))
                            for indHit = find(hitsTrial(indEpoch,:) == 1)
                                currtrangeResponse = obj.behaviour{indPP}.indReaction(indEpoch,indHit) + obj.eeg.responseEpoch./(1/obj.eeg.SampleRate);
                                if currtrangeResponse(end) < size(ERP,2)
                                    ERPResponse(:, :, indEpoch, indHit) = ERP(:, currtrangeResponse, indEpoch);
                                end
                            end
                        end
                    end
                    
                    ERPResponse(repmat(all(ERPResponse == 0,2),1,size(ERPResponse,2),1,1)) = NaN;
                end
                
                %% ------------ Extracting False alarms -------------------
                
                if (indWhich == 0  || indWhich == 3) && ~isempty(obj.stim.timing)
                    FalseAlarms  = obj.behaviour{indPP}.indFalseAlarm;
                    for indEpoch = 1:length(FalseAlarms)
                        currFA = FalseAlarms(indEpoch,~isnan(FalseAlarms(indEpoch,:)));
                        
                        if isempty(currFA)
                            continue
                        end
                        
                        if ~isempty(obj.stim.timing)
                            for indFA = 1:length(currFA)
                                
                                if currFA(indFA) + obj.eeg.responseEpoch(1) < obj.eeg.epochPlot(1)
                                    continue
                                else
                                    
                                    currFAtrange = obj.eeg.epochPlot >= currFA(indFA)  + obj.eeg.responseEpoch(1) &...
                                        obj.eeg.epochPlot <= currFA(indFA) + obj.eeg.responseEpoch(end)+1;
                                    toCheckFalseAlarm = ERP(:,currFAtrange, indEpoch);
                                    try
                                        if ~(sum(~isnan(toCheckFalseAlarm),2) ~= size(ERPFalseAlarm(:,:,indEpoch, indFA),2))
                                            ERPFalseAlarm(:,:,indEpoch, indFA) = ERP(:, currFAtrange, indEpoch);
                                        end
                                    catch
                                        keyboard
                                    end
                                end
                                
                            end
                        end
                        
                    end
                end
            end
        end
        
        function [within, between, modelFun, allTrialMatrix] = getConditions(obj, plotComb, TimeOnTask, TimeOnExperiment)
            
            if ~exist('TimeOnTask', 'var'), TimeOnTask = 0; end
            if ~exist('TimeOnExperiment', 'var'), TimeOnExperiment = 0; end
            
            within.Var = {}; within.Name = {}; between.Var = {}; between.Name = {};  within.Par = [];
            acc = 1;
            for indVar = plotComb
                if indVar == 0
                    within.Var{end+1}  = unique(obj.behaviour{1}.binRT(~isnan(obj.behaviour{1}.binRT)))';
                    within.Name{end+1} = 'RTbins';
                    within.Par(acc) = 1;
                    %within.numCond = within.numCond*length(within.Var{end});
                    
                elseif length(unique(obj.behaviour{1}.trialMatrix(:,indVar))) > 1
                    within.Var{end+1}  = obj.conditions{indVar};
                    within.Name{end+1} = strrep(obj.figLayOut.legTitle{indVar},' ','');
                    within.Par(acc) = 1;
                    %within.numCond = within.numCond*length(within.Var{end});
                else
                    between.Var{end+1}  = obj.conditions{indVar};
                    between.Name{end+1} = strrep(obj.figLayOut.legTitle{indVar},' ','');
                    within.Par(acc) = 0;
                end
                acc = acc + 1;
            end
            
            if TimeOnTask
                within.Var{end+1} = [1 2]; %[earlier vs. later in a block e.g. numTrial/2]
                within.Name{end+1} = 'TimeOnTask';
                within.numCond = within.numCond*length(within.Var{end});
            end
            
            if TimeOnExperiment
                within.Var{end+1} = [1 2]; %[earlier vs. later in the Experiment e.g. max(numBlocks)/2]
                within.Name{end+1} = 'TimeOnExperiment';
                within.numCond = within.numCond*length(within.Var{end});
            end
            
            switch length(within.Var)
                case 1
                    modelFun = sprintf('%s',...
                        within.Name{1});  % main effects
                case 2
                    modelFun = sprintf('%s + %s + %s * %s',...
                        within.Name{1}, within.Name{2},... % main effects
                        within.Name{1}, within.Name{2});
                case 3
                    modelFun = sprintf('%s + %s + %s + %s * %s + %s * %s + %s * %s + %s * %s * %s',...
                        within.Name{1}, within.Name{2}, within.Name{3},... % main effects
                        within.Name{1}, within.Name{2},  within.Name{1}, within.Name{3},... % two-way interactions
                        within.Name{2}, within.Name{3},...
                        within.Name{1}, within.Name{2}, within.Name{3}); % three-way interactions
                case 4
                    modelFun = sprintf('%s + %s + %s + %s + %s * %s + %s * %s + %s * %s + %s * %s + %s * %s + %s * %s + %s * %s * %s  + %s * %s * %s  + %s * %s * %s  + %s * %s * %s * %s',...
                        within.Name{1}, within.Name{2}, within.Name{3}, within.Name{4},... % main effects
                        within.Name{1}, within.Name{2},  within.Name{1}, within.Name{3},  within.Name{1}, within.Name{4},... % two-way interactions
                        within.Name{2}, within.Name{3}, within.Name{2}, within.Name{4},...
                        within.Name{3}, within.Name{4},...
                        within.Name{1}, within.Name{2}, within.Name{3}, within.Name{1}, within.Name{2}, within.Name{4},...% three-way interactions
                        within.Name{2}, within.Name{3}, within.Name{4},...
                        within.Name{1}, within.Name{2}, within.Name{3},  within.Name{4}); % four-way interactions
                otherwise
                    keyboard
            end
            
            allTrialMatrix = nan(max(obj.numTrials.*obj.numBlocks), length(within.Name), length(obj.ppNames));
            forWithin = [];
            for indPP = 1:length(obj.ppNames)
                withinComb = plotComb(within.Par == 1);
                
                allTrialMatrix(1:length(obj.behaviour{indPP}.trialMatrix), withinComb ~= 0, indPP)  = obj.behaviour{indPP}.trialMatrix(:,withinComb(withinComb ~= 0));
                
                if any(withinComb == 0)
                    allTrialMatrix(1:length(obj.behaviour{indPP}.trialMatrix), withinComb == 0, indPP) = obj.behaviour{indPP}.binRT;
                end
                
                
                if TimeOnTask
                    tmpTimeOnTask = nan(length(allTrialMatrix),1);
                    tmpOrder      = obj.order{indPP};
                    for indBlock = 1:obj.numBlocks(indPP)
                        tmpBlock = find(tmpOrder == indBlock);
                        tmpTimeOnTask(tmpBlock(1:length(tmpBlock)/2)) = 1;
                        tmpTimeOnTask(tmpBlock(length(tmpBlock)/2+1:end)) = 2;
                    end
                    
                    allTrialMatrix(:, strcmp(within.Name, 'TimeOnTask')) = tmpTimeOnTask;
                    clear tmp*
                end
                
                if TimeOnExperiment
                    tmpOrder	= obj.order{indPP};
                    tmpEarlier  = tmpOrder <= (max(obj.numBlocks)./2);
                    tmpLater    = tmpOrder > (max(obj.numBlocks)./2);
                    tmpOrder(tmpEarlier) = 1; tmpOrder(tmpLater) = 2;
                    
                    allTrialMatrix(:,strcmp(within.Name, 'TimeOnExperiment')) = tmpOrder;
                    clear tmp*
                end
                
                if sum(~within.Par)
                    between.Table(indPP) = unique(obj.behaviour{indPP}.trialMatrix(:,find(~within.Par)));
                else
                    between.Table(indPP)  = 1;
                end
                forWithin = [forWithin; allTrialMatrix(:,:,indPP)];
            end
            
            within.Condition = unique(forWithin, 'rows');
            within.Condition(any(isnan(within.Condition),2),:) = [];
            within.Table = array2table(within.Condition, 'VariableNames', within.Name);
            
            within.numCond = size(within.Condition,1);
            
        end
        %% -----------------     plot functions     -----------------------
        % In this section you will find the different options to plot the
        % information. This includes individual and average timeseries and
        % topoplots. Topoplots can have several inputs. [TODO make the
        % Misses/False alarms optional.
        
        %----------------    Behavioural plots  ---------------------------
        function figInfo = plotAverageBehaviour(obj, plotName, rm, plotComb)
            %% plotAvaregeBehaviour
            % this function is used to flexibility make plots of the
            % behaviour and timeseries data. It simply requires a y
            % (dependent variable), x (independent variable that you want
            % to plot on the x axis), group (grouping e.g. different
            % colours) and sepPlot (if there are many we have to plot the
            % in subplots)
            
            % make the save folder
            figureFolder = fullfile(obj.figFolder,  'groupAverage/Behaviour/');
            if ~exist(figureFolder, 'dir')
                mkdir(figureFolder);
            end
            
            plotCombination = [];
            for indCond = 1:length(rm.WithinFactorNames)
                plotCombination = [plotCombination double(table2array(rm.WithinDesign(:,indCond)))];
            end
            
            for indCond = 1:length(rm.BetweenFactorNames)
                keyboard
                plotCombination = [plotCombination double(table2array(rm.BetweenDesign(:,indCond)))];
            end
            
            
            % get depending variable on the z-axis, group and seperated
            % plots but first check if they are significant. Delete if
            % no main effect or interaction effect.
            
            plotData = margmean(rm, cat(2,rm.WithinFactorNames, rm.BetweenFactorNames));
            plotData(isnan(plotData.Mean),:) = [];
            
            grouping = ones(size(plotCombination,1),1);
            if size(plotCombination,2) == 1
                plotCombination(:,end+1) = 1;
            elseif size(plotCombination) == 3
                grouping = unique(plotCombination(:,3));
            end
            
            groups = unique(grouping);
            
            plotSaveName = [plotName num2str(plotComb)];
            
            for indGroup = 1:length(groups)  % first set the xy/group/seperatedPlots
                figure('units','normalized','outerposition',[0 0 1 1]); hold on;
                currGroup = grouping == groups(indGroup);
                
                if length(groups) > 1
                    plotSaveName = [plotName num2str(plotComb) obj.figLayOut.legNames{plotComb(3)}{indGroup}];
                end
                
                [~,~,XLoc] = unique(plotCombination(currGroup,1));
                depVar    = plotData.Mean(currGroup,:);
                depVarCI  = (plotData.Upper(currGroup,:) - plotData.Lower(currGroup,:))/2;
                
                currCond = plotCombination(currGroup,2);
                disPlot = -0.1:0.2/(length(unique(currCond))-1):0.1;
                
                acc = 1;
                for indCond = 1:length(unique(currCond))
                    plotCond = plotCombination(:,2) == indCond;
                    
                    plotX = XLoc(plotCond)+ disPlot(acc);
                    meanY = depVar(plotCond);
                    yCI95 = depVarCI(plotCond);
                    
                    lineInfo(acc) = errorbar(plotX, meanY, yCI95);
                    lineInfo(acc).Color     = obj.figLayOut.colours(acc,:);
                    lineInfo(acc).LineWidth = obj.figLayOut.lineWidth;
                    
                    acc = acc + 1;
                end
                
                figInfo = gca;
                set(figInfo,'FontSize', obj.figLayOut.letterSize);
                set(figInfo,'FontName', obj.figLayOut.letterType);
                
                ylabel(plotName);
                
                if figInfo.YLim(2) > 150
                    ylim([0  figInfo.YLim(2)])
                    yticks(0:0.2:figInfo.YLim(2))
                elseif contains(plotName, '(\muV/m^2)')
                else
                    ylim([0 105])
                    yticks(0:20:100)
                end
                
                xlim([1 size(unique(XLoc),1)] + [-0.2 0.2])
                xlabel(obj.figLayOut.legTitle{plotComb(1)});
                xticks(1:length(obj.conditions{plotComb(1)}))
                xticklabels(obj.figLayOut.legNames{plotComb(1)});
                
                plotSave(figInfo, [plotSaveName '.png'], figureFolder, obj.figLayOut.saveDim);
                
                if length(plotComb) > 1
                    leg = legend(obj.figLayOut.legNames{plotComb(2)}, 'Location', 'southwest');
                    leg.Title.String = obj.figLayOut.legTitle{plotComb(2)};
                    plotSave(figInfo, [plotSaveName '_Legend.png'], figureFolder, obj.figLayOut.saveDim);
                end
            end
        end
        
        function [figHandle, onlyLegend] = plotRThistogram(obj, plotComb, grouping)
            %% plot the RT distributions
            figureFolder = fullfile(obj.figFolder,  'groupAverage/Behaviour/');
            if ~exist(figureFolder, 'dir')
                mkdir(figureFolder);
            end
            if ~exist('grouping', 'var'); grouping = 0; end
            
            %% set up model parameters
            % as there are different parameters that you might wanna add,
            % like time dependency. For now it just simple time on task and
            % time within the experiment.
            [within, between, modelFun, allTrialMatrix] = getConditions(obj, plotComb);
            
            % means per participant per conditions
            reactionTimes = nan(obj.numTrials*max(obj.numBlocks), length(obj.ppNames), within.numCond);
            correct       = nan(obj.numTrials*max(obj.numBlocks), length(obj.ppNames), within.numCond);
            errors        = nan(obj.numTrials*max(obj.numBlocks), length(obj.ppNames), within.numCond);
            
            for indPP = 1:length(obj.ppNames)
                
                posCombination = unique(allTrialMatrix(:,:, indPP), 'rows');
                posCombination(any(isnan(posCombination),2),:) = [];
                
                RT          = reshape(obj.behaviour{indPP}.RT,[],1);
                
                for indCond = 1:length(posCombination)
                    currCond  = sum(allTrialMatrix(1:length(RT),:,indPP) == posCombination(indCond,:),2) == size(posCombination, 2) & ~isnan(RT);
                    currRT	  = RT(currCond,:);
                    currCondHit = ~isnan(currRT);
                    
                    reactionTimes(1:length(currRT),indPP,indCond) = currRT;
                    
                    if obj.DetectOrDisc
                        correct(1:length(currRT),indPP,indCond) = obj.behaviour{indPP}.Correct(currCondHit,:) & currRT > obj.stim.RTCutOff & currRT < obj.stim.RTdeadLine(end);
                        errors(1:length(currRT),indPP,indCond)  = ~(obj.behaviour{indPP}.Correct(currCondHit,:)) & currRT > obj.stim.RTCutOff & currRT < obj.stim.RTdeadLine(end);
                    else
                        correct(1:length(currRT),indPP,indCond) = 1;
                    end
                end
            end
            
            correct(isnan(correct)) = 0; errors(isnan(errors)) = 0;
            numResponse = sum(~isnan(reactionTimes(:)));
            
            % put seperately.
            timeSeries = 0:0.050:obj.stim.RTdeadLine(end);
            
            acc = 1; forColor = 1; figAxis = [];
            if grouping ~= 0
                groups = unique(posCombination(:,plotComb == grouping));
            else
                groups = 1;
            end
            
            for indGroup = 1:length(groups)
                if grouping ~= 0
                    currInformation = unique(posCombination(posCombination(:, plotComb == grouping ) == groups(indGroup),:), 'rows');
                else
                    currInformation = posCombination;
                end
                
                figHandle{indGroup} = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
                for indBetween = 1:length(unique(between.Table))
                    for indCond = 1:size(currInformation,1)
                        currCond    =  find( sum(posCombination == currInformation(indCond,:),2) == size(currInformation, 2));
                        currRT      = reactionTimes(:,between.Table == indBetween,currCond);
                        currCorrect = currRT(logical(correct(:,between.Table == indBetween,currCond)));
                        
                        x = hist(currCorrect,timeSeries);
                        xgrid  = linspace(timeSeries(1),obj.eeg.targetEpoch(end),length(x))';
                        
                        if any(plotComb == 0)
                            onlyLegend(indCond) = line(xgrid,x/numResponse.*100, 'Color', obj.figLayOut.colours(indCond,:), 'LineWidth', 2, 'LineStyle', '-.');
                        else
                            onlyLegend(indCond) = line(xgrid,x/numResponse.*100, 'Color', obj.figLayOut.colours(indCond,:), 'LineWidth', 2, 'LineStyle', '-');
                        end
                        
                        if obj.DetectOrDisc
                            currError = currRT(logical(errors(:,between.Table == indBetween,indCond)));
                            x      = hist(currError,timeSeries);
                            xgrid  = linspace(timeSeries(1),obj.eeg.targetEpoch(end),length(x))';
                            line(xgrid, x/numResponse.*100, 'Color', obj.figLayOut.colours(indCond,:), 'LineWidth', 1, 'LineStyle', ':');
                        end
                        meanRT(acc) = nanmedian(currCorrect); %
                        acc = acc +1;
                    end
                    tmpAxis = gca;
                    figAxis(end+1) = tmpAxis.YLim(2);
                end
            end
            
            figAxis = max(figAxis);
            acc = 0;
            for indGroup = 1:length(figHandle)
                currInformation( any(isnan( currInformation),2),:) = [];
                figure(figHandle{indGroup})
                figInfo = gca;
                line([0 0],  [0 figAxis], 'Color', 'k', 'LineWidth', 1.5)
                
                ylim([0 figAxis])
                ylabel('Proportion (%)')
                
                xstuff = [fliplr(0:-0.5:obj.eeg.targetEpoch(1)) 0.5:0.5:obj.eeg.targetEpoch(end)];
                
                xlim([obj.eeg.targetEpoch(1) obj.eeg.targetEpoch(end)])
                xticks(xstuff)
                xticklabels(xstuff);
                xlabel('RT relative to target onset (sec)')
                
                set(gca,'FontSize', obj.figLayOut.letterSize);
                set(gca,'FontName', obj.figLayOut.letterType);
                plotSave(gcf, ['RTPlot' num2str(indGroup) '.png'], figureFolder, [obj.figLayOut.saveDim]);
                
                leg = legend(onlyLegend, num2str(posCombination), 'Location', 'northwest');
                if length(plotComb) == 1
                    leg.Title.String = obj.figLayOut.legTitle{plotComb};
                end
                plotSave(gcf, ['RTPlot' num2str(indGroup) num2str(plotComb) 'Legend.png'], figureFolder, [obj.figLayOut.saveDim]);
            end
        end
        
        %----------------    Individual EEG plots.   ----------------------
        function figHandle = plotIndividualTimeSeries(obj, indPP, Target, Response, plotName, plotComb, grouping)
            % plot individual timeseries. This can either be by reaction
            % times (0) or choicing an 'extra' condition. Outerwise it uses
            % the first condition in obj.conditions to group the data.
            if ~exist('plotComb', 'var'); plotComb = []; end
            if ~exist('grouping', 'var'); grouping = 0; end
            
            currColor     = obj.figLayOut.colours;
            
            
            [within, between, ~, allTrialMatrix] = getConditions(obj, plotComb);
            
            if grouping ~= 0
                groups = unique(within.Condition(:,grouping))';
            else
                groups = 1;
            end
            groups(isnan(groups)) = [];
            
            for indGroup = 1:length(groups)
                figHandle{indGroup} = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
                
                % 2) split trials up according to the withinMatrix.
                posCombination = unique(allTrialMatrix(:,:, indPP), 'rows');
                posCombination(any(isnan(posCombination),2),:) = [];
                
                if grouping ~= 0
                    numCond = sum(within.Condition(:,grouping) == groups(indGroup));
                else
                    numCond = size(within.Condition,1);
                end
                
                %% plot individual timeseries.
                for indCond = 1:numCond
                    % plot evidence onset
                    subplot(1,3,1); hold on;
                    currCond = sum(allTrialMatrix(:,:,indPP) == posCombination(indCond,:),2) == size(posCombination,2);
                    meanPlot = nanmean(Target(:,:,currCond),3);
                    SMEPlot  = nanstd(Target(:,:,currCond), [],3)./sqrt(sum(currCond));
                    
                    currPlot = shadedErrorBar(obj.eeg.targetEpoch, meanPlot, SMEPlot);
                    currPlot.mainLine.Color  = currColor(indCond,:); currPlot.mainLine.LineWidth = 2;
                    
                    if sum(plotComb == 0) && posCombination(indCond,1) == 2
                        currPlot.mainLine.LineStyle = '-.';
                    end
                    
                    currPlot.patch.FaceColor = currColor(indCond,:); currPlot.edge(1).LineStyle = 'none'; currPlot.edge(2).LineStyle = 'none';
                    legendThis(indCond) = currPlot.mainLine;
                    
                    % plot response onset
                    subplot(1,3,2); hold on;
                    currCond = sum(allTrialMatrix(:,:,indPP)  == posCombination(indCond,:),2)  == size(posCombination,2);
                    meanPlot = nanmean(nanmean(Response(:,:,currCond,:),3),4);
                    SMEPlot  = nanstd(reshape(Response(:,:,currCond,:),size(Response,1),size(Response,2),[]), [],3)./sqrt(sum(currCond));
                    
                    currPlot = shadedErrorBar(obj.eeg.responseEpoch, meanPlot, SMEPlot);
                    currPlot.mainLine.Color  = currColor(indCond,:); currPlot.mainLine.LineWidth = 2;
                    if sum(plotComb == 0) && posCombination(indCond,1) == 2
                        currPlot.mainLine.LineStyle = '-.';
                    end
                    currPlot.patch.FaceColor = currColor(indCond,:); currPlot.edge(1).LineStyle   = 'none'; currPlot.edge(2).LineStyle = 'none';
                end
                
                % Make it look nice! (Evidence onset)
                subplot(1,3,2); hold on; figInfo2 = gca; h = findobj(gca,'Type','line');
                subplot(1,3,1); hold on; figInfo1 = gca;
                
                % calculate the ylim
                yMaxPlots = max([figInfo2.YLim(2), figInfo1.YLim(2)]);
                yMinPlots = min([figInfo2.YLim(1), figInfo1.YLim(1)]);
                
                % add lines, e.g. onset and mean fast and slow RTs
                
                line([0 0], [yMinPlots yMaxPlots], 'Color', 'k', 'LineWidth', 2)
                
                if sum(plotComb == 0)
                    lineSpecs = {'-', '-.'};
                    for indCond = 1:size(unqCond,1)
                        line([nanmean(meanRTs(indCond,:)) nanmean(meanRTs(indCond,:))], currAx.YLim, 'Color', obj.figLayOut.colours(indCond,:),...
                            'LineStyle', lineSpecs{unqCond(indCond,1)} ,'LineWidth', obj.figLayOut.lineWidth)
                    end
                end
                
                if obj.DetectOrDisc
                    line([nanmedian(obj.behaviour{indPP}.RT(:))  nanmedian(obj.behaviour{indPP}.RT(:))], [yMinPlots yMaxPlots], 'Color', 'k', 'LineWidth', 1);
                elseif any(plotComb == 0)
                    line([nanmedian(obj.behaviour{indPP}.RT(obj.behaviour{indPP}.binRT == 1)) nanmedian(obj.behaviour{indPP}.RT(obj.behaviour{indPP}.binRT == 1))], [yMinPlots yMaxPlots], 'Color', 'k', 'LineWidth', 1);
                    line([nanmedian(obj.behaviour{indPP}.RT(obj.behaviour{indPP}.binRT == 2))  nanmedian(obj.behaviour{indPP}.RT(obj.behaviour{indPP}.binRT == 2))], [yMinPlots yMaxPlots], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.');
                end
                
                title('Evidence onset')
                % add axes and shape
                xlim([obj.eeg.targetEpoch(1) obj.eeg.targetEpoch(end)])
                ylim([yMinPlots yMaxPlots])
                p = get(gca, 'Position');
                p(4) = p(4) - (p(4) * 0.5);
                p(2) = p(2) + p(4)/2;
                set(gca, 'Position', p);
                if obj.eeg.applyCSD
                    ylabel([plotName ' Amplitude (\muV/m^2)'])
                else
                    ylabel([plotName ' Amplitude (\muV)'])
                end
                
                % Make it look nice! (Response locked)
                subplot(1,3,2); hold on; figInfo2 = gca;
                
                % add lines, e.g. onset and mean fast and slow RTs
                line([0 0], [yMinPlots yMaxPlots], 'Color', 'k', 'LineWidth', 2)
                % add axes and shape
                xlim([-0.6 0.2])
                xticks([-0.6 0 0.2])
                ylim([yMinPlots yMaxPlots])
                figInfo2.YAxis.Visible = 'off';   % remove y-axis
                
                p = get(gca, 'Position');
                p(4) = p(4) - (p(4) * 0.5);
                p(2) = p(2) + p(4)/2;
                
                p_diff = p(3) * (length(-0.6:0.2)/diff([obj.eeg.targetEpoch(1) obj.eeg.targetEpoch(end)]));
                p(3) = p(3) - p_diff;
                p(1) = p(1) + p_diff/2;
                set(gca, 'Position', p);
                title('Response onset')
            end
        end
        
        function plotIndividualTopo(obj, topoData, channels)
            %% plot topoplots of hits and false alarms
            titleNames = {'Topoplots of hits'};
            
            %% Topoplots
            % 1) create topoplot hits
            subplot(size(topoData,1),3,3)
            topoFigInd = topoplot(double(nanmean(topoData(1,:,:),3)),obj.eeg.chanlocs, 'electrodes', 'labels');
            title(titleNames)
            
            if ~isempty(channels)
                ssvepChanPlot = {obj.eeg.ChannelsName{channels(~isnan(channels))}};
                keep = [];
                for kk = 1:length(topoFigInd.Parent.Children)
                    for indChan = 1:length(ssvepChanPlot)
                        try
                            tmp = topoFigInd.Parent.Children(kk).String(find(~isspace(topoFigInd.Parent.Children(kk).String)));
                            
                            if ~strcmpi(tmp, [ssvepChanPlot{indChan}]) %length([ssvepChanPlot{indChan}])
                                keep(kk, indChan) = 0;
                            else
                                keep(kk, indChan) = 1;
                            end
                        end
                    end
                end
                
                topoFigInd.Parent.Children( find(~sum(keep,2))').delete;
            end
            
        end
        
        %----------------    Averages EEG plot    -------------------------
        function fig = plotAverageTimeseries(obj, Target, Response, falseAlarms, baselineCorrect, plotName, plotComb, grouping, plotThis)
            %% plotAllAverage(obj)
            % this function is to replace the startChan that are used in the
            % plotERP to be able to get individualized the best electrode
            % per person. Based on average CPP WITH responses, also using the
            % CPP in the epoch locked on the response.
            
            % This can either be by reaction times (0) or choicing an 'extra' condition.
            % Outerwise it uses the first condition in obj.conditions to group the data.
            if ~exist('plotComb', 'var'); plotComb = []; end
            
            
            %% set up model parameters
            % as there are different parameters that you might wanna add,
            % like time dependency. For now it just simple time on task and
            % time within the experiment.
            
            [within, between, ~, allTrialMatrix] = getConditions(obj, plotComb);
            
            %% preallocated all the matrixes
            % extracted the means per participant per conditions
            meanTarget    = nan(size(Target, 1), size(Target, 2), within.numCond, length(obj.ppNames));
            meanMisses    = nan(size(Target, 1), size(Target, 2), within.numCond, length(obj.ppNames));
            meanResponse  = nan(size(Response, 1), size(Response, 2), within.numCond, length(obj.ppNames));
            meanFA        = nan(size(falseAlarms, 1), size(falseAlarms, 2), within.numCond, length(obj.ppNames));
            
            % extracted the means per participant per conditions
            stdTarget    = nan(size(Target, 1), size(Target, 2), within.numCond, length(obj.ppNames));
            stdMisses    = nan(size(Target, 1), size(Target, 2), within.numCond, length(obj.ppNames));
            stdResponse  = nan(size(Response, 1), size(Response, 2), within.numCond, size(Response, 4), length(obj.ppNames));
            stdFA        = nan(size(falseAlarms, 1), size(falseAlarms, 2), within.numCond, length(obj.ppNames));
            
            allFA     = zeros(within.numCond,length(obj.ppNames));
            allMisses = zeros(within.numCond,length(obj.ppNames));
            
            faEpoch = obj.eeg.epochPlot (obj.eeg.epochPlot >= obj.eeg.responseEpoch(1) & obj.eeg.epochPlot <= obj.eeg.responseEpoch(end)+1);
            
            %% for-loop to extract mean/std for plotting
            for indPP = 1:length(obj.ppNames)
                
                if plotThis == 6
                    % plot first derivative. Specificily interesting for CPP
                    % or any potential decision-making signals as it can
                    % show potential 'leaky accumulation'.
                    Target(:,:,:,indPP)        = [diff(Target(:,:,:,indPP)) nan(1,1,size(Target,3))];
                    Response(:,:,:,indPP)      = [diff(Response(:,:,:,:,indPP))  nan(1,1,size(Response,3),size(Response,4))];
                    falseAlarms(:,:,:,:,indPP) = [diff(falseAlarms(:,:,:,:,indPP))  nan(1,1,size(falseAlarms,3),size(falseAlarms,4))];
                end
                
                % 1) start of with removing the mean to get the standard
                % deviations, e.g. this will be a different Matrix.
                forStdTarget    = Target(:,:,:,indPP)        - nanmean(reshape(Target(:,:,:,indPP),[],1));
                forStdResponse  = Response(:,:,:,:,indPP)    - nanmean(reshape(Response(:,:,:,:,indPP),[],1));
                forStdFA        = falseAlarms(:,:,:,:,indPP) - nanmean(reshape(falseAlarms(:,:,:,:,indPP),[],1));
                
                % 2) split trials up according to the withinMatrix.
                posCombination = unique(allTrialMatrix(:,:, indPP), 'rows');
                posCombination(any(isnan(posCombination),2),:) = [];
                
                for indCond  = 1:size(posCombination,1)
                    
                    currCond = sum(allTrialMatrix(1:length(obj.behaviour{indPP}.RT),:,indPP) == posCombination(indCond,:),2) == size(posCombination, 2);
                    
                    currRT = reshape(obj.behaviour{indPP}.RT(currCond,:),[],1);
                    meanRTs(indCond) = nanmedian(currRT);
                    
                    % remove all conditions that didn't have a
                    % response, e.g. the nans.
                    
                    meanResponse(:,:, indCond, indPP) = nanmean(reshape(Response(:,:,currCond,:, indPP), size(Response,1),size(Response,2),[]),3);
                    
                    if baselineCorrect == 2
                        tmpBaseline = nanmean(meanResponse(:,obj.eeg.responseEpoch >= obj.eeg.baseline(1) & obj.eeg.responseEpoch <= obj.eeg.baseline(2),indCond,indPP),2);
                    else
                        tmpBaseline = 0;
                    end
                    
                    meanResponse(:,:, indCond, indPP) = nanmean(reshape(Response(:,:,currCond,:, indPP), size(Response,1),size(Response,2),[]) - tmpBaseline,3);
                    stdResponse(:,:, indCond, indPP)  = nanstd(reshape(forStdResponse(:,:,currCond,:), size(Response,1),size(Response,2),[]) - tmpBaseline,[],3);
               
                    
                    if plotThis == 2 || plotThis == 4
                        currCondHits = sum(allTrialMatrix(1:length(obj.behaviour{indPP}.RT),:,indPP) == posCombination(indCond,:),2) ==...
                            size(posCombination, 2) &...
                            ~obj.behaviour{indPP}.Misses;
                        meanTarget(:,:, indCond, indPP)	  = nanmean(Target(:,:, currCondHits, indPP) - tmpBaseline,3);
                        stdTarget(:,:, indCond, indPP)	  = nanstd(forStdTarget(:,:, currCondHits) - tmpBaseline,[],3);
                        
                        currCondMisses = sum(allTrialMatrix(1:length(obj.behaviour{indPP}.RT),:,indPP) == posCombination(indCond,:),2) == size(posCombination, 2) &...
                            obj.behaviour{indPP}.Misses;
                        
                        if sum(currCondMisses) > 5
                            meanMisses(:,:, indCond, indPP)	  = nanmean(Target(:,:, currCondMisses, indPP) - tmpBaseline,3);
                            stdMisses(:,:, indCond, indPP)	  = nanstd(forStdTarget(:,:, currCondMisses) - tmpBaseline,[],3);
                        end
                        allMisses(indCond, indPP) = sum(currCondMisses);
                    else
                        meanTarget(:,:, indCond, indPP)	  = nanmean(Target(:,:, currCond, indPP) - tmpBaseline,3);
                        stdTarget(:,:, indCond, indPP)	  = nanstd(forStdTarget(:,:, currCond) - tmpBaseline,[],3);
                    end
                    
               
                        
                    allFA(indCond, indPP) = allFA(indCond) + sum(sum(obj.behaviour{indPP}.FalseAlarm(currCond,:)));
                    
                    if sum( sum(sum(obj.behaviour{indPP}.FalseAlarm(currCond,:)))) > 1
                        try
                            if baselineCorrect
                                forBaseline = faEpoch > (nanmean(currRT)*-1)-100 & faEpoch < (nanmean(currRT)*-1);
                                meanFA(:,:, indCond, indPP) = nanmean(reshape(falseAlarms(:,:,currCond,:,indPP) - repmat(nanmean(falseAlarms(:,forBaseline,currCond,:,indPP),2), 1, size(falseAlarms,2),1,1),size(falseAlarms,2),[]),2);
                            else
                                meanFA(:,:, indCond, indPP) = nanmean(reshape(falseAlarms(:,:,currCond,:,indPP), size(falseAlarms,1),size(falseAlarms,2),[]),3);
                            end
                            
                            stdFA(:,:, indCond, indPP)  = nanstd(reshape(forStdFA(:,:,currCond,:), size(falseAlarms,1),size(falseAlarms,2),[]),[],3);
                        end
                    end
                end
            end
            
            %% %%%%%%%%%%%%%%%%%%%%%    timeSeries    %%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot timeSeries seperately for fast/slow and no response for
            % both condition both Epoch locked at stimulus onset and
            % response.
            if ~isempty(obj.figLayOut.targetLim)
                xstuff.Target = obj.figLayOut.targetLim;
            else
                xstuff.Target = [fliplr(0:-0.5:obj.eeg.targetEpoch(1)) 0.5:0.5:obj.eeg.targetEpoch(end)];
            end
            
            if ~isempty(obj.figLayOut.responseLim)
                xstuff.Response = obj.figLayOut.responseLim;
            else
                xstuff.Response = [-0.5 0 0.2];
            end
            
            if grouping ~= 0
                groups = unique(within.Condition(:,grouping))';
            else
                groups = 1;
            end
            groups(isnan(groups)) = [];
            
            acc = 1;
            for indGroup = 1:length(groups)
                
                fig.Handle{indGroup} = figure('units','normalized','outerposition',[0 0 1 1]);
                
                % create panel, when putting response and target-locked
                % together.
                if plotThis < 4
                    fig.Info{indGroup} = panel();
                    fig.Info{indGroup}.pack(1, {65 35});
                    fig.Info{indGroup}(1,1).marginright = 4;
                    fig.Info{indGroup}(1,2).marginleft  = 4;
                    
                    fig.Info{indGroup}.marginbottom = 10;
                    fig.Info{indGroup}.marginright  = 10;
                    fig.Info{indGroup}.fontsize = obj.figLayOut.letterSize;
                    fig.Info{indGroup}.fontname = obj.figLayOut.letterType;
                else
                    fig.ResponseHandle{indGroup} = figure('units','normalized','outerposition',[0 0 1 1]);
                end
                
                if obj.stim.FA && plotThis == 2
                    fig.FAHandle{indGroup} = figure('units','normalized','outerposition',[0 0 1 1]);
                    fig.FAInfo{indGroup} = panel();
                    fig.FAInfo{indGroup}.pack(1, {37 63});
                    fig.FAInfo{indGroup}(1,1).marginright = 4;
                    fig.FAInfo{indGroup}(1,2).marginleft  = 4;
                    
                    fig.FAInfo{indGroup}.marginbottom = 10;
                    fig.FAInfo{indGroup}.marginright  = 10;
                    fig.FAInfo{indGroup}.fontsize = obj.figLayOut.letterSize;
                    fig.FAInfo{indGroup}.fontname = obj.figLayOut.letterType;
                end
                
                if plotThis == 2
                    fig.MissesHandle{indGroup} = figure('units','normalized','outerposition',[0 0 1 1]);
                    fig.MissesInfo{indGroup} = panel();
                    fig.MissesInfo{indGroup}.pack(1, {65 35});
                    fig.MissesInfo{indGroup}(1,1).marginright = 4;
                    fig.MissesInfo{indGroup}(1,2).marginleft  = 4;
                    
                    fig.MissesInfo{indGroup}.marginbottom = 10;
                    fig.MissesInfo{indGroup}.marginright  = 10;
                    fig.MissesInfo{indGroup}.fontsize = obj.figLayOut.letterSize;
                    fig.MissesInfo{indGroup}.fontname = obj.figLayOut.letterType;
                end
                
                for indBetween = unique(between.Table)
                    if grouping ~= 0
                        numCond = sum(within.Condition(:,grouping) == groups(indGroup));
                    else
                        numCond = size(within.Condition,1);
                    end
                    
                    for indCond = 1:numCond
                        % this allows us to check when it is signifanctly
                        % different from baseline, can be use in the plots
                        % but right now not
                        
                        plotTarget   = nanmean(meanTarget(:,:,acc, between.Table == indBetween),4);
                        plotResponse = nanmean(meanResponse(:,:,acc, between.Table == indBetween),4);
                        plotFA	     = nanmean(meanFA(:,:,acc, between.Table == indBetween),4);
                        
                        CITarget   = (nanstd(stdTarget(:,:,acc, between.Table == indBetween), [], 4)/sqrt(sum(between.Table == indBetween)));
                        CIResponse = (nanstd(stdResponse(:,:,acc, between.Table == indBetween), [], 4)/sqrt(sum(between.Table == indBetween)));
                        CIFA	   = (nanstd(stdFA(:,:,acc, between.Table == indBetween), [], 4)/sqrt(sum(between.Table == indBetween)));
                        
                        if plotThis == 2
                            plotMisses   = nanmean(meanMisses(:,:,acc, between.Table == indBetween),4);
                            CIMisses   = (nanstd(stdMisses(:,:,acc, between.Table == indBetween), [], 4)/sqrt(sum(between.Table == indBetween)));
                        end
                        
                        if plotThis ~= 4 &&  plotThis ~= 5
                            YStuff(acc,1) = min([plotTarget-CITarget plotResponse-CIResponse]);
                            YStuff(acc,2) = max([plotTarget+CITarget plotResponse+CIResponse]);
                        else
                            YStuff(acc,1) = min([plotTarget-CITarget]);
                            YStuff(acc,2) = max([plotTarget+CITarget]);
                        end
                        
                        %% plot Target-locked epochs
                        figure(fig.Handle{indGroup})
                        if plotThis < 4
                            fig.Info{indGroup}(1, 1).select(); hold on;
                        end
                        
                        % plot evidence onset
                        h = shadedErrorBar(obj.eeg.targetEpoch, plotTarget, squeeze(CITarget));
                        h.mainLine.Color = obj.figLayOut.colours(acc,:);  h.patch.FaceColor = obj.figLayOut.colours(acc,:);
                        h.edge(1).Color  = obj.figLayOut.colours(acc,:);  h.edge(2).Color   = obj.figLayOut.colours(acc,:);
                        h.mainLine.LineWidth = obj.figLayOut.lineWidth; h.edge(1).LineStyle = 'none'; h.edge(2).LineStyle = 'none';
                        h.patch.FaceAlpha = obj.figLayOut.plotCI;
                        
                        h.mainLine.LineStyle = obj.figLayOut.lineType{acc};
                        if any(plotComb == 0) && posCombination(mod(indCond-1, 4)+1, 1) == 2
                            h.mainLine.LineStyle = '-.';
                        end
                        legendThis(acc) = h.mainLine;
                        grid on
                        
                        xlim([xstuff.Target(1) xstuff.Target(end)])
                        xticks(xstuff.Target)
                        xlabel ('Time from target (s)')
                        ax1 = gca;
                        
                        if obj.eeg.applyCSD
                            ylabel([plotName ' Amplitude (\muV/m^2)'])
                        else
                            ylabel([plotName ' Amplitude (\muV)'])
                        end
                        
                        if ~isempty(obj.figLayOut.plotRT) ||  sum(plotComb == 0)
                            line([nanmean(meanRTs(acc)) nanmean(meanRTs(acc))], ax1.YLim, 'Color', obj.figLayOut.colours(acc,:),...
                                'LineStyle', obj.figLayOut.lineType{acc},'LineWidth', obj.figLayOut.lineWidth)
                        end
                        %% plot response-locked epochs
                        if plotThis < 4
                            fig.Info{indGroup}(1, 2).select(); hold on;
                        else
                            figure(fig.ResponseHandle{indGroup})
                            
                        end
                        % plot evidence onset
                        h = shadedErrorBar(obj.eeg.responseEpoch, plotResponse, CIResponse);
                        h.mainLine.Color = obj.figLayOut.colours(acc,:);  h.patch.FaceColor = obj.figLayOut.colours(acc,:);
                        h.mainLine.LineWidth = obj.figLayOut.lineWidth;
                        h.edge(1).Color = obj.figLayOut.colours(acc,:);  h.edge(2).Color = obj.figLayOut.colours(acc,:);
                        h.edge(1).LineStyle = 'none'; h.edge(2).LineStyle = 'none';
                        h.patch.FaceAlpha = obj.figLayOut.plotCI;
                        
                        h.mainLine.LineStyle = obj.figLayOut.lineType{acc};
                        if any(plotComb == 0) && posCombination(mod(indCond-1, 4)+1, 1) == 2
                            h.mainLine.LineStyle = '-.';
                        end
                        
                        ax2 = gca;
                        
                        linkaxes([ax1 ax2], 'y')
                        
                        grid on
                        xlabel ('from response(s)')
                        xlim([ xstuff.Response(1)  xstuff.Response(end)])
                        xticks( xstuff.Response)
                        
                        %% plot FA-locked epochs
                        if obj.stim.FA && plotThis == 2
                            figure(fig.FAHandle{indGroup})
                            fig.FAInfo{indGroup}(1, 2).select(); hold on;
                            
                            if all(~isnan(plotFA))
                                h = shadedErrorBar(faEpoch, plotFA, CIFA);
                                h.mainLine.Color = obj.figLayOut.colours(acc,:);  h.patch.FaceColor = obj.figLayOut.colours(acc,:);
                                h.mainLine.LineWidth = obj.figLayOut.lineWidth;
                                h.edge(1).Color = obj.figLayOut.colours(acc,:);  h.edge(2).Color = obj.figLayOut.colours(acc,:);
                                h.edge(1).LineStyle = 'none'; h.edge(2).LineStyle = 'none';
                                h.patch.FaceAlpha = 0;%obj.figLayOut.plotCI;
                                h.mainLine.LineStyle = '-.';
                            end
                            xlim([-0.8 1])
                            xticks([-0.8 0:0.4:1])
                            grid on
                            xlabel ('before false alarm(s)')
                        end
                        
                        %% plot Target-locked epochs
                        if plotThis == 2 && all(~isnan(plotMisses))
                            
                            figure(fig.MissesHandle{indGroup})
                            fig.MissesInfo{indGroup}(1, 1).select(); hold on;
                            
                            % plot evidence onset
                            h = shadedErrorBar(obj.eeg.targetEpoch, plotMisses, squeeze(CIMisses));
                            h.mainLine.Color = obj.figLayOut.colours(acc,:);  h.patch.FaceColor = obj.figLayOut.colours(acc,:);
                            h.edge(1).Color  = obj.figLayOut.colours(acc,:);  h.edge(2).Color   = obj.figLayOut.colours(acc,:);
                            h.mainLine.LineWidth = obj.figLayOut.lineWidth; h.edge(1).LineStyle = 'none'; h.edge(2).LineStyle = 'none';
                            h.patch.FaceAlpha = obj.figLayOut.plotCI;
                            
                            h.mainLine.LineStyle = obj.figLayOut.lineType{acc};
                            if any(plotComb == 0) && posCombination(mod(indCond-1, 4)+1, 1) == 2
                                h.mainLine.LineStyle = '-.';
                            end
                            legendThis(acc) = h.mainLine;
                            
                            grid on
                            tmpXstuff = [fliplr(0:-0.5:obj.eeg.targetEpoch(1)) 0.5:0.5:obj.eeg.targetEpoch(end)];
                            xlim([obj.eeg.targetEpoch(1) obj.eeg.targetEpoch(end)])
                            xticks(tmpXstuff)
                            xlabel ('Time from target (s)')
                            ax1 = gca;
                            
                            if obj.eeg.applyCSD
                                ylabel([plotName ' Amplitude (\muV/m^2)'])
                            else
                                ylabel([plotName ' Amplitude (\muV)'])
                            end
                        end
                        acc = acc + 1;
                        
                    end
                end
                
            end
            
            %% layouting
            YStuff = [min(YStuff(:,1)) max(YStuff(:,2))];
            if obj.eeg.applyCSD & plotThis ~= 5 & any(abs(YStuff) > 10)
                limits = [floor(YStuff(:,1)/5)*5 ceil(YStuff(:,2)/5)*5];
                if limits(1) == 0; limits(1) = -5; end
                if limits(2) == 0; limits(2) = 5; end
            else
                limits = [floor(YStuff(:,1)/0.1)*0.1 ceil(YStuff(:,2)/0.1)*0.1];
                if limits(1) == 0; limits(1) = -1; end
                if limits(2) == 0; limits(2) = 1; end
            end
            
            for indFig = 1:length(fig.Handle)
                figure(fig.Handle{indFig})
                
                set(gca,'FontSize', obj.figLayOut.letterSize);
                set(gca,'FontName', obj.figLayOut.letterType);
                
                if plotThis < 4
                    fig.Info{indFig}(1, 1).select(); hold on;
                end
                
                ylim([limits]);
                line([0 0], limits, 'Color', 'k', 'LineWidth', 2)
                if plotThis < 4
                    
                    fig.Info{indFig}(1, 2).select(); hold on;
                    currAx = gca;
                    currAx.YAxis.Visible = 'off';
                    ylim(limits);
                    line([0 0], limits, 'Color', 'k', 'LineWidth', 2)
                end
                
                if plotThis == 2
                    if obj.stim.FA
                        figure(fig.FAHandle{indFig})
                        fig.FAInfo{indFig}(1, 2).select(); hold on;
                        currAx = gca;
                        %                         currAx.YAxis.Visible = 'off';
                        %                         ylim(limits);
                        %                         line([0 0], limits, 'Color', 'k', 'LineWidth', 2)
                    end
                    
                    
                    figure(fig.MissesHandle{indFig})
                    fig.MissesInfo{indFig}(1, 1).select(); hold on;
                    currAx = gca;
                    currAx.YAxis.Visible = 'off';
                    ylim(limits);
                    line([0 0], limits, 'Color', 'k', 'LineWidth', 2)
                end
            end
            
            %save(fullfile(obj.outputFolder, ['meanTarget' plotName '.mat']), 'meanTarget')
        end
        
        
        function figInfo = plotAverageSurface(obj, Target, plotName, grouping, plotComb)
            %% %%%%%%%%%%%% surface plots %%%%%%%%%%%%%%%%%%%%%%%%%
            % Surface plot CPP sorted along the y-axis by RT
            % z-scored inside participant, interchange interval duration,
            % and contrast change direction to facilitate the pooling of data
            % across participants and exclude the influence of the aforementioned
            % experimental factors
            if ~exist('grouping', 'var'); grouping = 0; end
            
            % preset figure save folder
            currFolder = fullfile(obj.figFolder,  'groupAverage', ['HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)], [plotName '/']);
            if ~exist(currFolder, 'dir'); mkdir(currFolder); end
            keyboard
            % 1) first get 'remove between participant variablity'.
            % get the individual means and the grandAverage
            indivMeanTarget   = repmat(nanmean(Target,3),1,1,size(Target,3),1,1);
            grandMeanTarget   = repmat(nanmean(nanmean(Target,3),4),1,1,size(Target,3),size(Target,4));
            
            for indPP = 1:length(obj.ppNames)
                for indEpoch = 1:size(Target,3)
                    removeTarget(:,:,indEpoch, indPP)   = Target(:,:,indEpoch,indPP)   - indivMeanTarget(:,:,indEpoch,indPP)   + grandMeanTarget(:,:,indEpoch,indPP);
                end
            end
            
            % here seperated for up and down trails.
            % 1) reshape into table with rt per participant per condition
            withinVar = {}; withinName = {}; betweenVar = {}; betweenName = {};  withinPar = [];
            for indCond = plotComb
                if length(unique(obj.behaviour{1}.trialMatrix(:,indCond))) > 1
                    withinVar{end+1}  = obj.conditions{indCond};
                    withinName{end+1} = strrep(obj.figLayOut.legTitle{indCond},' ','');
                    withinPar(end+1) = 1;
                else
                    betweenVar{end+1}  = obj.conditions{indCond};
                    betweenName{end+1} = strrep(obj.figLayOut.legTitle{indCond},' ','');
                    withinPar(end+1) = 0;
                end
            end
            
            allRT       = []; allRTzScore = [];  allTarget   = []; allCond = [];
            
            counter = 0;
            for indPP = 1:length(obj.ppNames)
                goodTrials = ~isnan(obj.behaviour{indPP}.RT) & obj.behaviour{indPP}.RT > obj.stim.RTCutOff;
                
                if grouping
                    plotComb = plotComb(plotComb ~= 0);
                    VartrialMatrix  = obj.behaviour{indPP}.trialMatrix(:,plotComb(find(withinPar)));
                else
                    VartrialMatrix  = obj.behaviour{indPP}.trialMatrix;
                end
                
                posCombination = unique(VartrialMatrix(~isnan(obj.behaviour{indPP}.RT(:,1)),:), 'rows');
                
                for indCond = 1:length(posCombination)
                    currCond = sum(VartrialMatrix == posCombination(indCond,:),2) == size(posCombination, 2) & goodTrials;
                    
                    for indResponse = find(any(currCond) == 1)
                        currRT = obj.behaviour{indPP}.RT(currCond(:,indResponse));
                        
                        allCond(counter+1:counter+length(currRT),1:size(VartrialMatrix,2)) = VartrialMatrix(currCond(:,indResponse),:);
                        
                        if sum(~withinPar)
                            allCond(counter+1:counter+length(currRT),size(VartrialMatrix,2) + 1)  = unique(obj.behaviour{indPP}.trialMatrix(:,plotComb(find(~withinPar))));
                        end
                        
                        allRT(counter+1:counter+length(currRT),:)       = currRT;
                        allRTzScore(counter+1:counter+length(currRT),:) = zscore(currRT);
                        allTarget(:,counter+1:counter+length(currRT))   = squeeze(removeTarget(:,:, currCond(:,indResponse), indPP));
                        counter = counter + length(currRT);
                    end
                end
            end
            
            if grouping
                groups = unique(posCombination(:,grouping))';
            else
                groups = 1;
            end
            
            acc = 1;
            for indCond = groups
                if grouping
                    currCond = allCond(:,grouping) == indCond;
                else
                    currCond = 1:length(allCond);
                end
                
                currZscore  = allRTzScore(currCond);
                currRT      = allRT(currCond);
                currTarget  = allTarget(:,currCond);
                
                % find trials to be deletede due to artifacts or
                % no-data/blinks
                indBadEEG = all(isnan(currTarget));
                currZscore(indBadEEG) = [];
                currRT(indBadEEG) = [];
                currTarget(:,indBadEEG) = [];
                
                fig{acc} = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
                
                % 2.2) reshape rt into vector and sort and allRT to get line.
                [~, indRTs] = sort(currZscore);
                currRT      = sgolayfilt(currRT(indRTs),1,101);
                
                currTarget = currTarget(:,indRTs)';
                currTarget = sgolayfilt(currTarget,1,101);
                
                
                contourf(obj.eeg.targetEpoch,1:size(currTarget,1),currTarget,40,'linecolor','none');
                colormap('jet')
                
                if strcmp(plotName, 'CPP')
                    ylimAxis(acc,:)  = [0  max(nanmean(currTarget))];
                else
                    ylimAxis(acc,:) = round([min(nanmean(currTarget)) max(nanmean(currTarget))]);
                end
                hold on, plot(currRT, 1:length(currRT), 'k', 'LineWidth', obj.figLayOut.lineWidth);
                xlabel ('Time from target (s)')
                
                figInfo{acc} = gca;
                set(figInfo{acc} ,'FontSize', obj.figLayOut.letterSize);
                set(figInfo{acc} ,'FontName', obj.figLayOut.letterType);
                acc = acc+1;
            end
            
            ylimAxis = [min(ylimAxis(:)) max(ylimAxis(:))];
            xstuff = [fliplr(0:-0.5:obj.eeg.targetEpoch(1)) 0.5:0.5:obj.eeg.targetEpoch(end)];
            for indFig = 1:length(figInfo)
                figure(fig{indFig})
                
                if obj.eeg.applyCSD
                    limits = [floor(ylimAxis(:,1)/5)*5 ceil(ylimAxis(:,2)/5)*5];
                    if limits(1) == 0; limits(1) = -5; end
                    if limits(2) == 0; limits(2) = 5; end
                else
                    limits = [floor(ylimAxis(:,1)/1)*5 ceil(ylimAxis(:,1)/5)*5];
                    if limits(1) == 0; limits(1) = -1; end
                    if limits(2) == 0; limits(2) = 1; end
                end
                
                set(figInfo{indFig}, 'ydir','normal','xlim',[xstuff(1) obj.eeg.targetEpoch(end)], 'XTick', xstuff,...
                    'YTick', xstuff, 'clim', limits)
                tmp = plotComb(find(withinPar));
                
                if grouping
                    title(obj.figLayOut.legNames{tmp(grouping)}(indFig));
                    plotSave(figInfo{indFig}, [plotName 'SurfacePlot_' num2str(indFig) '.png'], currFolder,  [obj.figLayOut.saveDim(1)./2 obj.figLayOut.saveDim(2)]);
                else
                    plotSave(figInfo{indFig}, [plotName 'SurfacePlot_' num2str(indFig) '.png'], currFolder,  [obj.figLayOut.saveDim(1)./2 obj.figLayOut.saveDim(2)]);
                end
                
                
                cbar = colorbar;
                
                if obj.eeg.applyCSD
                    cbar.Ticks = limits(1):5:limits(2);
                else
                    cbar.Ticks = limits(1):1:limits(2);
                end
                
                if obj.eeg.applyCSD
                    cbar.Label.String = {[plotName ' (\muV/m^2)']};
                else
                    cbar.Label.String = {[plotName ' (\muV)']};
                end
                ylabel('Trials')
                
                plotSave(figInfo{indFig}, ['SurfacePlot_Legend_' num2str(indFig) '.png'], currFolder, [obj.figLayOut.saveDim]);
            end
        end
        
        function figInfo = plotAverageTopo(obj, topoData, plotName, channels, grouping)
            %% topoplot
            if nargin < 4
                channels = [];
            end
            
            currFolder = fullfile(obj.figFolder,  'groupAverage', ['HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)], [plotName{1} '/']);
            if ~exist(currFolder, 'dir'); mkdir(currFolder); end
            
            %% set up model parameters
            % as there are different parameters that you might wanna add,
            % like time dependency. For now it just simple time on task and
            % time within the experiment.
            if grouping ~= 0
                [within, between, ~, allTrialMatrix] = getConditions(obj, grouping);
            else
                allTrialMatrix = ones(max(obj.numBlocks.*obj.numTrials),1, length(obj.ppNames));
                between.Table = ones(length(obj.ppNames),1);
                within.numCond = 1;
            end
            
            plotTopo = nan(size(topoData,1), size(topoData,2), within.numCond, length(obj.ppNames));
            if size(topoData,1) > 1
                keyboard
                % TODO for left and right stuff.
            end
            
            for indPP = 1:length(obj.ppNames)
                
                posCombination = unique(allTrialMatrix(:,:, indPP), 'rows');
                posCombination(any(isnan(posCombination),2),:) = [];
                
                for indCond  = 1:size(posCombination,1)
                    currCond = sum(allTrialMatrix(:,:,indPP) == posCombination(indCond,:),2) == size(posCombination, 2);
                    
                    % remove all conditions that didn't have a
                    % response, e.g. the nans.
                    currTopo      = topoData(1,:,currCond,indPP);
                    
                    % get mean reaction times
                    plotTopo(:,:,indCond, indPP) = nanmean(currTopo,3);
                end
                
            end
            
            acc = 0;
            for indTopo = 1:size(plotTopo,1)
                
                currTopoData(1, :,size(plotTopo,3),:) = plotTopo(indTopo,:,:,:);
                for indGroup = unique(between.Table)'
                    for indCond = 1:size(plotTopo,3)
                        acc = acc+1;
                        
                        figure('units','normalized','outerposition',[0 0 1 1]); hold on;
                        currFig = panel();
                        currFig.pack(1, 1);
                        currFig(1,1).marginright = 4;
                        currFig(1,1).marginleft  = 4;
                        currFig.marginbottom = 4;
                        currFig.marginright  = 4;
                        currFig.fontsize = obj.figLayOut.letterSize;
                        currFig.fontname = obj.figLayOut.letterType;
                        plotFig{acc} = gcf;
                        
                        if isempty(channels)
                            topoFigInd{acc} = topoplot(squeeze(nanmean(currTopoData(1,:,indCond,between.Table == indGroup),4)),obj.eeg.chanlocs,...
                                'electrodes', 'labels');
                        else
                            topoFigInd{acc} = topoplot(squeeze(nanmean(currTopoData(1,:,indCond,between.Table == indGroup),4)),obj.eeg.chanlocs,...
                                'electrodes', 'labels');
                            
                            ChanPlot = {obj.eeg.ChannelsName{channels}};
                            keep = [];
                            for kk = 1:length(topoFigInd{acc}.Parent.Children)
                                for indChan = 1:length(ChanPlot)
                                    try
                                        if ~strcmpi(topoFigInd{acc}.Parent.Children(kk).String(~isspace(topoFigInd{acc}.Parent.Children(kk).String)), ChanPlot{indChan})
                                            keep(kk, indChan) = 0;
                                        else
                                            keep(kk, indChan) = 1;
                                        end
                                    end
                                end
                            end
                            
                            topoFigInd{acc}.Parent.Children( find(~sum(keep,2))').delete
                            for indChan = 1:length(find(sum(keep,2)))
                                topoFigInd{acc}.Parent.Children(indChan).String = '*';
                                topoFigInd{acc}.Parent.Children(indChan).FontSize = obj.figLayOut.letterSize;
                            end
                        end
                        getLim(:,acc) = topoFigInd{acc}.Parent.CLim;
                    end
                end
            end
            
            % set Clim
            for indPlot = 1:length(plotFig)
                figInfo{indPlot} = figure(plotFig{indPlot});
                if obj.eeg.applyCSD && any(getLim > 10)
                    limits = [floor(min(getLim(1,:))/5)*5 ceil(max(getLim(2,:))/5)*5];
                    if limits(1) == 0; limits(1) = -5; end
                    if limits(2) == 0; limits(2) = 5; end
                else
                    limits = [floor(min(getLim(1,:))/1)*1 ceil(max(getLim(2,:))/1)*1];
                    if limits(1) == 0; limits(1) = -1; end
                    if limits(2) == 0; limits(2) = 1; end
                end
                topoFigInd{indPlot}.Parent.CLim = limits;
                
                cb = colorbar('eastoutside');
                cb.Ticks = [limits(1) 0 limits(2)];
                
                if ~ischar(grouping) && grouping ~= 0
                    title(obj.figLayOut.legNames{grouping(1)}(indPlot));
                end
                
                plotSave(gcf, [strcat(plotName{:}) '_Topo_' num2str(indPlot) '.png'], currFolder, [obj.figLayOut.saveDim(1) obj.figLayOut.saveDim(1)]);
            end
        end
        
        function [figHandle, figInfo] = plotDifferenceTopo(obj, plotName, Elec, plotComb, trangeTopo, TargetOrResponse,...
                BaselineCorrect, MorletOrsFFT, freqRange, PeakMeanOrMax, plotpValues)
            %% topoplot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -------------- PRE-SET PARAMETERS  -------------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % check for parameters.
            if ~exist('TargetOrResponse', 'var');   TargetOrResponse = 2;  end  % standard on response-locked
            if ~exist('MorletOrsFFT', 'var');       MorletOrsFFT = 0;      end  % standard on average
            if ~exist('plotpValues', 'var');        plotpValues = 0;      end  % standard on average
            
            trangeBaseline = obj.eeg.epochPlot > obj.eeg.baseline(1) & obj.eeg.epochPlot < obj.eeg.baseline(2);
            
            % preset folders
            currOutput = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName);
            if ~exist(currOutput, 'dir'); mkdir(currOutput); end
            
            % preset .mat file folder
            currTopo = fullfile(currOutput, [plotName 'Topo_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD) '.mat']);
            
            if ischar(Elec)
                tmpChan = split(Elec); clear Elec;
                for indChan = 1:length(tmpChan)
                    Elec(indChan) = find(strcmp(tmpChan{indChan}, obj.eeg.ChannelsName));
                end
            elseif iscell(Elec)
                for indPP = 1:length(Elec)
                    tmpChan = split(Elec{indPP});
                    for indChan = 1:length(tmpChan)
                        indElec(indPP, indChan) = find(strcmp(tmpChan{indChan}, obj.eeg.ChannelsName));
                    end
                end
                clear tmp*
                Elec = unique(indElec)';
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -------------- CREATE AVERAGE TOPOPLOT  --------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for electrode selection, we first create an average topoplot.
            % This allows us to visualize and properly chooise the cluster
            % of electrodes the check.
            % pre-set filename
            if TargetOrResponse == 1
                trangeTopo = obj.eeg.targetEpoch > trangeTopo(1) & obj.eeg.targetEpoch < trangeTopo(2);
            elseif TargetOrResponse == 2
                trangeTopo = obj.eeg.responseEpoch > trangeTopo(1) & obj.eeg.responseEpoch < trangeTopo(2);
            end
            
            if ~exist(currTopo, 'file')
                fprintf('Extract data for avarage topoplot to determine best %s channels.\n', plotName)
                
                % 1) extract all topoplot per participant
                ERPTopo = nan(1,obj.eeg.NumberOfChannels, max(obj.numBlocks)*obj.numTrials , length(obj.ppNames));
                
                for indPP = 1:length(obj.ppNames)
                    clear ERP tmp*
                    fprintf(['Now processing participant ' obj.ppNames{indPP} ' to get topoplot\n'])
                    currInput  = fullfile(obj.outputFolder,'EEG data', obj.ppNames{indPP});
                    
                    % load the EEG data.
                    if obj.eeg.applyCSD
                        load(fullfile(currInput,[obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.stim.timing '.mat']), 'csdERP');
                        ERP = csdERP;
                    else
                        load(fullfile(currInput,[obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.stim.timing '.mat']), 'ERP');
                    end
                    
                    % Baseline-correct the data for the target and
                    % response. Again not for the ERPwhole
                    ERP = ERP - repmat(nanmean(ERP(:, trangeBaseline, :),2), 1, size(ERP,2), 1);
                    
                    % calculated ERP topography
                    if TargetOrResponse == 1
                        [~, tmpERP] = sortERPs(obj, ERP, indPP, 1, 1);
                    elseif TargetOrResponse == 2
                        [~, ~, tmpERP] = sortERPs(obj, ERP, indPP, 2, 1);
                    end
                    
                    tmpERP = tmpERP(:,trangeTopo,:,1);
                    
                    if MorletOrsFFT ~= 0
                        clear down*
                        % SSMVEP
                        
                        Fs = obj.eeg.SampleRate; % Sampling frequency
                        fftlen  = size(tmpERP,2); % samples
                        dF = Fs/fftlen;
                        
                        F = 0:dF:Fs-dF + (dF/2)*mod(fftlen,2);                                      % frequency in hertz
                        
                        intF = F  >= freqRange(1) & F <= freqRange(end);
                        
                        % remove SSVEP
                        if isfield(obj.stim, 'freqSSVEP')
                            intF(F  > obj.stim.freqSSVEP -1 & F  < obj.stim.freqSSVEP +1) = 0;
                        end
                        
                        fourTrans	  = abs(fft(tmpERP,[],2));%./(winLength/2);
                        clear tmpERP
                        tmpERP = fourTrans(:,  F >= freqRange(1) & F <= freqRange(end),:);
                        
                        if BaselineCorrect ~= 0
                            ERP = ERP - repmat(nanmean(ERP(:, trangeBaseline, :),2), 1, size(ERP,2), 1);
                        end
                    end
                    
                    if PeakMeanOrMax == 1
                        ERPTopo(1,:,:,indPP) = squeeze(nanmean(tmpERP,2));
                    elseif PeakMeanOrMax == 2
                        ERPTopo(1,:,:,indPP) = squeeze(max(tmpERP,[],2));
                    end
                    
                end
                
                % temporary dock this figure as to select electrodes
                currFig = gcf;
                set(currFig,'WindowStyle','docked');
                
                save(currTopo, 'ERPTopo');
                
            else
                load(currTopo, 'ERPTopo');
            end
            
            if nargin < 4
                channels = [];
            end
            
            currFolder = fullfile(obj.figFolder,  'groupAverage', ['HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD)], [plotName '/']);
            if ~exist(currFolder, 'dir'); mkdir(currFolder); end
            
            [within, between, ~, allTrialMatrix] = getConditions(obj,plotComb);
            
            for indPP = 1:length(obj.ppNames)
                posCombination = unique(allTrialMatrix(~isnan(obj.behaviour{indPP}.RT),:, indPP), 'rows');
                posCombination(any( isnan(posCombination),2),:) = [];
                
                for indCond = 1:size(posCombination,1)
                    currCond = all(allTrialMatrix(:,:, indPP) == posCombination(indCond,:), 2);
                    newTopo(indCond,:,indPP) = nanmean(ERPTopo(:,:, currCond, indPP),3);
                end
                
                newTopo(:,:,indPP) = reshape(nanzscore(reshape(newTopo(:,:,indPP),1,[])),size(newTopo,1),size(newTopo,2));
            end
            
            %% plot differences, first is going to be the party of interest,
            % eg. the plots that are shown, while the second is going to be
            % the subtraction!
            
            if length(plotComb) == 1
                posCombination = [ones(size(posCombination(:,1))) posCombination(:,1)];
            end
            
            AllSub = unique(posCombination(:,1));
            acc = 0;
            
            for indGroup = unique(between.Table)
                for indSub = 1:size(AllSub,1)
                    try
                        indCond = find(posCombination(:,1) == AllSub(indSub))';
                        Topo1 = squeeze(newTopo(indCond(1),:,:));
                        Topo2 = squeeze(newTopo(indCond(2),:,:));
                        
                        % diffTopo =  nanmean(Topo1,2) -  nanmean(Topo2,2)
                        
                        diffTopo(indSub,:,1:length(obj.ppNames)) = Topo1- Topo2;
                        diffTopo = squeeze(nanmean(diffTopo(indSub,:,:),3));
                        
                        if plotpValues
                            [test, ~, ~, stats] = ttest((Topo1-Topo2)');
                            stats.tstat(~logical(test)) = 0;
                            diffTopo = stats.tstat;
                        end
                        
                        acc = acc+1;
                        
                        h{acc} = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
                        
                        if isempty(Elec)
                            topoFigInd{acc} = topoplot(diffTopo,obj.eeg.chanlocs, 'electrodes', 'off');%, 'maplimits', 'minmax');
                        else
                            topoFigInd{acc} = topoplot(diffTopo,obj.eeg.chanlocs, 'electrodes', 'labels');
                            
                            ChanPlot = {obj.eeg.ChannelsName{Elec}};
                            keep = [];
                            for kk = 1:length(topoFigInd{acc}.Parent.Children)
                                for indChan = 1:length(ChanPlot)
                                    try
                                        if ~strcmpi(topoFigInd{acc}.Parent.Children(kk).String(~isspace(topoFigInd{acc}.Parent.Children(kk).String)), ChanPlot{indChan})
                                            keep(kk, indChan) = 0;
                                        else
                                            keep(kk, indChan) = 1;
                                        end
                                    end
                                end
                            end
                            
                            topoFigInd{acc}.Parent.Children( find(~sum(keep,2))').delete
                            for indChan = 1:length(find(sum(keep,2)))
                                topoFigInd{acc}.Parent.Children(indChan).String = '*';
                                topoFigInd{acc}.Parent.Children(indChan).FontSize = obj.figLayOut.letterSize;
                            end
                        end
                        
                        getLim(:,acc) = topoFigInd{acc}.Parent.CLim;
                        
                    catch
                        continue
                    end
                end
            end
            
            % set Clim
            acc = 0;
            for indGroup = unique(between.Table)
                figure(h{indGroup});
                for indSub = 1:size(AllSub,1)
                    
                    try
                        
                        acc = acc + 1;
                        figure(h{acc})
                        if obj.eeg.applyCSD && any(getLim(:) > 5)
                            
                            limits = [round(min(getLim(1,:))/5)*5 round(max(getLim(2,:))/5)*5];
                            if limits(1) == 0; limits(1) = -5; end
                            if limits(2) == 0; limits(2) = 5; end
                            
                        elseif obj.eeg.applyCSD && any(getLim(:) < 5)
                            
                            limits = [round(min(getLim(1,:))/0.1)*0.1 round(max(getLim(2,:))/0.1)*0.1];
                            if limits(1) == 0; limits(1) = -0.1; end
                            if limits(2) == 0; limits(2) = 1; end
                            
                            %                                                                     limits = [-0.6 0.6];
                        else
                            limits = [round(min(getLim(1,:))/1)*1 round(max(getLim(2,:))/0.1)*1];
                            if limits(1) == 0; limits(1) = -1; end
                            if limits(2) == 0; limits(2) = 1; end
                        end
                        
                        currPlot = gca;
                        currPlot.CLim = limits;
                        
                        cb = colorbar('eastoutside');
                        cb.Ticks = [limits(1) 0 limits(2)];
                        if length(plotComb) == 1
                            titlePlot = obj.figLayOut.legTitle{plotComb(1)};
                            saveName = sprintf('%s - %s', obj.figLayOut.legNames{plotComb(1)}{1},  obj.figLayOut.legNames{plotComb(1)}{2});
                        else
                            titlePlot = obj.figLayOut.legNames{plotComb(1)}{indSub};
                            saveName = sprintf('%s - %s', obj.figLayOut.legNames{plotComb(2)}{1},  obj.figLayOut.legNames{plotComb(2)}{2});
                        end
                        title(titlePlot);
                        plotSave(gcf, [plotName saveName titlePlot '.png'], currFolder, [obj.figLayOut.saveDim(1) obj.figLayOut.saveDim(1)]);
                    catch
                        keyboard
                    end
                end
            end
        end
        
        %% %%%%%%%%%%%% frequency analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % two different methods possible. Either short-time fourier
        % transformation or complex morlet wavelength
        
        % Power was measured in a sliding boxcar window of 240 ms with 20 ms step size.
        function [shortFFT, timeSeries] = shortfft(obj, ERP, intFreq)
            
            if isfield(obj.stim, 'freqSSVEP')
                windowSize = (4 * (1/obj.stim.freqSSVEP));
            else
                windowSize = 0.2;
            end
            
            Fs   = obj.eeg.SampleRate;       % Sampling frequency
            L    = floor(windowSize/(1/Fs)); % hertz per sample
            f    = (0:L-1)/L*Fs;
            intF = f >= intFreq(1) & f <= intFreq(end);
            
            if isfield(obj.stim, 'freqSSVEP')
                intF(f > obj.stim.freqSSVEP - 1 & f < obj.stim.freqSSVEP + 1) = 0;
            end
            
            stepSize    = ceil(0.05/(1/Fs)); % 100 ms steps
            windowWidth = floor(L/2);
            timeSeries  = windowWidth+1:stepSize:length(obj.eeg.epochPlot)-windowWidth;

            ind = 1;
            for tt = timeSeries
                fourTrans = abs(fft(ERP(:, tt-windowWidth:tt+windowWidth, :),[],2));%./(winLength/2);
                shortFFT(:,:, ind,:) = fourTrans(:, intF, :)/L;
                ind = ind+1;
            end
        end
        
        function output = morletWavelength(obj, ERP, frex, channels, numCycles)
            
            if ~exist('numCycles', 'var')
                numCycles = repmat(14, length(frex),1);
            elseif length(numCycles) == 1 && length(frex) > 1
                numCycles = repmat(numCycles,length(frex),1);
            end
            
            normalFreq   = 8:30; % mu and ERD range
            range_cycles = [15 22];
            
            % other wavelet parameters
            num_cycles  = logspace(log10(range_cycles(1)),log10(range_cycles(end)),length(normalFreq));
            [~,indFreq] = intersect(normalFreq, frex);
            numCycles  = num_cycles(indFreq);
            % numCycles = repmat(40, length(frex),1);
            
            if ~exist('channels', 'var') | isempty(channels)
                channels = 1:size(ERP,1);
            end
            
            % other wavelet parameters
            time = -2:1/obj.eeg.SampleRate:2;
            half_wave = (length(time)-1)/2;
            
            % FFT parameters
            nKern = length(time);
            
            try
                nData = sum(sum(~isnan(ERP(1,:,:))));
            catch
                keyboard
            end
            nConv = nKern+nData-1;
            
            % initialize output time-frequency data
            output = zeros(length(channels), length(frex), size(ERP,2), size(ERP,3));
            
            
            indAcc = 1;
            for indChan = channels
                % FFT of data (doesn't change on frequency iteration)
                % reshape ERP and remove nans
                catERP   = reshape(ERP(indChan,:,:),1,[]);
                tempNaNs = nan(size(catERP));
                
                indNaN = ~isnan(catERP);
                catERP(~indNaN) = [];
                
                dataX = fft(catERP,nConv);
                
                % loop over frequencies
                for fi=1:length(frex)
                    
                    % create wavelet and get its FFT
                    s = numCycles(fi)/(2*pi*frex(fi));
                    
                    cmw  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
                    cmwX = fft(cmw,nConv);
                    cmwX = cmwX./max(cmwX);
                    
                    % run convolution, trim edges, and reshape to 2D (time X trials)
                    as = ifft(cmwX.*dataX, nConv);
                    as = as(half_wave+1:end-half_wave);
                    
                    tempNaNs(indNaN) = as;
                    as = reshape(tempNaNs,size(ERP,2),size(ERP,3));
                    
                    % put power data into big matrix
                    output(indAcc, fi, :, :) = abs(as).^2;
                end
                
                indAcc = indAcc + 1;
            end
            
            if length(frex) == 1
                output = squeeze(output);
            end
        end
        
    end
    
    methods(Static)
        function ins = getInstance()
            persistent instance;
            
            if( ~strcmpi(class(instance), 'dataAnalysis') )
                instance = dataAnalysis();
            end
            
            ins = instance;
        end
    end
    
    methods(Access = private)
        function obj = dataAnalysis()
            % Pre-set all variables.
            obj.inputFolder     = {};
            obj.outputFolder    = {};
            obj.logFolder       = {};
            obj.figFolder       = {};
            obj.dataFiles       = {};
            
            obj.ppNames         = {};
            obj.conditions      = {};
            obj.condNames       = {};
            obj.numBlocks       = [];
            obj.numTrials       = [];
            obj.order           = {};
            
            
            obj.analysisEEG     = 0;
            obj.analysisEyelink	= 0;
            obj.analysisEOG     = 0;
            obj.analysisEMG     = 0;
            
            obj.stim            = [];
            obj.triggers        = [];
            obj.eeg             = [];
            obj.system          = {};
            obj.event           = {};
            obj.experiment      = {};
            obj.behaviour       = {};
            obj.modelBehaviour  =[];
            
            obj.DetectOrDisc   = [];
            obj.eyelink        = [];
        end
    end
end
