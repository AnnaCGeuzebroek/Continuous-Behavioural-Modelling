classdef dataModelling < handle
    
    properties
        % set up parameters
        inputFolder
        outputFolder
        logFolder
        figFolder
        
        ppNames
        conditions
        condNames
        numBlocks
        numTrials
        order
        
        stim
        modelBehaviour
        options
        behaviour
        DetectOrDisc
        
        % plotting information
        figLayOut
    end
    
    methods(Access = public)
        
        function [bestpmFits, EE] = applyModelling(obj, modelName, numRepeats,  modelPar, initialGuess, jitter, freeParams)
            %% applyModelling
            % set the intial parameters, for now it is just an extensive
            % search within the bounds given in modelPar. After which it
            % will run through a parfor to optimize the intial search
            % parameters.
            % Use as:
            %   modelName   =   Just in order to save stuff (such as
            %                   parameters fits).
            %   numRepeats  =   Number use to sample the parameters space.
            %                   otherwise just use the intialGuess to
            %                   refine.
            %   intialGuess =   when left empty 'generateRandParameters' is
            %                   used to search whole parameters space given
            %                   by the limits in modelPar.
            %                   Otherwise, it will take the intially
            %                   crudely fitted parameters to start of with.
            % There are also a couple of ways to use intial guesses to
            % jitter around or created completely new parameters.
            %	jitter      =	Sometimes a new estimated can be informed
            %                   by previously fitted parameters. This
            %                   can use intialGuess.
            %	freeParams  =   Create completely new parameters with
            %                   modelPar which will be added to 
            
            
            saveFolder = fullfile(obj.outputFolder, 'Modelling');
            if ~exist(saveFolder, 'dir'); mkdir(saveFolder); end
            
            if ~exist('initialGuess', 'var');   initialGuess   = []; end
            if ~exist('freeParams', 'var');     freeParams = []; end
            if ~exist('jitter', 'var');         jitter = []; end
            
            if isempty(initialGuess) || ~isempty(freeParams)
                % for initial modelling, we start with generate a random matrix with
                % starting parameter vectors the specific range given in
                % modelPar.
                if ~isempty(freeParams)
                   
                    index  = 1:size(modelPar,1);
                    
                    % preset newParams, with 'little' var around
                    % intialGuess that are 'fixed'.
                    newParams = nan(numRepeats, size(modelPar,1));
                    newParams(:, ~ismember(index, freeParams)) = initialGuess(randi(size(initialGuess,1),numRepeats,1),:) + (randn(numRepeats, size(initialGuess,2)).*mean(initialGuess).*0.1);
                    
                    for indRow = 1:size(freeParams,1)
                        if isempty(jitter) ||(indRow == 1 && size(freeParams,1) > 1)
                            newParams(:,ismember(index, freeParams(indRow,:)))  = generateRandParameters(obj, modelPar(freeParams(indRow,:),:), numRepeats);
                        else
                            newParams(:,ismember(index, freeParams(indRow,:)))  = jitter(randi(size(jitter,1),numRepeats,1),:) +  (randn(numRepeats, size(jitter,2)).*mean(jitter).*0.15);
                        end
                    end
                    initialGuess = newParams;
                else
                    initialGuess = generateRandParameters(obj, modelPar, numRepeats);
                end
                
                saveFile = ['paramfits_' modelName '_1'];
            else
                % when model is already run, we are going to reuse the
                % refitted initialGuesss and refine them. We here save the
                % files as refine.
                saveFile = ['paramfits_' modelName '_2'];
            end
            
            % include check if it already exist.
            if ~exist(fullfile(saveFolder, [saveFile '.mat']), 'file')
                EE  = nan(size(initialGuess,1),1);
                pmfits = nan(size(initialGuess));
                
                numInitialGuess = length(initialGuess);
                modelNames = modelPar.Names;
                lowerBound = modelPar.Lower';
                upperBound = modelPar.Upper';
                
                noise = obj.modelBehaviour.noiseSTD*randn(length(obj.modelBehaviour.datsum.TOW), max(obj.modelBehaviour.datsum.maxln), obj.modelBehaviour.simulateMoreX); 
               
                % move through the fit parameter space and 'optimize'
                % parameters with SIMPLEX algorithm. 
                parfor indPM = 1:numInitialGuess
                    fprintf('start pt %i:', indPM)
                    [pm, Err] = fminsearchbnd(@(pm) diffusionModel(obj, pm, modelNames, noise), initialGuess(indPM,:), lowerBound, upperBound,  obj.options);
                    
                    % save them:
                    EE(indPM)       = Err;
                    pmfits(indPM,:) = pm';
                    
                    % only plot if error semi-decent
                    if Err <= 100
                        fprintf('start pt %i: E = %0.2f, pm = %s\n', indPM, EE(indPM), num2str(pmfits(indPM,:)))
                    end
                end
                
                save(fullfile(saveFolder, saveFile), 'EE', 'pmfits', 'initialGuess')
            else
                load(fullfile(saveFolder, saveFile), 'EE', 'pmfits')
            end
            
            % Get best parameter fits as output
            [EE, BestEE] = sort(EE);
            bestpmFits   = pmfits(BestEE, :);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ---------------  modelling functions     ----------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %-----------------  Leaky accumulation  ---------------------------
        function [err, simdat, pred, DecValue] = diffusionModel(obj, pm, parNames, noise, getsimDV)
            %% diffusionModel
            % Implements the sequential sampling modelling based on
            % accumulation where several parameters can be set to be free.
            % Function can be used to optimize within SIMPLEX function (see
            % obj.applyModelling), as well as just on each own, allowing
            % you to extract a prediction of a certain parameters
            % combination and the simulated decision value.  Model based on 
            % Ossmy, et al 2013 (https://pubmed.ncbi.nlm.nih.gov/23684972/)
            % Parameters can be set free: 
            %   bound   =   Threshold determining the amount of accumulate
            %               before commiting to a decision.
            %   drift   =   average amount of evidence accumulated per unit
            %               time, this can be an index of task difficulty
            %               or on participants ability. When having several
            %               difficutly levels, 1) this can either be set as
            %               one parameters which will be scaled accordingly
            %               2) one for each difficulty level or 3)
            %               indiviudally fitted drift per difficulty level
            %               per condition. 
            %   leak    =   How much previous accumulated evidence (t-1) is
            %               weighted to current incoming evidence (t). Set
            %               between 0 = perfect accumulation and 1 = full
            %               leak e.g. decide based on current evidence
            %               only. 
            %   ndt     =   non-decision time accounting for motor response
            %               delays
            %   boost   =   possibilty boosting the drift rate, if we
            %               expect learning effects --> index to which
            %               condition in the first column of the
            %               trialmatrix.
            %   criteria =  reference on momentary evidence, determine what is see as
            %               signal or noise (similar to signal detection theory)
            %
            %
            % inputParameters include
            %   1) pm       =   which are either estimated in fminsearchbnd using
            %                   going through the fit parameters space or a
            %                   single set of parameters to get the output. 
            %   2) parNames =   Names corresponding to the pm, should be
            %                   ones as described above and for each pm
            %                   one! 
            %   3) noise    =   set noise levels, e.g. the stocastic part.
            %   4) getsimDV =   0 - don't, 1 - get simulated DV. 
            %   
            
            if ~exist('getsimDV','var'); getsimDV = 0; end
            
            TOW    = obj.modelBehaviour.datsum.TOW(1,:);
            bt     = obj.modelBehaviour.datsum.bt;
            trialMatrix = obj.modelBehaviour.datsum.trialMatrix;
            reflectingBound = obj.modelBehaviour.reflectingBound;

            % get timing parameters
            dt = 1/obj.stim.refreshRate;      % sample period
            maxRT = obj.stim.RTdeadLine(end); % the upper limit for RTs to be classed as hits
            minRT = obj.stim.RTCutOff;        % minium allowable RT to call a response a 'hit', in msec. There was only one RT in the S blocks across all subjects that was shorter than this, and next shortest was nearly 100 ms later
            
            
            % don't start accumulating until this number of samples after response. This is simply a reasonable guess.
            % HOWEVER, it might be that targets appear together and this needs to be adjusted!!
            % TODO make this independent on refresh rate and rather on
            % timing!!
            postRespPause = 1*obj.stim.refreshRate; 
            
            % Without loss of generality, and to facilitate DV simulation and comparison with NI models later, will assign 67 ms of
            % the nondecision time to the motor end, and the free part will be for 'PRE' decisional NDT:
            
            MT = 4*dt;
            continAccDur = obj.stim.refreshRate/10; % for how many samples (100 ms) should the CPP keep accumulating after reaching commitment? From R-locked CPPs it looks like around 100 ms

            targetPlot = floor((obj.stim.targetEpoch(1))/dt):ceil((obj.stim.targetEpoch(end))/dt);
            responsePlot = floor((obj.stim.responseEpoch(1))/dt):ceil((obj.stim.responseEpoch(end))/dt);

            
            % preset and extract parameters.                                         
            bound = []; drift = []; leak = []; ndt = []; boost = []; criteria = []; trialnoise = [];
            for indParName = 1:size(parNames,1)
                eval(sprintf('%s(end+1) = pm(:, %i);', lower(parNames(indParName,:)), indParName))
            end
            
            % Some parameters are not going to be used, they do need to be
            % defined, therefore to save time in the for loop, we set them
            % here.
            
            if isempty(bound),      bound = 1; end
            if isempty(leak),       leak = 0; end
            if isempty(criteria),   criteria = 0; end
            if isempty(boost),      boost = 1; end
            if isempty(trialnoise), trialnoise = 1; end

            % adjust drift accordingly.
            if size(drift,2) == 1
                % if only one drift rate is fitted we are going to assume that the
                % drift rate scale depending on the coherence levels.
                
                % however, this is how it is implemented now, making it more
                % transfable (for example with more coherence levels):
                drift = drift*(obj.modelBehaviour.cohs./min(obj.modelBehaviour.cohs)); % scale up the drift rate parameters for other conditions
            elseif size(drift,2) > length(obj.modelBehaviour.cohs)
                drift = reshape(drift, length(obj.modelBehaviour.cohs), []);
            elseif size(boost,2) > 1
                drift = repmat(drift, 2,1).* boost';
            end
            
            % preallocated parameters.
            simdat  = [];
            
            acc = 1;
            hit = 1; miss = 2; false_alarm = 3;     % codes for response types
            
            % set-up for-loop
            for indNoise = 1:size(noise,3) % for robustness, it can help to simulate more trials than there are in the real data. This is in line with 'several' participants
                
                % loop through block trials TODO now only continuous, might
                % need to add a trial-loop instead of block for discrete
                % experiments.
                for indBlock = 1:length(TOW)
                    
                    % get the appropiated fit parametrs per conditions.
                    % When there is only 1 it will be fixed accross
                    % conditions. If the user set several it will find the
                    % appropiated current pm for that condition. Currently,
                    % this is just per block (e.g. bt(indBlock, 1)), TODO
                    % switch to the condtion using getCondition ect. 
                    
                    if length(bound) > 1;      currBound = bound(bt(indBlock,1));         else; currBound = bound; end
                    if length(leak) > 1;       currLeak = leak(bt(indBlock,1));           else; currLeak = leak; end
                    if length(criteria) > 1;   currCriteria = criteria(bt(indBlock,1));	  else; currCriteria = criteria; end 
                    if length(trialnoise) > 1; currTrialNoise = trialnoise(bt(indBlock)); else; currTrialNoise = trialnoise; end

                    if obj.modelBehaviour.learnBoost == bt(indBlock,1)
                        currBoost = boost; 
                    else
                        currBoost = 1;
                    end

                    if size(drift,1) > 1
                        if obj.modelBehaviour.learnBoost == bt(indBlock,1)
                            currDrift = drift(2,:);
                        else
                            currDrift = drift(1,:);
                        end
                    else
                        currDrift = drift.*currBoost;
                    end
                    
                    % make sensory evidence waveform - the target-on pulse train plus
                    % Gaussian noise   changed this to make it more transferable, 
                    % e.g. get the timeline which should be indexing the
                    % difficulty levels to get the right drift rate
                    timeLine = TOW{indBlock}(:,1)';
                    timeLine(timeLine ~= 0) = currDrift(abs(timeLine(timeLine ~= 0)));
                    
                    if obj.DetectOrDisc
                        direction = TOW{indBlock}(:,1)';
                        timeLine(direction < 1) = timeLine(direction < 1)*-1;
                    end
                    
                    sensEv = timeLine + currTrialNoise.*noise(indBlock,1:length(TOW{indBlock}), indNoise);
                    
                    % initalization of parameters. 
                    DV(1:round(ndt/dt)) = 0;   % initialize DV to a single 0 for the first time points of the block, up until the pre-decision nondecision time
                    lastresp = -postRespPause; % we only accumulate if it's a certain time since last response. Initialise to this value so accumulation begins right away at start of block
                    
                    RespT   = []; % keep track of all responses times
                    targT   = []; % keep track of all targets onset times
                    
                    if obj.DetectOrDisc
                        Resp = [];  % keep track of all actual response
                    end
                    
                    DVendval = nan; 
                    
                    % now simulate the block by looping through all sample points:
                    for indTime = round(ndt/dt)+1:length(sensEv)
                        
                        if indTime <= lastresp+continAccDur || indTime > lastresp + postRespPause % don't accumulate unless it has been a sufficient time since last response
                            
                            % this gives us the 'reflecting bound'
                            if reflectingBound
                                DV(indTime) = max((1-currLeak)*DV(indTime-1) + sensEv(indTime-round(ndt/dt)) - currCriteria, 0);% There is the main model equation! (see e.g. Ossmy et al 2013)
                                % TODO REFLECTING BOUND FOR
                                % DISCRIMINATION!!!
                            else
                                DV(indTime) = (1-currLeak)*DV(indTime-1) + sensEv(indTime-round(ndt/dt)) - currCriteria;   % There is the main model equation! (see e.g. Ossmy et al 2013)
                            end
                        else
                            % As we are using the empiraical CPP as an
                            % indeicator of the evidence accumulation, we
                            % want to see how the decision variable, e.g.
                            % DV, encodes and compares with this. In the
                            % future we can also use DV and CPP to contrain
                            % our models. However, this requires us to not
                            % have a 'absolute off' from the DV. 
                            
                            % Continue accumulation for 
                            if isnan(DVendval) % I'm using this as a way to know when it is time to linearly ramp the DV down to zero
                                DVendval = DV(indTime-1);
                            end
                            
                            DV(indTime) = max(0, DVendval*(1-(indTime-(lastresp+continAccDur))/25)); % linear decrease to a floor of zero
                            if DV(indTime) == 0, DVendval = nan; end % when the CPP has reached back down to zero, turn off the linear ramp-down
                        end
                        
                        % Detect target transitions
                        if abs(TOW{indBlock}(indTime)) > abs(TOW{indBlock}(indTime-1))
                            targT = [targT indTime*dt];
                        end
                        
                        % Detect responses: This is when the DV crosses the
                        % bound. 
                        if indTime > lastresp+postRespPause 
                            if obj.DetectOrDisc
                                if abs(DV(indTime)) > currBound     % TODO for disc task < -currBound beside RespT check if correct or mistake.
                                    RespT    = [RespT indTime*dt + MT];    % log the response after adding the non-decision time
                                    lastresp = indTime;                 % and now this is the last response that happened, at sample n
                                   
                                    if sign(DV(indTime)) == sign(TOW{indBlock}(indTime))
                                        Resp = [Resp 1];
                                    else
                                        Resp = [Resp 0];
                                    end
                                end
                            else
                                if DV(indTime) > currBound 
                                    RespT    = [RespT indTime*dt+MT];   % log the response after adding the non-decision time
                                    lastresp = indTime;                 % and now this is the last response that happened, at sample n
                                end
                            end
                        end
                    end
                    
                    % Now make the output matrix that has all the RTs w.r.t. target onset
                    % and false alarms. This has to be equivalent to how the real data were
                    % analysed! Deals in seconds. FOR FUTHER, I want to
                    % make the simdat more readable, and add it in a
                    % structure more like behavior as in obj.behaviour!
                    
                    for indTarget = 1:length(targT)
                        nextrespind = find(RespT > targT(indTarget) + minRT & RespT < targT(indTarget) + maxRT, 1); % find the index of the next response which is within the allowable @hit@ window
                        if ~isempty(nextrespind) 
                            if obj.DetectOrDisc
                                simdat  = [simdat; trialMatrix{indBlock}(indTarget,:) hit RespT(nextrespind)-targT(indTarget) Resp(nextrespind)]; % 1 = hit
                            else
                                simdat  = [simdat; trialMatrix{indBlock}(indTarget,:) hit RespT(nextrespind)-targT(indTarget)]; % 1 = hit
                            end
                        else % if there WAS no next response, set the response parameters for this trial as 'not a number'
                            if obj.DetectOrDisc
                                simdat  = [simdat; trialMatrix{indBlock}(indTarget,:) miss nan nan]; % 2 = miss
                            else
                                simdat  = [simdat; trialMatrix{indBlock}(indTarget,:) miss nan];
                            end
                        end
                        
                        if getsimDV == 1
                            if (targT(indTarget)/dt) + ceil((obj.stim.targetEpoch(end))/dt) < length(DV)
                                DecValue.Target(1:length(targetPlot),acc) = DV(round(targT(indTarget)/dt) + targetPlot);
                                
                                if ~isempty(nextrespind)
                                    DecValue.Response(1:length(responsePlot),acc) = DV(round(RespT(nextrespind)/dt) + responsePlot);
                                else
                                    DecValue.Response(1:length(responsePlot),acc) = nan;
                                end
                                
                            else
                                DecValue.Target(1:length(targetPlot),acc) = nan;
                                if ~isempty(nextrespind)
                                    DecValue.Response(1:length(responsePlot),acc) = nan;
                                end
                            end
                        end
                        acc = acc + 1;
                    end
                    
                    % TODO check what to do with FA in disc task
                    ITIstartT = [0 targT+obj.stim.duration];
                    
                    for indITI = 1:length(ITIstartT)-1
                        nexttargind = find(targT > ITIstartT(indITI),1); % index of next target
                        % from this establish the end of the window starting from the current ITI start where we will check for FAs
                        if ~isempty(nexttargind)
                            endtime = targT(nexttargind) + minRT;
                        else
                            endtime = length(DV)*dt-.125; % like in real data, if we didn't find a next target then this must be the end of the block, so check up as far as we ould possibly extract an erpr
                        end
                        % now find responses in this ITI window:
                        nextrespind = find(RespT > ITIstartT(indITI)+maxRT-(obj.stim.duration) & RespT < endtime); % find indices of responses the ITI window, ruling out any at the very start that are within the hit window from the previous target. Target duration is 1 sec, so fs in sample points
                        
                        for m = 1:length(nextrespind)
                            if obj.DetectOrDisc
                                simdat  = [simdat; trialMatrix{indBlock}(indITI,:) false_alarm RespT(nextrespind(m))-endtime  Resp(nextrespind(m))];
                            else
                                simdat  = [simdat; trialMatrix{indBlock}(indITI,:) false_alarm RespT(nextrespind(m))-endtime];
                            end
                            
                            if getsimDV == 1
                                DecValue.Target(1:length(targetPlot),acc)   = nan;
                                
                                try
                                    DecValue.Response(1:length(responsePlot),acc) = DV(round(RespT(nextrespind(m))/dt) + responsePlot);
                                catch
                                    DecValue.Target(1:length(responsePlot),acc)   = nan;
                                end
                                
                            end
                            acc = acc + 1;
                        end
                        
                    end
                end
            end
            
            % It's possible that with certain parameters the DV NEVER crosses the bound, so need a quick fix so no error happens:
            % this happened once by total fluke - there was no 'E' below, so there must have been a false alarm before the first target (cond=nan) and nothing else!
            if isempty(simdat)
                tmpCond = unique(reshape([trialMatrix{:}]', size(trialMatrix{1},2), [])', 'rows');
                simdat = [tmpCond repmat(3, size(tmpCond,1),1) repmat(10, size(tmpCond,1),1)];
            end
            
            if obj.modelBehaviour.ChiOrG == 1
                [err,pred] = Chisquared(obj, simdat);
            elseif obj.modelBehaviour.ChiOrG == 2
                [err,pred] = Gsquared(obj,simdat);
            end
            
        end
        
        %----------------- swarm approach to the parameters ---------------
        function randStartPar = generateRandParameters(obj, modelPar, numRepeats)
            if istable(modelPar)
                modelPar = [modelPar.Lower modelPar.Upper]';
            else
                modelPar = modelPar';
            end
            randStartPar = rand(numRepeats, size(modelPar,2)).*diff(modelPar,[],1) + modelPar(1,:);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% --------------- plotting functions    --------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [pred, currPar] = plotModellingResults(obj, bestpm, modelName, parNames, err, simdat, DV)
           
            % preset saving folder. 
            figureFolder = fullfile(obj.figFolder,  'groupAverage/Modelling/', [modelName '/']);
            if ~exist(figureFolder, 'dir')
                mkdir(figureFolder);
            end
            
            % preset parameters to created figures. 
            dt  = 1/obj.stim.refreshRate;
            datsum = obj.modelBehaviour.datsum;

            targetEpoch   = obj.stim.targetEpoch(1):dt:obj.stim.targetEpoch(end)+dt;
            responseEpoch = obj.stim.responseEpoch(1):dt:obj.stim.responseEpoch(end)+(dt*2);
            trangeBaseline = targetEpoch >= obj.stim.baselineCorrect(1) & targetEpoch <= obj.stim.baselineCorrect(2);

            % Get model output.
            % for NI models these are given as an input currently.
            % Otherwise, fit them here. 
            if ~(exist('err', 'var') && exist('simdat', 'var'))         
                noise = obj.modelBehaviour.noiseSTD*randn(length(obj.modelBehaviour.datsum.TOW), max(obj.modelBehaviour.datsum.maxln), obj.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
                [err, simdat, ~, DV] = diffusionModel(obj, bestpm, parNames, noise, 1);
            end
            
            % get goodness-of-fit for model comparison.
            [AIC, BIC] = obj.penalize(parNames,err);
            
            currPar = table(parNames, bestpm', 'VariableNames', {'Parameter', 'Value'});
            addGOFString = char(['err' repmat(' ',  1, length(currPar.Parameter(1,:))-3)],...
                ['AIC' repmat(' ',  1, length(currPar.Parameter(1,:))-3)],...
                ['BIC' repmat(' ',  1, length(currPar.Parameter(1,:))-3)]);
            currPar = [currPar; table(addGOFString, [err, AIC, BIC]', 'VariableNames', {'Parameter', 'Value'})];
            
            
            %% ---------- extract behaviour ------------------------------
            if obj.DetectOrDisc
                cond       = simdat(:,1:end-3);
                outcome    = simdat(:,end-2);
                RT         = simdat(:,end-1);
                Response   = simdat(:,end);
                Response(isnan(Response)) = 0;
            else
                cond       = simdat(:,1:end-2);
                outcome    = simdat(:,end-1);
                RT         = simdat(:,end);
            end
            
            uniCond = unique(cond(all(cond ~= 0,2),:), 'rows');
            qps     = obj.modelBehaviour.qps;
            numBins = length(qps) + 1; furtherIndex = numBins+1;

            % preallocated sim. decision value.
            meanTarget    = nan(size(DV.Target, 1),   length(uniCond));
            meanResponse  = nan(size(DV.Response, 1), length(uniCond));
                    
            for indCond = 1:size(uniCond,1)
                    
                currCond    = all(cond == uniCond(indCond,:), 2);
                currOutcome = outcome(currCond);    
                currMisses  = currOutcome == 2;
                currFA      = currOutcome == 3;
                
                currFART    = RT(currCond & outcome == 3);
                currRT      = RT(currCond & outcome == 1);
                
                currNumTrials = sum(currOutcome ~= 3);
                
                % first create simulated DV.
                meanTarget(:,indCond)   = nanmean(DV.Target(:, currCond & outcome == 1),2);
                meanResponse(:,indCond) = nanmean(DV.Response(:, currCond & outcome == 1),2);
             
                if obj.DetectOrDisc  % get only correct responses for the cummeliatve response function.
                    currCorrect = all(cond == uniCond(indCond,:), 2)  & Response;
                    currRT      = RT(currCorrect);
                    
                    currError   = all(cond == uniCond(indCond,:), 2)  & ~Response;
                    currRTError = RT(currError);
                end
                
                % 1) we'll get the average RT quantiles and numbers of trials in
                % each quantile, including FAs, because this is the behavioural
                % summary we wll try to fit. This is for Chi^2
                % e.g. TODO which optimiziation function
                q(indCond,:) = [0 quantile(currRT,qps)-0.5/60 max(currRT)];  % subtracted half of a refresh cycle just because response can only happen on each refresh, and setting the boundaries between quantile bins is then safer so > and >= don't give different results
                    
                for indBin = 1:numBins
                    % number of trials in each of the 6 quantile bins
                    qn(indCond,indBin) = sum(currRT > q(indCond,indBin) & currRT <= q(indCond,indBin+1));
                    % median RT within each quantile bin
                    pred.qm(indCond,indBin) = median(currRT(currRT > q(indCond,indBin) & currRT <= q(indCond,indBin+1)));
                    
                    % 2) Now convert the counts data into proportions so we can compute G^2
                    % conveniently, ALL subjects did exactly the sae number of targets, 288 in
                    % total, 72 for each of the four main conditions. So we can get proportions
                    % from the counts we derived above. Let's count the false alarms as the
                    % proportion of all 2-sec ITI intervals that contained a button click. The
                    % ITIs were randomly 2,4,6,8 sec equally likely.
                    % We'll use 'qn' from above, to get the per-subject counts
                    pred.pij(indCond,indBin) = qn(indCond,indBin)/(currNumTrials);
                end
                                                
                if obj.DetectOrDisc
                    keyboard % TODO GET THE THREE BINS?
                    for indQ = 1:3
                        % number of trials in each of the 3 bins
                        pred.pij(indCond,numBins+indQ)  =  sum(currRTError>datsum.q(indCond,numBins+1+indQ) & currRTError<=datsum.q(indCond,numBins+1+indQ+1))/currNumTrials;
                    end
                    if numBins ==  furtherIndex
                        furtherIndex = size(pred.pij,2)+1;
                    end
                end

                pred.pij(indCond,furtherIndex) = sum(currMisses)/currNumTrials;

                % then add a final bin that counts the false alarms:
                pred.pij(indCond,furtherIndex+1) = sum(currFA)/(datsum.nShortestITI(indCond)*obj.modelBehaviour.simulateMoreX);

                if ~isempty(obj.stim.FACutOff)
                    currFA = currRT(currOutcome == 3);
                    currFA = currFA > obj.stim.FACutOff;
                    pred(indCond,furtherIndex+2) = sum(currFA)/currNumTrials;
                end

                % note the min(1, .. limiter makes sure the proportion doesn't go above 1, which causes complex G-squared
                % values. I *think* this also might mean the starting vectors that are really far off in that they produce way
                % too many false alarms, are given up on earlier in the fminsearch procedure, which is a good thing - no point
                % in wasting time on them.
            end
            
            %% Quantile probability of the proportion plots
            % TODO on the moment only plotting for G squared!!
            figure; hold on;
            
            for indCond = 1:length(uniCond)
                
                % first plot prediction
                legendThis = plot(pred.qm(indCond,:),cumsum(pred.pij(indCond,1:end-2)),...
                    'LineStyle', obj.figLayOut.lineType{1},...
                    'Color',obj.figLayOut.colours(indCond,:),...
                    'MarkerSize', 3, 'LineWidth', 1);
                
                % Additionally add the empirically extracted behavioural
                % results. 
                CI95 = 1.96.*(datsum.pijstd(indCond,1:end-2)./sqrt(length(obj.ppNames)));
                
                h = errorbar(datsum.qm(indCond,:), cumsum(datsum.pij(indCond,1:end-2)), CI95, 'CapSize', 0);
                h.LineWidth = 1;
                h.Marker = 'o'; h.MarkerSize = 3;
                h.LineStyle = 'none';
                h.Color = [156 156 156]/255;
                
                h.MarkerEdgeColor = obj.figLayOut.colours(indCond,:);
                h.MarkerFaceColor = [1 1 1];
                legendThis(indCond) = h;
            end
            
            figInfo = gca;
            ylim([0 1.01])
            ylabel('Proportion')
            
            line([0 0],  [0 1.01], 'Color', 'k', 'LineWidth', 1.5)
            
            if obj.figLayOut.saveDim == [4 6]
                xlim([0 1.25])
                set(gca,'XTick', 0:0.5:obj.stim.RTdeadLine(end))
            else
                xlim([0.1 1.25])
                set(gca,'XTick', 0.2:0.2:obj.stim.RTdeadLine(end))
            end
            xlabel('Reaction times (sec)');

            set(gca,'FontSize', obj.figLayOut.letterSize);
            set(gca,'FontName', obj.figLayOut.letterType);
           
            plotSave(gcf, ['QuantileProb' modelName '.png'], figureFolder,  obj.figLayOut.saveDim);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -----------------     misses -------------------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure, hold on
            for indCond = 1:length(uniCond)
                CI95 = 1.96.*(datsum.pijstd(indCond,furtherIndex)./sqrt(length(obj.ppNames)));
                
                bb = bar(indCond, pred.pij(indCond,furtherIndex));
                set(bb,'FaceColor', obj.figLayOut.colours(indCond,:), 'FaceAlpha', 0.3)
                
                plot(indCond, datsum.pij(indCond,furtherIndex),'o',...
                    'MarkerEdgeColor', [1 1 1],...
                    'MarkerFaceColor', obj.figLayOut.colours(indCond,:),...
                    'LineWidth', obj.figLayOut.lineWidth)
               errorbar(indCond, datsum.pij(indCond,furtherIndex), CI95,'k','LineWidth',obj.figLayOut.lineWidth)

            end
            
            ylim([0 figInfo.YLim(end)])
            ylabel({'Proportion' 'misses'})
            title(sprintf('%s (%0.1f)', modelName, err))
            
            set(gca,'FontSize', obj.figLayOut.letterSize);
            set(gca,'FontName', obj.figLayOut.letterType);
            
            % proportional to the whole figure, e.g.
            plotSave(gcf, ['Misses' modelName '.png'], figureFolder,  [obj.figLayOut.saveDim(1) obj.figLayOut.saveDim(1)]);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -----------------     False Alarms --------------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure, hold on
            for indCond = 1:length(uniCond)
                CI95 = 1.96.*(datsum.pijstd(indCond,furtherIndex+1)./sqrt(length(obj.ppNames)));
                
                bb = bar(indCond, pred.pij(indCond,furtherIndex+1));
                set(bb,'FaceColor', obj.figLayOut.colours(indCond,:), 'FaceAlpha', 0.3)
                
                plot(indCond, datsum.pij(indCond,furtherIndex+1),'o',...
                    'MarkerEdgeColor', [1 1 1],...
                    'MarkerFaceColor', obj.figLayOut.colours(indCond,:),...
                    'LineWidth', obj.figLayOut.lineWidth)
                    
                errorbar(indCond, datsum.pij(indCond,furtherIndex+1), CI95,'k','LineWidth',obj.figLayOut.lineWidth)
            end
            
            ylim([0 figInfo.YLim(end)])
            ylabel({'Proportion' 'False Alarms'})
            title(sprintf('%s (%0.1f)', modelName, err))
            
            set(gca,'FontSize', obj.figLayOut.letterSize);
            set(gca,'FontName', obj.figLayOut.letterType);
            
            % proportional to the whole figure, e.g.
            plotSave(gcf, ['FA' modelName '.png'], figureFolder,  [obj.figLayOut.saveDim(1) obj.figLayOut.saveDim(1)]);
          
            %% CPP plotting
            figHandle   = figure('units','normalized','outerposition',[0 0 1 1]); hold on
            figResponse = figure('units','normalized','outerposition',[0 0 1 1]); hold on

            acc = 1;
            for indCond = 1:size(meanTarget,2)
                % this allows us to check when it is signifanctly
                % different from baseline, can be use in the plots
                % but right now not
                plotTarget = meanTarget(:,indCond);
                plotTarget(plotTarget == 0,:) = NaN;
                
                baseline = nanmean(plotTarget(trangeBaseline,:),1);

                plotResponse = meanResponse(:,indCond);
                plotResponse(plotResponse == 0,:) = NaN;
                
                plotTarget   = plotTarget   - baseline;
                plotResponse = plotResponse - baseline;
                
                YStuff(acc,1) = min([min(plotTarget) min(plotResponse)]); 
                YStuff(acc,2) = max([max(plotTarget) max(plotResponse)]); 
                
                %% plot Target-locked epochs
                figure(figHandle)
                
                % plot evidence onset
                h = plot(targetEpoch, plotTarget);
                h.Color = obj.figLayOut.colours(acc,:);
                h.LineWidth = obj.figLayOut.lineWidth;
                h.LineStyle = obj.figLayOut.lineType{acc};
                
                legendThis(acc) = h;
                
                %% plot response-locked epochs
                figure(figResponse)
                % plot evidence onset
                h = plot(responseEpoch(responseEpoch < 3*dt*-1), plotResponse(responseEpoch < 3*dt*-1));
                h.Color = obj.figLayOut.colours(acc,:);
                h.LineWidth = obj.figLayOut.lineWidth;
                h.LineStyle = obj.figLayOut.lineType{acc};
                
                acc = acc + 1;
            end
            limits = [min(YStuff(:,1)) max(YStuff(:,2)).*1.02];
            
            % layouting Target
            figure(figHandle)
            ylim(limits);
            line([0 0], limits, 'Color', 'k', 'LineWidth', 2)
         
            xlim([0 1])
            xticks(0:0.5:1)
            
            grid on
            xstuff = [fliplr(0:-0.5:obj.stim.targetEpoch(1)) 0.5:0.5:obj.stim.targetEpoch(end)];
            xlim([obj.stim.targetEpoch(1) obj.stim.targetEpoch(end)])
            xticks(xstuff)
            xlabel ('Time from target (s)')
            
            set(gca,'FontSize', obj.figLayOut.letterSize);
            set(gca,'FontName', obj.figLayOut.letterType);
     
            plotSave(gcf, ['simCPP' modelName '.png'], figureFolder,  [4 6]);
                   
            % layouting Response
            figure(figResponse)
            currAx = gca;
            ylim(limits);
            line([0 0], limits, 'Color', 'k', 'LineWidth', 2)
            xlabel(' ')
            
            if ~isempty(obj.figLayOut.responseLim)
                xstuffResponse = obj.figLayOut.responseLim;
            else
                xstuffResponse = [-0.5 0 0.2];
            end
            
            xlim([xstuffResponse(1) 0])
            xticks(xstuffResponse)
            
            yticklabels('')
            set(gca,'FontSize', obj.figLayOut.letterSize);
            set(gca,'FontName', obj.figLayOut.letterType);
            plotSave(gcf, ['simCPP' modelName 'Response.png'], figureFolder,  [4 2.6]);
            
            %% plot first derivative
            acc = 1;
            figHandle = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
           
            % get first derivative and add 
            diffmeanTarget = [diff(meanTarget); nan(1,size(meanTarget,2))];
            
            % normalize
            diffmeanTarget = (diffmeanTarget -  min(diffmeanTarget(:)))/(max(diffmeanTarget(:)) - min(diffmeanTarget(:)));
            
            for indCond = 1:size(meanTarget,2)
                
                % this allows us to check when it is signifanctly
                % different from baseline, can be use in the plots
                % but right now not    
                diffmeanTarget(:,indCond) = diffmeanTarget(:,indCond) -...
                    nanmean(diffmeanTarget(targetEpoch > -0.05 & targetEpoch < 0.05,indCond),1);
                plotTarget = diffmeanTarget(:,indCond);

                %% plot Target-locked epochs
                h = plot(targetEpoch, plotTarget);
                h.Color = obj.figLayOut.colours(acc,:);
                h.LineWidth = obj.figLayOut.lineWidth;
                h.LineStyle = obj.figLayOut.lineType{acc};
                
                legendThis(acc) = h;
               
                acc = acc + 1;
            end
            
            grid on
            xstuff = [fliplr(0:-0.5:obj.stim.targetEpoch(1)) 0.5:0.5:obj.stim.targetEpoch(end)];
            xlim([obj.stim.targetEpoch(1) obj.stim.targetEpoch(end)])
            xticks(xstuff)
            xlabel ('Time from target (s)')
            
            % layouting
            limits = [-0.5 1];
            ylim(limits);
            yticks(-0.5:0.5:1)
            line([0 0], limits, 'Color', 'k', 'LineWidth', 2)
            
            xlim([0 1])
            xticks(0:0.5:1)
            
            set(gca,'FontSize', obj.figLayOut.letterSize);
            set(gca,'FontName', obj.figLayOut.letterType);
            
            plotSave(gcf, ['CPP' modelName 'firstDeriv.png'], figureFolder,  [2.5 6]);
        end
        
        function plotParameters(obj, pm, modelName, parNames, parLimits)
            
            figureFolder = fullfile(obj.figFolder,  'groupAverage/Modelling/');
            if ~exist(figureFolder, 'dir')
                mkdir(figureFolder);
            end
            
            numPlot = size(pm,1);
            figure, hold on
            for indPlot = 2:size(pm,2)
                subplot(1,size(pm,2)-1,indPlot-1)
                scatter(pm(:,1), pm(:,indPlot), [], 1:numPlot);
                xlabel({parNames{1}})
                ylabel({parNames{indPlot}})
                axis square
            end
            linkaxes
            
            colormap( flipud(jet(numPlot)))
            
            axis square
            %{
            c = colorbar('Ticks', [1 40], 'TickLabels', {'Low', 'High'});
            if  obj.modelBehaviour.ChiOrG == 1
                c.Label.String = 'Chi^{2}';
            elseif  obj.modelBehaviour.ChiOrG == 2
                c.Label.String = 'G^{2}';
            end
            %}
            
            set(gca,'FontSize', obj.figLayOut.letterSize);
            set(gca,'FontName', obj.figLayOut.letterType);
            
            plotSave(gcf, ['currParameters' modelName '.png'], figureFolder,  [obj.figLayOut.saveDim(1) obj.figLayOut.saveDim(1)].*2);
            
        end
        
        function datsum = plotRTquantiles(obj, plotComb)
            %% plot the RT quantiles and get output for the model.
            % Now we'll get the average RT quantiles and numbers of trials
            % in each quantile, including FAs, because this is the behavioural
            % summary we wll try to fit
            
            
            figureFolder = fullfile(obj.figFolder,  'groupAverage/Modelling/');
            if ~exist(figureFolder, 'dir')
                mkdir(figureFolder);
            end
            
            %% set up model parameters
            % as there are different parameters that you might wanna add,
            % like time dependency. For now it just simple time on task and
            % time within the experiment.
            
            [~, ~, ~, allTrialMatrix] = getConditions(obj, plotComb);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                       Plot behavioural data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This is for illustration and get a quick figure with all results together.
            % This is not nescecary as the signalProcessing code also provides this,
            % but just quickly get everything together here (TODO cancel the plotting
            % by giving an input!)
            %     function extractBehavioural()
            
            qps = obj.modelBehaviour.qps; acc = 1;
            numBins = length(qps) + 1; furtherIndex = numBins+1; 
            for indPP = 1:length(obj.ppNames)
                
                % pool the data across all subjects as if it is one person, and plot histograms:
                currTrials = allTrialMatrix(:,:, indPP);
                uniCond = unique(currTrials, 'rows');
                uniCond(isnan(uniCond(:,1)),:) = [];
                
                obj.behaviour{indPP}.indFalseAlarm(obj.behaviour{indPP}.FalseAlarm == 0) = NaN;

                for indCond = 1:length(uniCond)
                    numberOfTrials(indCond, indPP) = sum(all(currTrials == uniCond(indCond,:), 2));
                    nShortestITI(indCond,indPP)    = sum(obj.behaviour{indPP}.ITI(all(currTrials == uniCond(indCond,:), 2)))./min(obj.behaviour{indPP}.ITI);

                    currCond   = find(all(currTrials == uniCond(indCond,:), 2));
                    
                    currMisses = obj.behaviour{indPP}.Misses(currCond);
                    currFA     = obj.behaviour{indPP}.FalseAlarm(currCond,:);
                    currRT     = obj.behaviour{indPP}.RT(currCond);

                    if obj.DetectOrDisc  % get only correct responses for the cummeliatve response function. 
                        currCorrect = all(currTrials(1:length(obj.behaviour{indPP}.Correct),:) == uniCond(indCond,:), 2)  & obj.behaviour{indPP}.Correct;
                        currRT      = obj.behaviour{indPP}.RT(currCorrect);
                        
                        currError   = all(currTrials(1:length(obj.behaviour{indPP}.Correct),:) == uniCond(indCond,:), 2)  & ~obj.behaviour{indPP}.Correct;
                        currRTError = obj.behaviour{indPP}.RT(currError);
                    end
                    
                    % 1) we'll get the average RT quantiles and numbers of trials in
                    % each quantile, including FAs, because this is the behavioural
                    % summary we wll try to fit. This is for Chi^2
                    % e.g. TODO which optimiziation function
                    q(indCond,1:numBins+1,indPP) = [0 quantile(currRT,qps)-0.5/60 max(currRT)];  % subtracted half of a refresh cycle just because response can only happen on each refresh, and setting the boundaries between quantile bins is then safer so > and >= don't give different results
                    
                    for indQ = 1:numBins
                        % number of trials in each of the 6 quantile bins
                        qn(indCond,indQ,indPP) = sum(currRT > q(indCond,indQ,indPP) & currRT <= q(indCond,indQ+1,indPP));
                        % median RT within each quantile bin
                        qm(indCond,indQ,indPP) = median(currRT(currRT > q(indCond,indQ,indPP) & currRT <= q(indCond,indQ+1,indPP)));
                    end
                    
                    if obj.DetectOrDisc
                        % 20% bins for errors to the conditional accuracy functions (CAFs)
                        % of each individual data set.
                        % CAFs represent accuracy as a function of time. In
                        % the Simon task, the incompatible condition is typically associated
                        % with an early drop of accuracy, resulting in a
                        % concave CAF shape (see Figure 5A, top). This particular
                        % pattern implies that incorrect responses are faster than
                        % correct responses. CAFs were constructed by sorting
                        % the RT data into five bins of equal size; the proportion
                        % of errors in each RT bin was then computed, providing
                        % the error data considered in the fitting procedure
                        q(indCond,numBins+2:numBins+2+3, indPP) = [0 quantile(currRTError, [1/3 2/3])-0.5/60  max(currRTError)];
                        
                        for indQ = 1:3
                            % number of trials in each of the 5 bins
                            qn(indCond,numBins+indQ,indPP)  = sum(currRTError > q(indCond,numBins+1+indQ,indPP) & currRTError <= q(indCond,numBins+1+indQ+1,indPP));
                            
                            % median RT within each quantile bin
                            qm(indCond,numBins+indQ,indPP)  = median(currRTError(currRTError > q(indCond,numBins+1+indQ,indPP) & currRTError <= q(indCond,numBins+1+indQ+1,indPP)));
                        end
                        
                        if numBins ==  furtherIndex
                            furtherIndex = length(qn)+1;
                        end
                    end
                    
                    % add a  bin that counts the Misses:
                    qn(indCond,furtherIndex,indPP)   = sum(isnan(currRT));
                    
                    % add a final bin that counts the false alarms:
                    qn(indCond,furtherIndex+1,indPP) = sum(currFA(:));
                       
                    if ~isempty(obj.stim.FACutOff)                  
                        limitedFalseAlarm = obj.behaviour{indPP}.FalseAlarm;                 
                        limitedFalseAlarm(obj.behaviour{indPP}.indFalseAlarm < obj.stim.FACutOff) = 0;      
                        qn(indCond,furtherIndex+2,indPP) = sum(sum(limitedFalseAlarm(currCond,:)));
                    end

                    % 2) Now convert the counts data into proportions so we can compute G^2
                    % conveniently, ALL subjects did exactly the sae number of targets, 288 in
                    % total, 72 for each of the four main conditions. So we can get proportions
                    % from the counts we derived above. Let's count the false alarms as the
                    % proportion of all 2-sec ITI intervals that contained a button click. The
                    % ITIs were randomly 2,4,6,8 sec equally likely.
                    % We'll use 'qn' from above, to get the per-subject counts
                    for indQ = 1:size(qn,2) % for the target quantiles (leaving out the last column which is false alarms)
                        % qn(c,i,s) is the number of trials in each of the 6 quantile bins
                        % Get proportions:
                        pij(indCond,indQ,indPP) = qn(indCond,indQ,indPP)/numberOfTrials(indCond, indPP);
                    end
                    
                    % ALSO record the total number of min ITI "trials" in total across all
                    % blocks of each type (so we can get the correct proportion of false alarms in the SIMULATED DATA)
                    pij(indCond, furtherIndex+1, indPP) = qn(indCond,furtherIndex+1,indPP)/nShortestITI(indCond,indPP);
                    
                    % ------------ addition of evidence timeline ----------
                    tmpMomEvidence = reshape(obj.behaviour{indPP}.MomEvidence(currCond,:,1)', [],1);
                    
                    if size(obj.behaviour{indPP}.MomEvidence,3) == 2
                        tmpMomEvidence2 = reshape(obj.behaviour{indPP}.MomEvidence(currCond,:,2)', [],1);
                        tmpMomEvidence = [tmpMomEvidence tmpMomEvidence2];
                    end
                    
                    tmpMomEvidence(isnan(tmpMomEvidence(:,1)),:) = [];
                    
                    datsum.TOW{acc}     = tmpMomEvidence;
                    datsum.maxln(acc)   = length(tmpMomEvidence);
                    datsum.bt(acc,:)    = uniCond(indCond,:);
                    datsum.trialMatrix{acc}  = currTrials(all(currTrials == uniCond(indCond,:), 2),:);
                    acc = acc + 1;
                    
                end
            end
            
            % Average across subjects and put in a structure:
            datsum.q  = nanmean(q,3);  % so q are the RT values separating the bins
            datsum.qn = nanmean(qn,3); % and qn are the average number of trials in each bin
            datsum.qm = nanmean(qm,3);
            datsum.numberOfTrials  = nanmean(numberOfTrials,2);

            datsum.nShortestITI    = sum(nShortestITI,2);

            datsum.pij = datsum.qn./datsum.numberOfTrials;
            datsum.pij(:,furtherIndex+1) = datsum.qn (:,furtherIndex+1)./nanmean(nShortestITI,2);
       
            % Average across subjects and put in the structure:
            %{   
            datsum.pij = nanmean(pij,3);

            if obj.DetectOrDisc
                datsum.qCAF  = nanmean(qCAF,3); % so q are the RT values separating the bins
                datsum.qnCAF  = nanmean(qnCAF,3); % so q are the RT values separating the bins
                datsum.pCAF  = nanmean(pCAF,3); % so q are the RT values separating the bins
                datsum.mCAF  = nanmean(mCAF,3);  % and qn are the average number of trials in each bin
            end
            %}
           
            % standard deviation.
            for indPP = 1:length(obj.ppNames)
                qstd(:,:,indPP)     = (q(:,:,indPP)  - nanmean(reshape(q(:,:,indPP),1,[])));
                qnstd(:,:,indPP)    = qn(:,:,indPP)  - nanmean(reshape(qn(:,:,indPP),1,[]));
                qmstd(:,:,indPP)    = qm(:,:,indPP)  - nanmean(reshape(qm(:,:,indPP),1,[]));
                pijstd(:,:,indPP)   = pij(:,:,indPP) - nanmean(reshape(pij(:,:,indPP),1,[]));
                %{
                if obj.DetectOrDisc
                    pCAFstd(:,:,indPP) = pCAF(:,:,indPP) - nanmean(reshape(pCAF(:,:,indPP),1,[])); % so q are the RT values separating the bins
                    mCAFstd(:,:,indPP) = mCAF(:,:,indPP) - nanmean(reshape(mCAF(:,:,indPP),1,[]));   % and qn are the average number of trials in each bin
                end
               %}
            end
            
            datsum.qstd  = nanstd(qstd,[],3);  % so q are the RT values separating the bins
            datsum.qnstd = nanstd(qnstd,[],3); % and qn are the average number of trials in each bin
            datsum.qmstd = nanstd(qmstd,[],3);
           
            % Average across subjects and put in the structure:
            datsum.pijstd = nanstd(pijstd,[],3);
           
            %{
            if obj.DetectOrDisc
                datsum.pCAFstd = nanstd(pCAFstd,[],3); % so q are the RT values separating the bins
                datsum.mCAFstd = nanstd(mCAFstd,[],3);   % and qn are the average number of trials in each bin
            end
            %}
            
            obj.modelBehaviour.datsum = datsum;
            
            %% Quantile probability plots
            if obj.modelBehaviour.ChiOrG == 1
                %{
                figure; hold on;
                for indCond = 1:size(datsum.q, 1)
                    % plot(datsum.qm(indCond,:),cumsum(datsum.qn(indCond,1:end-1)),'*',...
                    % 'Color',obj.figLayOut.colours(indCond,:),'LineWidth', obj.figLayOut.lineWidth)
                    keyboard
                    legendThis(indCond) = errorbar(datsum.qm(indCond,:), cumsum(datsum.qn(indCond,1:end-1)), datsum.qnstd(indCond,1:end-1), 'CapSize',0);
                    legendThis(indCond).Color = obj.figLayOut.colours(indCond,:);
                    legendThis(indCond).LineWidth = 1;
                    legendThis(indCond).LineStyle = obj.figLayOut.lineType{indCond};
                end
                
                xlim([0 obj.stim.RTdeadLine(end)+0.1])
                set(gca,'XTick',[0:0.2:obj.stim.RTdeadLine(end)])
                set(gca,'XTickLabel',[0:0.2:obj.stim.RTdeadLine(end)]);
                xlabel('Reaction times (sec)');
                
                figInfo = gca;
                line([0 0],  [0 figInfo.YLim(end)], 'Color', 'k', 'LineWidth', 1.5)
                
                ylim([0 figInfo.YLim(end)])
                ylabel('Response count (#)')
                
                xstuff = [fliplr(0:-500:obj.eeg.targetEpoch(1)) 500:500:obj.eeg.targetEpoch(end)]./1000;
                
                xticks(xstuff)
                xticklabels(xstuff);
                xlabel('RT relative to target onset (sec)')
                
                set(gca,'FontSize', obj.figLayOut.letterSize);
                set(gca,'FontName', obj.figLayOut.letterType);
                %}
            elseif  obj.modelBehaviour.ChiOrG == 2
                %% Quantile probability of the proportion plots
                fig1 = figure; hold on;
                
                for indCond = 1:size(datsum.q, 1)
                    CI95 = 1.96.*(datsum.pijstd(indCond,1:numBins)./sqrt(length(obj.ppNames)));
                    h = errorbar(datsum.qm(indCond,1:numBins), cumsum(datsum.pij(indCond,1:numBins)),CI95,...
                        'Color', obj.figLayOut.colours(indCond,:), 'LineStyle', obj.figLayOut.lineType{indCond},...
                        'LineWidth', obj.figLayOut.lineWidth);
                    
                    legendThis(indCond) = h;
                end
                
                limits = max(datsum.qm(:, numBins)) + 0.1;
                xlim([0 limits])
                set(gca,'XTick',[0:0.4:limits limits])
                set(gca,'XTickLabel', [num2str((0:0.4:limits)'); '   '; '   ']);
                xlabel('Reaction times (sec)');
                line([0 0], [0 1.05], 'Color', 'k', 'LineWidth', 1.5)
                ylim([0 1.05])
                
                set(gca,'FontSize', obj.figLayOut.letterSize);
                set(gca,'FontName', obj.figLayOut.letterType);
                
                figMisses = figure; 
                figFA = figure; 
                if ~isempty(obj.stim.FACutOff)
                    figFAcut = figure;
                end
                
                dist = -0.1:0.2/(size(datsum.q, 1)-1):0.1;
                for indCond = 1:size(datsum.q, 1)
                    figure(figMisses);hold on
                    bb = bar(1 + dist(indCond), datsum.pij(indCond, furtherIndex),0.05);
                    bb.FaceColor = obj.figLayOut.colours(indCond,:);
                    CI95 = 1.96.*(datsum.pijstd(indCond,furtherIndex)./sqrt(length(obj.ppNames)));
                    errorbar(1 + dist(indCond),datsum.pij(indCond, furtherIndex), CI95, 'Color', [0 0 0])
                   
                    figure(figFA);hold on
                    bb2 = bar(1 + dist(indCond), datsum.pij(indCond, furtherIndex+1),0.05);
                    bb2.FaceColor = obj.figLayOut.colours(indCond,:);
                    CI95 = 1.96.*(datsum.pijstd(indCond,furtherIndex+1)./sqrt(length(obj.ppNames)));
                    errorbar(1 + dist(indCond), datsum.pij(indCond, furtherIndex+1), CI95, 'Color', [0 0 0])
                    
                    if ~isempty(obj.stim.FACutOff)
                        figure(figFAcut);hold on
                        bb2 = bar(1 + dist(indCond), datsum.pij(indCond, furtherIndex+2),0.05);
                        bb2.FaceColor = obj.figLayOut.colours(indCond,:);
                        CI95 = 1.96.*(datsum.pijstd(indCond,furtherIndex+2)./sqrt(length(obj.ppNames)));
                        errorbar(1 + dist(indCond), datsum.pij(indCond, furtherIndex+2), CI95, 'Color', [0 0 0])
                    end
                end
                
                figure(figMisses);hold on
                ylim([0 1.05])
                xlim([1+(dist(1)-0.2) 1+(dist(end)+0.2)])
                xticks([])
                
                set(gca,'FontSize', obj.figLayOut.letterSize);
                set(gca,'FontName', obj.figLayOut.letterType);
               
                figure(figFA);hold on
                
                ylim([0 1.05])
                
                ylim([0 0.5]);
                yticks(0:0.25:0.5);

                xlim([1+(dist(1)-0.2) 1+(dist(end)+0.2)])
                xticks([])
                
                set(gca,'FontSize', obj.figLayOut.letterSize);
                set(gca,'FontName', obj.figLayOut.letterType);
                if ~isempty(obj.stim.FACutOff)
                    figure(figFAcut);hold on
                    ylim([0 1.05])
                    
                    ylim([0 0.5]);
                    yticks(0:0.25:0.5);
                    ax1 = gca;
                    ax1.YAxis.Visible = 'off'; 
                    xlim([1+(dist(1)-0.2) 1+(dist(end)+0.2)])
                    xticks([])
                    
                    set(gca,'FontSize', obj.figLayOut.letterSize);
                    set(gca,'FontName', obj.figLayOut.letterType);
                end
                
                if obj.DetectOrDisc
                    figCAF = figure; hold on;
                    for indCond = 1:size(datsum.q, 1)
                        CI95 = 1.96.*(datsum.pijstd(indCond,numBins+1:numBins+3)./sqrt(length(obj.ppNames)));
                        h = errorbar(datsum.qm(indCond,numBins+1:numBins+3), (datsum.pij(indCond,numBins+1:numBins+3)),CI95,...
                            'Color', obj.figLayOut.colours(indCond,:), 'LineStyle', obj.figLayOut.lineType{indCond},...
                            'LineWidth', obj.figLayOut.lineWidth);
                        legendThis(indCond) = h;
                    end
                    
                    line([0 0],  [0 1.05], 'Color', 'k', 'LineWidth', 1.5)
                    
                    ylim([0 0.6]);%1.05])
                    yticks([0:0.2:0.6])
                    %ylabel('Proportion')
                    
                    xlim([0 limits])
                    set(gca,'XTick',[0:0.4:limits])
                    set(gca,'XTickLabel', [num2str([0:0.4:limits]')]);
                    xlabel('Reaction times (sec)');
                    
                    set(gca,'FontSize', obj.figLayOut.letterSize);
                    set(gca,'FontName', obj.figLayOut.letterType);
                end
                %}
            end
            
            figure(fig1)
            plotSave(gca, ['QuantileProb' num2str(plotComb) '.png'], figureFolder, obj.figLayOut.saveDim);
            
            legend(legendThis, obj.figLayOut.legends, 'location', 'southeast')
            plotSave(gca, ['QuantileProb' num2str(plotComb) 'Legend.png'], figureFolder, obj.figLayOut.saveDim);

            figure(figMisses)
            plotSave(gca, ['Misses' num2str(plotComb) '.png'], figureFolder, [3.3 3]);
            
            figure(figFA)
            plotSave(gca, ['FalseAlarms' num2str(plotComb) '.png'], figureFolder, [3.3 3]);
            if ~isempty(obj.stim.FACutOff)
                figure(figFAcut)
                plotSave(gca, ['FalseAlarmsCutOff' num2str(plotComb) '.png'], figureFolder, [3.3 3]);
            end 
            
            %{
            if obj.DetectOrDisc
                figure(figCAF)
                plotSave(gca, ['CAF' num2str(plotComb) '.png'], figureFolder, obj.figLayOut.saveDim);
            end
            %}
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% --------------- Supportive functions'    --------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% --------------- optimalization functions    --------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [err, pred] = Chisquared(obj, simdat)
            % Now count the number of simulated trials falling in each of the quantile bins of the real data.
            % 'E' here refers to the "Expected value" from the typical chi-squared formulation.
            % the below exactly mimics how we summarised the Real data, except that
            % here we regard the simulated data as one subject.
            
            datsum     = obj.modelBehaviour.datsum;
            cond       = simdat(:,1:end-2); uniComb = unique(cond, 'rows');
            outcome    = simdat(:,end-1);
            RT         = simdat(:,end);
            
            for indCond = 1:size(uniComb,1)
                currCond    = find(all(cond == uniComb(indCond,:),2));
                currOutcome = outcome(currCond);
                %                 currRT      = RT(currCond);
                %
                %                 currNumTrials = sum(currOutcome ~= 3);
                
                for indBin = 1:size(datsum.q,2)-1
                    E(indCond,indBin) = sum(currOutcome == 1 & RT > datsum.q(indCond,indBin) & RT <= datsum.q(indCond,indBin+1));
                end
                % then add a final bin that counts the false alarms:
                E(indCond,indBin+1) = length(find(cond==indCond & outcome==3));
            end
            
            if any(cond == 3)
                indCond = 3;
                for indCoh = 1:length(cohs)
                    for indBin=1:size(datsum.q,2)-1
                        pred(indCond+indCoh-1,indBin) = length(find(cond==indCond & evidence == cohs(indCoh) & outcome==1 &...
                            RT>datsum.q(indCond+indCoh-1,indBin) & RT<=datsum.q(indCond+indCoh-1,indBin+1)));
                    end
                    
                    % then add a final bin that counts the false alarms:
                    pred(indCond:indCond+1,indBin+1) = length(find(cond==indCond & outcome==3))./2;
                end
            end
            
            pred = pred/obj.modelBehaviour.simulateMoreX/length(obj.ppNames); % now scale down of course, by the factor by which we simulate more blocks than any one subject does
            
            % compute chi-square
            err = 0;
            for indCond = 1:size(pred,1)
                for indBin = 1:size(pred,2)
                    if isnan(pred(indCond,indBin)), continue; end % just because there is one cell of the real and simulated data that is empty - there is no such thing as false alarm number for a "second" coherence level
                    if pred(indCond,indBin) >= 1
                        err = err + (datsum.qn(indCond,indBin)-pred(indCond,indBin))^2/pred(indCond,indBin);
                    else
                        err = err + (datsum.qn(indCond,indBin)-pred(indCond,indBin))^2/1;
                    end
                end
            end
        end
        
        function [err, pred] = Gsquared(obj,simdat)
            % Now count the number of simulated trials falling in each of 
            % the quantile bins of the real data.
            % Express as proportions to compute G^2
            % the below exactly mimics how we summarised the Real data, except that
            % here we regard the simulated data as one subject.
            % check out Linking Theoretical Decision-making Mechanisms in
            % the Simon Task with Electrophysiological Data: A Model-based
            % Neuroscience Study in Humans, by Servant et al.
            
            % TODO change datsum/sim set up and Discrimination
            datsum     = obj.modelBehaviour.datsum;
           
            if obj.DetectOrDisc
                cond       = simdat(:,1:end-3); 
                outcome    = simdat(:,end-2);
                RT         = simdat(:,end-1);
                Response   = simdat(:,end);
                Response(isnan(Response)) = 0;
            else
                cond       = simdat(:,1:end-2);
                outcome    = simdat(:,end-1);
                RT         = simdat(:,end);
            end
            
            uniCond = unique(cond(all(cond ~= 0,2),:), 'rows');
            qps     = obj.modelBehaviour.qps;
            numBins = length(qps) + 1; furtherIndex = numBins+1;

            for indCond = 1:size(uniCond,1)
                currCond    = all(cond == uniCond(indCond,:), 2);
                currOutcome = outcome(currCond);    
                currMisses  = currOutcome == 2;
                currFA      = currOutcome == 3;
                
                currFART    = RT(currCond & outcome == 3);
                currRT      = RT(currCond & outcome == 1);
                
                currNumTrials = sum(currOutcome ~= 3);
                
                if obj.DetectOrDisc  % get only correct responses for the cummeliatve response function.
                    currCorrect = all(cond == uniCond(indCond,:), 2)  & Response;
                    currRT      = RT(currCorrect);
                    
                    currError   = all(cond == uniCond(indCond,:), 2)  & ~Response;
                    currRTError = RT(currError);
                end
                
                    
                for indBin = 1:numBins
                    % "predicted" proportion for this parameter vector - estimated from the simulation
                    pred(indCond,indBin) = min(1, sum(currRT>datsum.q(indCond,indBin) & currRT<=datsum.q(indCond,indBin+1))/currNumTrials);
                end
                                                
                if obj.DetectOrDisc
                    for indQ = 1:3
                        keyboard
                        % number of trials in each of the 5 bins
                        pred(indCond,numBins+indQ)  =  min(1, sum(currRTError>datsum.q(indCond,numBins+1+indQ) & currRTError<=datsum.q(indCond,numBins+1+indQ+1))/currNumTrials);
                    end
                    if numBins ==  furtherIndex
                        furtherIndex = length(pred)+1;
                    end
                end

                pred(indCond,furtherIndex) = min(1, sum(currMisses)/currNumTrials);

                % then add a final bin that counts the false alarms:
                pred(indCond,furtherIndex+1) = min(1, sum(currFA)/(datsum.nShortestITI(indCond)*obj.modelBehaviour.simulateMoreX));

                if ~isempty(obj.stim.FACutOff)
                    currFA = currRT(currOutcome == 3);
                    currFA = currFA > obj.stim.FACutOff;
                    pred(indCond,furtherIndex+2) = min(1, sum(currFA)/currNumTrials);
                end

                % note the min(1, .. limiter makes sure the proportion doesn't go above 1, which causes complex G-squared
                % values. I *think* this also might mean the starting vectors that are really far off in that they produce way
                % too many false alarms, are given up on earlier in the fminsearch procedure, which is a good thing - no point
                % in wasting time on them.
            end
            
            err = [];
            datsum.pij(datsum.pij == 0) = 0.001;
            for indCond = 1:size(pred,1)
                % hits and misses
                for indBin = 1:numBins+1 
                    err = [err datsum.qn(indCond,indBin)*log(datsum.pij(indCond,indBin)/(pred(indCond,indBin)+0.0001))]; % small constant is there to avoid 0 in denominator
                end
                
                if obj.DetectOrDisc
                    keyboard
                    for indBin = 1:size(predCAF,2)
                        err = [err datsum.qCAF(indCond,indBin)*log(datsum.pCAF(indCond,indBin)/(predCAF(indCond,indBin)+0.0001))];
                    end
                end
             
                nShortestITI = datsum.nShortestITI(indCond)/length(obj.ppNames);
                
                indFA = furtherIndex + 1;
               
                err = [err datsum.qn(indCond,indFA)*log(datsum.pij(indCond,indFA)/(pred(indCond, indFA)+0.0001))]; % small constant is there to avoid 0 in denominator
                err = [err (nShortestITI*(1-datsum.pij(indCond,indFA)))*log((1-datsum.pij(indCond,indFA))/((1-pred(indCond, indFA))+0.0001))]; % small constant is there to avoid 0 in denominator
               
                if ~isempty(obj.stim.FACutOff)
                    indFA = furtherIndex + 2;
                    numtarg = datsum.numberOfTrials(indCond);

                    err = [err datsum.qn(indCond,indFA)*log(datsum.pij(indCond,indFA)/(pred(indCond,indFA)+0.0001))]; % small constant is there to avoid 0 in denominator
                    err = [err numtarg*(1-datsum.pij(indCond,indFA))*log((1-datsum.pij(indCond,indFA))/(1-pred(indCond,indFA)+0.0001))]; % small constant is there to avoid 0 in denominator
                end
            end
            
            err = nansum(err)*2;
        end
        
        function [AIC, BIC] = penalize(obj, pm, err)
            f = size(pm,1);
            AIC = err + 2*f;
            BIC = err + f*log(nanmean(obj.numBlocks.*obj.numTrials));
        end
    end
    
    methods(Static)
        function ins = getInstance(transfer)
            persistent instance;
            
            if( ~strcmpi(class(instance), 'dataModelling') )
                instance = dataModelling();
            end
            
            if nargin == 1
                allProp = properties(instance);
                for indProp = 1:length(allProp)
                    if isprop(transfer, allProp{indProp})
                        eval(sprintf('instance.%s = transfer.%s;', allProp{indProp}, allProp{indProp}));
                    end
                end
            end
            
            
            ins = instance;
        end
    end
    
    methods(Access = private)
        function obj = dataModelling()
            % Pre-set all variables.
            obj.inputFolder     = {};
            obj.outputFolder    = {};
            obj.logFolder       = {};
            obj.figFolder       = {};
            
            obj.ppNames         = {};
            obj.conditions      = {};
            obj.condNames       = {};
            obj.numBlocks       = [];
            obj.numTrials       = [];
            obj.order           = {};
            
            obj.stim            = [];
            obj.modelBehaviour  = [];
            obj.options         = [];
            obj.behaviour       = [];
            
            obj.DetectOrDisc   = [];
        end
    end
end
