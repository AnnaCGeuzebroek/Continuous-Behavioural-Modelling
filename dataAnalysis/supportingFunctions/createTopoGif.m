function createTopoGif(obj, spacing, plotCondition)
% function is created to make a gif of the topoplot and the corresponding
% time traces. We use this to compare the N2 and CPP, in order to get an
% idea of the overlap between the two signals.
keyboard
% get ERP data as the topoplot will be the same for both
trangeBaseline = obj.eeg.epochPlot > obj.eeg.baseline(1) & obj.eeg.epochPlot < obj.eeg.baseline(2);
for applyCSD = [1]
    %% For-loop to extract the topoplots (later obviously also the timeseries
    % for the different conditions)
    % get the while-loop going to get the topoplots average across
    % spacing in steps of 50 sample rates
    window = [-1*spacing/2 spacing/2];
    spacingTarget = window(1):50:1250-window(2);
    spacingResponse = -500:50:200;
    
    for indCPP = 1:length(obj.eeg.cppChan)
        CPPTarget{indCPP}  = nan(1, length(obj.eeg.targetEpoch),   max(obj.numTrials*obj.numBlocks), length(obj.ppNames));
        CPPResponse{indCPP} = nan(1, length(obj.eeg.responseEpoch), max(obj.numTrials*obj.numBlocks), length(obj.ppNames));
    end
    
    N2Target    = nan(1, length(obj.eeg.targetEpoch),   max(obj.numTrials*obj.numBlocks), length(obj.ppNames));
    N2Response  = nan(1, length(obj.eeg.responseEpoch), max(obj.numTrials*obj.numBlocks), length(obj.ppNames));
    topoTarget   = nan(length(spacingTarget), obj.eeg.NumberOfChannels, max(obj.numTrials*obj.numBlocks), length(obj.ppNames));
    topoResponse = nan(length(spacingResponse), obj.eeg.NumberOfChannels, max(obj.numTrials*obj.numBlocks), length(obj.ppNames));
    for indPP = 1:length(obj.ppNames)
        clear ERP
        fprintf(['Now processing participant ' obj.ppNames{indPP} ' to get topoplot\n'])
        currInput  = fullfile(obj.outputFolder,'EEG data', obj.ppNames{indPP});
        
        % load the EEG data.
        if applyCSD
            load(fullfile(currInput,[obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.eeg.timing '.mat']), 'csdERP');
            ERP = csdERP;
        else
            load(fullfile(currInput,[obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.eeg.timing '.mat']), 'ERP');
        end
        
        % Baseline-correct the data for the target and
        % response. Again not for the ERPwhole
        ERP = ERP - repmat(nanmean(ERP(:, trangeBaseline, :),2), 1, size(ERP,2), 1);
        
        % calculated ERP topography
        [target, response] = sortERPs(obj, ERP, indPP, 0, 1);
        
        
        % first the target
        for indTime = 1:length(spacingTarget)
            currRange  = spacingTarget(indTime) + window;
            getTopo = obj.eeg.targetEpoch > currRange(1) & obj.eeg.targetEpoch < currRange(2);
            topoTarget(indTime, :,:,indPP) = squeeze(nanmean(target(:,getTopo,:,1),2));
        end
        
        % first the response
        for indTime = 1:length(spacingResponse)
            currRange  = spacingResponse(indTime) + window;
            getTopo = obj.eeg.responseEpoch > currRange(1) & obj.eeg.responseEpoch < currRange(2);
            topoResponse(indTime, :,:,indPP) = nanmean(response(:,getTopo,:,1),2);
        end
        for indCPP = 1:length(obj.eeg.cppChan)
            CPPTarget{indCPP}(:,:,:,indPP)   = nanmean(target(obj.eeg.cppChan{indCPP},:,:),1);
            CPPResponse{indCPP}(:,:,:,indPP) = nanmean(response(obj.eeg.cppChan{indCPP},:,:),1);
        end
        
        N2Target(:,:,:,indPP)   = nanmean(target(obj.eeg.n2Chan,:,:),1);
        N2Response(:,:,:,indPP) = nanmean(response(obj.eeg.n2Chan,:,:),1);
        
    end
    
    %% topoplot
    if nargin < 4
        channels = [];
    end
    
    currFolder = fullfile(obj.figFolder,  'groupAverage', ['TopoGif_CSD' num2str(applyCSD) '/']);
    if ~exist(currFolder, 'dir'); mkdir(currFolder); end
    
    numCond = 1; plotComb = plotCondition;
    withinVar = {}; withinName = {}; betweenVar = {}; betweenName = {};  withinPar = [];
    
    for indVar = 1:length(plotComb)
        if plotComb(indVar) == 0
            
        elseif length(unique(obj.behaviour{1}.trialVase(:,plotComb(indVar)))) > 1
            withinVar{end+1}  = obj.conditions{plotComb(indVar)};
            withinName{end+1} = strrep(obj.figLayOut.legTitle{plotComb(indVar)},' ','');
            withinPar(indVar) = 1;
            numCond = numCond*length(withinVar{end});
        else
            betweenVar{end+1}  = obj.conditions{plotComb(indVar)};
            betweenName{end+1} = strrep(obj.figLayOut.legTitle{plotComb(indVar)},' ','');
            withinPar(indVar) = 0;
        end
    end
    
    
    if ~isempty(plotCondition)
        for indPP = 1:length(obj.ppNames)
            
            currCond = unique(obj.behaviour{indPP}.trialVase(:,plotComb), 'rows'); %(withinPar)
            for indCond = 1:length(currCond)
                theseCond = all(obj.behaviour{indPP}.trialVase(:,plotComb) == currCond(indCond,:),2);
                newTarget(indCond,:,:,indPP)   = nanmean(topoTarget(:,:,theseCond,indPP),3);
                newResponse(indCond,:,:,indPP) = nanmean(topoResponse(:,:,theseCond,indPP),3);
                
                for indCPP = 1:length(obj.eeg.cppChan)
                plotCPPTarget{indCPP}(indCond,:,indPP)    = squeeze(nanmean(CPPTarget{indCPP}(:,:,theseCond,indPP),3));
                plotCPPResponse{indCPP}(indCond,:,indPP)  = squeeze(nanmean(CPPResponse{indCPP}(:,:,theseCond, indPP),3));
                end
                plotN2Target(indCond,:,indPP)     = squeeze(nanmean(N2Target(:,:,theseCond, indPP),3));
                plotN2Response(indCond,:,indPP)   = squeeze(nanmean(N2Response(:,:,theseCond,indPP),3));
                
                if sum(~withinPar)
                    betweenTable(indPP) = unique(obj.behaviour{indPP}.trialVase(:,plotComb(~withinPar)));
                else
                    betweenTable(indPP)  = 1;
                end
            end
        end
    else
        keyboard
        newTarget   = squeeze(nanmean(topoTarget,3));
        newResponse = squeeze(nanmean(topoTarget,3));
        
        betweenTable = ones(1, length(obj.ppNames));
    end
    
    limits = [min(reshape(nanmean(newTarget(:,:,:,:),4),[],1))...
        max(reshape(nanmean(newTarget(:,:,:,:),4),[],1))];
    keyboard
    limitsDiff = [min(reshape((nanmean(nanmean(newTarget([1 3],:,:,:),4),1) - ...
        nanmean(nanmean(newTarget([2 4],:,:,:),4),1)),[],1))...
        max(reshape((nanmean(nanmean(newTarget([1 3],:,:,:),4),1) - ...
        nanmean(nanmean(newTarget([2 4],:,:,:),4),1)),[],1))];
    for indGroup = unique(betweenTable)
        %for indCond = 1:size(newResponse,1)
        h1 = figure;
        fileName1 = fullfile(currFolder, ['topoResponse_Group' num2str(indGroup) '.gif']);
        figName1  = fullfile(currFolder, ['topoResponse_Group' num2str(indGroup)]);
        
        axis tight manual % this ensures that getframe() returns a consistent size
        h2 = figure;
        figName2  = fullfile(currFolder, ['topoResponseDiff1_Group' num2str(indGroup)]);
        axis tight manual % this ensures that getframe() returns a consistent size
         h3 = figure;
        figName3  = fullfile(currFolder, ['topoResponseDiff2_Group' num2str(indGroup)]);
        axis tight manual % this ensures that getframe() returns a consistent size
        
        for indTime = 1:size(newResponse,2)
            figure(h1)
            subplot(2,2,[1 3])
            topoplot(nanmean(nanmean(newResponse(:,indTime,:,betweenTable == indGroup),4),1),...
                obj.eeg.chanlocs, 'style', 'map', 'electrodes', 'off',...
                'maplimits', limits);
            
            subplot(2,2, 2);
            
            lineshere = plot(obj.eeg.responseEpoch, nanmean(plotCPPResponse{1}(:, :, betweenTable == indGroup),3)); hold on
            for indLine = 1:length(lineshere)
                lineshere(indLine).Color = obj.figLayOut.colours(indLine,:);
            end
            currAx = gca;
            area = fill([spacingResponse(indTime) + window fliplr(spacingResponse(indTime) + window)],...
                [currAx.YLim(1) currAx.YLim(1) currAx.YLim(2) currAx.YLim(2)], 'k');
            area.LineStyle = 'none';
            area.FaceColor = [0.7 0.7 0.7];
            area.FaceAlpha = 0.3;
            xlim([spacingResponse(1) spacingResponse(end)])
            
            subplot(2,2, 4);
            lineshere = plot(obj.eeg.responseEpoch, nanmean(plotN2Response(:, :, betweenTable == indGroup),3)); hold on
            for indLine = 1:length(lineshere)
                lineshere(indLine).Color = obj.figLayOut.colours(indLine,:);
            end
            currAx = gca;
            area = fill( [spacingResponse(indTime) + window fliplr(spacingResponse(indTime) + window)],...
                [currAx.YLim(1) currAx.YLim(1) currAx.YLim(2) currAx.YLim(2)], 'k');
            area.LineStyle = 'none';
            area.FaceColor = [0.7 0.7 0.7];
            area.FaceAlpha = 0.3;
            xlim([spacingResponse(1) spacingResponse(end)])
            
            % Capture the plot as an image
            frame = getframe(h1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if indTime == 1
                imwrite(imind,cm,fileName1,'gif', 'Loopcount',inf, 'DelayTime', 0.5);
            else
                imwrite(imind,cm,fileName1,'gif','WriteMode','append', 'DelayTime', 0.5);
            end
            
           if spacingResponse(indTime) == 0
                imwrite(imind,cm,[figName1 '_' num2str(spacingResponse(indTime)) '.jpg'],'jpg');
                
                figure(h2); clf;
                topoplot(nanmean(nanmean(newResponse(1,indTime,:,betweenTable == indGroup),4),1) - ...
                    nanmean(nanmean(newResponse(2,indTime,:,betweenTable == indGroup),4),1),...
                    obj.eeg.chanlocs, 'style', 'map', 'electrodes', 'off',...
                    'maplimits', limitsDiff);
              
                % Capture the plot as an image
                frame = getframe(h2);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                imwrite(imind,cm,[figName2 '_2' num2str(spacingResponse(indTime)) '.jpg'],'jpg')  
                
                figure(h3); clf;
                topoplot(nanmean(nanmean(newResponse(3,indTime,:,betweenTable == indGroup),4),1) - ...
                    nanmean(nanmean(newResponse(4,indTime,:,betweenTable == indGroup),4),1),...
                    obj.eeg.chanlocs, 'style', 'map', 'electrodes', 'off',...
                    'maplimits', limitsDiff);
                   % Capture the plot as an image
                frame = getframe(h3);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                imwrite(imind,cm,[figName3 '_2' num2str(spacingResponse(indTime)) '.jpg'],'jpg') 
            end
            figure(h1);
            clf
        end
        
        %end
    end
    
    %
    for indGroup = unique(betweenTable)
        %for indCond = 1:size(newResponse,1)
        h1 = figure;
        fileName1 = fullfile(currFolder, ['topoTarget_Group' num2str(indGroup) '.gif']);
        figName1  = fullfile(currFolder, ['topoTarget_Group' num2str(indGroup)]);
        
        axis tight manual % this ensures that getframe() returns a consistent size
        h2 = figure;
        figName2  = fullfile(currFolder, ['topoTargetDiff_Group' num2str(indGroup)]);
        axis tight manual % this ensures that getframe() returns a consistent size
        
        for indTime = 1:size(newTarget,2)
            figure(h1)
            subplot(2,2,[1 3])
            topoplot(nanmean(nanmean(newTarget(:,indTime,:,betweenTable == indGroup),4),1),...
                obj.eeg.chanlocs, 'style', 'map', 'electrodes', 'off',...
                'maplimits', limits);
            
            subplot(2,2, 2);
            lineshere = plot(obj.eeg.targetEpoch, nanmean(plotCPPTarget{1}(:, :, betweenTable == indGroup),3)); hold on
            for indLine = 1:length(lineshere)
                lineshere(indLine).Color = obj.figLayOut.colours(indLine,:);
            end
            currAx = gca;
            area = fill([spacingTarget(indTime) + window fliplr(spacingTarget(indTime) + window)],...
                [currAx.YLim(1) currAx.YLim(1) currAx.YLim(2) currAx.YLim(2)], 'k');
            area.LineStyle = 'none';
            area.FaceColor = [0.7 0.7 0.7];
            area.FaceAlpha = 0.3;
            xlim([spacingTarget(1) spacingTarget(end)])
            
            subplot(2,2, 4);
            lineshere = plot(obj.eeg.targetEpoch, nanmean(plotN2Target(:, :, betweenTable == indGroup),3)); hold on
            for indLine = 1:length(lineshere)
                lineshere(indLine).Color = obj.figLayOut.colours(indLine,:);
            end
            currAx = gca;
            area = fill( [spacingTarget(indTime) + window fliplr(spacingTarget(indTime) + window)],...
                [currAx.YLim(1) currAx.YLim(1) currAx.YLim(2) currAx.YLim(2)], 'k');
            area.LineStyle = 'none';
            area.FaceColor = [0.7 0.7 0.7];
            area.FaceAlpha = 0.3;
            xlim([spacingTarget(1) spacingTarget(end)])
            
            % Capture the plot as an image
            frame = getframe(h1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if indTime == 1
                imwrite(imind,cm,fileName1,'gif', 'Loopcount',inf, 'DelayTime', 0.5);
            else
                imwrite(imind,cm,fileName1,'gif','WriteMode','append', 'DelayTime', 0.5);
            end
            
            if spacingTarget(indTime) >= 550 && spacingTarget(indTime) < 900
                imwrite(imind,cm,[figName1 '_' num2str(spacingTarget(indTime)) '.jpg'],'jpg');
                
                figure(h2); clf;
                topoplot(nanmean(nanmean(newTarget(1,indTime,:,betweenTable == indGroup),4),1) - ...
                    nanmean(nanmean(newTarget(2,indTime,:,betweenTable == indGroup),4),1),...
                    obj.eeg.chanlocs, 'style', 'map', 'electrodes', 'off',...
                    'maplimits', limitsDiff);
                % Capture the plot as an image
                frame = getframe(h2);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                imwrite(imind,cm,[figName2 '_' num2str(spacingTarget(indTime)) '.jpg'],'jpg')
            end
            figure(h1);
            clf
        end
        
        %end
    end
    %}
end