function adaptedFAPlot(obj, pred, modelName)

figureFolder = fullfile(obj.figFolder,  'groupAverage/Modelling/', [modelName '/']);
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

datsum = obj.modelBehaviour.datsum;

%% first plot and save FA
figure, hold on
% Quantile probability of the proportion plots

for indCond = 1:3 
    %% Quantile probability of the proportion plots
    if indCond < 3 
        
        CI95 = 1.96.*(datsum.pijstd(indCond,end)./sqrt(length(obj.ppNames)));
        
        bb = bar([indCond+0.2 ], [pred.pij(indCond,end) ], 0.4 ); % datsum.pij(indCond,end)
        set(bb,'FaceColor', obj.figLayOut.colours(indCond,:))
        
        % errorbar(indCond, datsum.pij(indCond,end), CI95,'k','LineWidth',obj.figLayOut.lineWidth)
        
        errorbar(indCond-0.2, datsum.pij(indCond,end),CI95, 'o',...
            'MarkerEdgeColor', [0 0 0],...
            'MarkerFaceColor', obj.figLayOut.colours(indCond,:),...
            'LineStyle', 'none','MarkerSize', 3, 'LineWidth',1,...
            'Color', [0 0 0]);
    else
        CI95 = 1.96.*(nanmean(datsum.pijstd(indCond:end,end))./sqrt(length(obj.ppNames)));
        
        bb = bar([indCond+0.2 ], [nanmean(pred.pij(indCond:end,end)) ] , 0.4);
        set(bb,'FaceColor', obj.figLayOut.colours(indCond,:))
        
        errorbar(indCond-0.2, nanmean(datsum.pij(indCond:end,end),1),CI95, 'o',...
            'MarkerEdgeColor', [0 0 0],...
            'MarkerFaceColor', obj.figLayOut.colours(indCond,:),...
            'LineStyle', 'none','MarkerSize', 3, 'LineWidth',1,...
            'Color', [0 0 0]);
    end
end

figInfoFA = gca;
ylim([0 0.14])
yticks([0 0.1 0.2]);

xticks(1:3);
xlim([0.5 3.5]);

fig2 = gca;
xtickangle( 45 )
xticks(1:indCond)
xticklabels({'Weak', 'Strong', 'Mixed'});

set(gca,'FontSize', obj.figLayOut.letterSize);
set(gca,'FontName', obj.figLayOut.letterType);
                   
% proportional to the whole figure, e.g.
plotSave(gcf, ['QuantileProbFA' modelName '.png'], figureFolder,  [2.8 2.4]);
