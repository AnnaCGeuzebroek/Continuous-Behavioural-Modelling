function plotSave(plotHandle, condition, savePath, dim, plotLegend)

if ~exist('plotLegend', 'var'); plotLegend = 0; end
set(plotHandle, 'HandleVisibility','on');


width  = dim(2); % 10     % Width in inches
height = dim(1);
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'centimeters');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print([savePath condition],'-dtiff','-r300');

end