function [txt] = myupdatefcn(trash,event)
global channel
pos = get(event,'Position');
dts = get(event.Target,'Tag');
display([dts, ['X: ',num2str(pos(1))], ['Y: ',num2str(pos(2))]]);
txt = dts;
Num = regexp(dts,'\d');
channel = str2num(dts(Num));
%        [' X: ',num2str(pos(1))],[' Y: ',num2str(pos(2))]};
end