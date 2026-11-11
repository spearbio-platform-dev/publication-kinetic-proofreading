function fun_PlotWithSimpleCI(hAxis,xVals,yVals,vfTop,vfBot,viColor)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
hold(hAxis,'on');
vfInBetweenX = [xVals,fliplr(xVals)];
vfInBetweenY = [vfTop,fliplr(vfBot)];

plot(hAxis,xVals,yVals,'Color',viColor,'LineWidth',2);
fill(hAxis,vfInBetweenX,vfInBetweenY,viColor,'EdgeColor','none','FaceAlpha',0.15);

end

