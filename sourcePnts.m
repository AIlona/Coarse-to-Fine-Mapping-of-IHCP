function [xSource,ySource]=sourcePnts(ax,bx,ay,by,NSx,NSy, delta)
% returns:
% xSource, ySource = vectors with x and y coordinates of 'source' points
% Inputs:
% [ax,bx] = length of the domain in x-direction
% [ay,by] = length of the domain in y-direction
% NSx, NSy = # of source points in x and y directions
% delta = distance between the physical domain and 'source contour'

% 08/16/2015

%%% Source points %%%
[xS,yS] = meshgrid(linspace(ax-delta,bx+delta,NSx),linspace(ay-delta,by+delta,NSy));
xSource = zeros(2*NSy+(NSx-2)*2,1);
xSource(1:NSy,1) = xS(1,1);
for k=1:NSx-2
    xSource(NSy+2*k-1,1)=xS(1,k+1);
    xSource(NSy+2*k,1)=xS(1,k+1);
end
xSource(NSy+2*(NSx-2)+1:end,1) =xS(1,end);

ySource =zeros(2*NSy+(NSx-2)*2,1);
ySource(1:NSy,1) = yS(:,1);
for k=1:NSx-2
    ySource(NSy+2*k-1,1) =yS(1,1);
    ySource(NSy+2*k,1) =yS(end,1);
end
ySource(NSy+2*(NSx-2)+1:end,1) = yS(:,1);
end