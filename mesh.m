function [M,NB,NS,xDom,yDom,xBound,yBound,xSource,ySource]=mesh(ax,bx,ay,by,NSx,NSy, NBx ,NBy , delta)
% returns:
% - M = # of internal points with known measurements
% - xDom, yDom = vectors with coordinates of internal points
% - xBound, yBound = vectors with coordinates of boundary points
% - xSource, ySource = vectors with coordinates of 'source' points
% Inputs:
% [ax,bx] = length in x-direction
% [ay,by] = length in y-direction
% NBx, NBy = # of points on each horizontal/vertical boundary, respectively
% NSx, NSy = # of source points in x and y directions
% delta = distance between the physical domain and 'source contour'

% 08/15/2015

hx = (bx-ax)/(NBx-1);        % distance between points in x-direction
hy = (by-ay)/(NBy-1);        % distance between points in y-direction
M  =  (NBx-2)*(NBy-2);       % number of interpolations points inside


%%% Domain (bndry + inside) %%%
[xCoord,yCoord] = meshgrid(linspace(ax,bx,NBx),linspace(ay,by,NBy));

%%% Points inside of the domain %%%
[xDom0, yDom0] = meshgrid(linspace(ax+hx,bx-hx,NBx-2),linspace(ay+hy,by-hy,NBy-2));
xDom=reshape(xDom0.',1,[])';
yDom=reshape(yDom0.',1,[])';


%%% Boundary points %%%
xBound =zeros(2*NBy+(NBx-2)*2,1);
xBound(1:NBy,1) = xCoord(1,1);
for k=1:NBx-2
    xBound(NBy+2*k-1,1)= xCoord(1,k+1);
    xBound(NBy+2*k,1)  = xCoord(1,k+1);
end
xBound(NBy+2*(NBx-2)+1:end,1) =xCoord(1,end);

yBound =zeros(2*NBy+(NBx-2)*2,1);
yBound(1:NBy,1) = yCoord(:,1);
for k=1:NBx-2
    yBound(NBy+2*k-1,1) = yCoord(1,1);
    yBound(NBy+2*k,1)   = yCoord(end,1);
end
yBound(NBy+2*(NBx-2)+1:end,1) = yCoord(:,1);
NB = length(xBound);



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
NS = length(xSource);
end