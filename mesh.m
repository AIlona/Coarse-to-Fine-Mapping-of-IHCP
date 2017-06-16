function [M,xDom,yDom,xBound,yBound,xSource,ySource]=mesh(NS, NB, delta)
% returns:
% - M = # of internal points with known measurements
% - xDom, yDom = vectors with coordinates of internal points
% - xBound, yBound = vectors with coordinates of boundary points
% - xSource, ySource = vectors with coordinates of 'source' points



hBound = 1/(NB-1);  % distance between points in each direction
M  =  (NB-2)^2; % number of interpolations points


%%% Domain %%%
[xCoord,yCoord] = meshgrid(linspace(0,1,NB),linspace(0,1,NB));

%%% Points inside %%%
[xDom0, yDom0] = meshgrid(linspace(hBound,1-hBound,NB-2),linspace(hBound,1-hBound,NB-2));
xDom = [xDom0(1,:),xDom0(2,:)]';  % make it to be a vector
yDom = [yDom0(1,:),yDom0(2,:)]';  

%%% Boundary points %%%
xBound =zeros(NB^2-(NB-2)^2,1);
xBound(1:NB,1) = xCoord(1,1);
for k=1:(NB-2)
    xBound(NB+1+(k-1)*(NB-2):NB+1+(k-1)*(NB-2)+NB-3,1)=xCoord(1,2:end-1);
end
xBound(3*NB-3:NB^2-(NB-2)^2,1) =xCoord(1,end);

yBound =zeros(NB^2-(NB-2)^2,1);
yBound(1:NB,1) = yCoord(:,1);
yBound((NB+1):(NB+1+NB-2),1) =yCoord(1,1)*ones();
yBound((NB+1+NB-2):3*NB-3,1) =yCoord(end,1)*ones();
yBound(3*NB-3:NB^2-(NB-2)^2,1) = yCoord(:,1);

%%% Source points %%%
[xS,yS] = meshgrid(linspace(0-delta,1+delta,NS),linspace(0-delta,1+delta,NS));
xSource =zeros(NS^2-(NS-2)^2,1);
xSource(1:NS,1) = xS(1,1);
for k=1:(NS-2)
    xSource(NS+1+(k-1)*(NS-2):NS+1+(k-1)*(NS-2)+NS-3,1)=xS(1,2:end-1);
end
xSource(3*NS-3:NS^2-(NS-2)^2,1) =xS(1,end);

ySource =zeros(NS^2-(NS-2)^2,1);
ySource(1:NS,1) = yS(:,1);
ySource((NS+1):(NS+1+NS-2),1) =yS(1,1)*ones();
ySource((NS+1+NS-2):3*NS-3,1) =yS(end,1)*ones();
ySource(3*NS-3:NS^2-(NS-2)^2,1) = yS(:,1);