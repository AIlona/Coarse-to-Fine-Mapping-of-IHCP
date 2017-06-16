% Main file, works in both cases: analytic function / measurements 
% 08/15/2015

 K    = 3;        % number of polynomials 
 delta = 0.1;     % distance between source point and 
                  % the boundary of the domain

                  
%%% Configuration of the domain. Use if points are not given %%%
%  ax = 0;         % LHS point in x-direction 
%  bx = 1;         % RHS point in x-direction
%  ay = 0;         % LHS point in y-direction 
%  by = 1;         % RHS point in y-direction
%  NSx = 4;        % number of source points in x-direction
%  NSy = 4;        % number of source points in y-direction
%  NBx = 4;        % number of BC points in x-direction
%  NBy = 4;        % number of BC points in y-direction
% [M,NB,NS,xDom,yDom,xBound,yBound,xSource,ySource]=mesh(ax,bx,ay,by,NSx,NSy,NBx,NBy,delta);


%%% Configuration of the domain and measurements %%%
%%% Use in case of given measurements %%%
%%% Load files if use given data %%%
load 'Moving2D_bottom_boundary.txt'
load 'Moving2D_interior_points.txt'
load 'Moving2D_right_boundary.txt'
load 'Moving2D_top_boundary.txt'
load 'Moving2D_left_boundary.txt'

% Particular case 1: 6X8 grid
% BndryMeas1 = [Moving2D_bottom_boundary(:,1),Moving2D_left_boundary,Moving2D_top_boundary(:,1),Moving2D_bottom_boundary(:,2),Moving2D_top_boundary(:,2),Moving2D_bottom_boundary(:,3),Moving2D_top_boundary(:,3),Moving2D_bottom_boundary(:,4),Moving2D_top_boundary(:,4),Moving2D_bottom_boundary(:,5),Moving2D_top_boundary(:,5),Moving2D_bottom_boundary(:,6),Moving2D_top_boundary(:,6),Moving2D_bottom_boundary(:,7),Moving2D_top_boundary(:,7),Moving2D_bottom_boundary(:,8),Moving2D_top_boundary(:,8),Moving2D_bottom_boundary(:,9),Moving2D_right_boundary,Moving2D_top_boundary(:,9)];
% BndryMeas =[BndryMeas1(:,1:2),BndryMeas1(:,4),BndryMeas1(:,6),BndryMeas1(:,8:15),BndryMeas1(:,18:25),BndryMeas1(:,27),BndryMeas1(:,29),BndryMeas1(:,31:end)];
% intPntsMeas1 = Moving2D_interior_points;
% intPntsMeas = [intPntsMeas1(:,1:3),intPntsMeas1(:,5:7),intPntsMeas1(:,15:17),intPntsMeas1(:,19:21),intPntsMeas1(:,29:31),intPntsMeas1(:,33:35),intPntsMeas1(:,43:45),intPntsMeas1(:,47:49)];
% M  = size(intPntsMeas,2);
% NB = size(BndryMeas,2);
% NS = NB;
% xDom = intPntsMeas(1,:)';
% yDom = intPntsMeas(2,:)';
% xBound = BndryMeas(1,:)';
% yBound = BndryMeas(2,:)';
% [xSource,ySource]=sourcePnts(0,1,0,1,8,6, delta);

% Particular case 2: 8X9 grid
BndryMeas1 = [Moving2D_bottom_boundary(:,1),Moving2D_left_boundary,Moving2D_top_boundary(:,1),Moving2D_bottom_boundary(:,2),Moving2D_top_boundary(:,2),Moving2D_bottom_boundary(:,3),Moving2D_top_boundary(:,3),Moving2D_bottom_boundary(:,4),Moving2D_top_boundary(:,4),Moving2D_bottom_boundary(:,5),Moving2D_top_boundary(:,5),Moving2D_bottom_boundary(:,6),Moving2D_top_boundary(:,6),Moving2D_bottom_boundary(:,7),Moving2D_top_boundary(:,7),Moving2D_bottom_boundary(:,8),Moving2D_top_boundary(:,8),Moving2D_bottom_boundary(:,9),Moving2D_right_boundary,Moving2D_top_boundary(:,9)];
BndryMeas =[BndryMeas1(:,1:4),BndryMeas1(:,6:27),BndryMeas1(:,29:end)];
intPntsMeas1 = Moving2D_interior_points;
intPntsMeas = [intPntsMeas1(:,1:24),intPntsMeas1(:,29:49)];
M  = size(intPntsMeas,2);
NB = size(BndryMeas,2);
NS = NB;
xDom = intPntsMeas(1,:)';
yDom = intPntsMeas(2,:)';
xBound = BndryMeas(1,:)';
yBound = BndryMeas(2,:)';
[xSource,ySource]=sourcePnts(0,1,0,1,9,8, delta);


%%% time discretization %%%   
theta  = 1;                  % parameter of discretization
Time   = 10;                 % total time
Ntime  = 1001;               % number of time steps
deltaT = Time/(Ntime-1);     % size of time step
tau=linspace(0,Time,Ntime)'; % all timesteps
lmbda = 1/sqrt(theta*deltaT);

%%% Initial condition %%%
ICValue = 0.0;

coeff =zeros(NB+K+M,1);
Tsol = zeros(Ntime-1,1);
Qsol = zeros(Ntime-1,1);

%%% For iterations %%%
InT  = ICValue;
InLapl = 0;

%%% Choice of testing point %%%
ksi=intPntsMeas1(1,25);
%ksi=0.2;
%eta=0.3;
eta=intPntsMeas1(2,25);

%%% Solution %%%
for s=1:Ntime-1
    %%% Assembling the matrix %%%
    A=LHS(K,NS,NB,M,xDom,yDom,xSource,ySource,xBound,yBound,lmbda);

    %%% Assembling of RHS %%%
    %b=RHS(tau(s),K,M,NB,xDom,yDom,xBound,yBound);   % true solution
    b=RHS_fromData(s,K,M,NB,BndryMeas,intPntsMeas);  % given data

    %%% Coefficients (solving the linear system) %%%
    coeff = A\b;

       
    %%% Solution at current time-step at given point %%%
    Tsol(s,1)=Sol(ksi,eta,lmbda,NS,M,coeff,xSource,ySource, xDom, yDom, K);
    [Qsol(s,1),InLapl]=Force(ksi,eta,InT,coeff,xDom,yDom,NS,M,K, deltaT, lmbda,theta,InLapl,Tsol(s,1));
    InT=Tsol(s,1);
end
  Tsol = [ICValue;Tsol];
 
    

