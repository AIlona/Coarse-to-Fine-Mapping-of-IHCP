% Implementation of method of fundamental solutions 
% and Laplace transform
% Domain is a unit square [0,1]X[0,1]


NS  = 4;        % number of source points
NB  = 4;        % number of BC points in  each direction
K   = 3;        % number of polynomials
delta = 0.1;    % distance between source point and 
                % the boundary of the domain

% time discretization                
Time = 100;      % total time
Ntime = 101;     % number of timesteps
deltaT= Time/(Ntime-1);
tau=linspace(0,Time,Ntime)'; % timesteps

syms x; syms y; syms t;

%%% True solution (function of x,y,t)
T = symfun(t*((x-6)^3 + (y-6)^3)/6, [x y t]);

%%% True RHS (function of x,y,t)
q = symfun((x^3+y^3)/6 -3*(x^2+y^2)+18*(x+y-4)-t*(x+y-12), [x y t]);

%%% Domain %%%
[M,xDom,yDom,xBound,yBound,xSource,ySource]=mesh(NS, NB, delta);

%%% Initial condition %%%
ICValue = T(x,y,0);

coeff =zeros(3*NB+K+M,1);
%Tlapl = zeros(Ntime,1);

Tsol = zeros(Ntime,1);
%Qsol = zeros(Ntime,1);
%Tin  = @(x,y)ICValue;

ksi=0.5;
eta=0.5;
for s=1:Ntime
[psiMONO1,psiMONO2,psiMONO3,psiMONO4,psiMONO5,psiMONO6,phiMONO1,phiMONO2,phiMONO3,phiMONO4,phiMONO5,phiMONO6] = monomials(deltaT);
    %%% Assembling the matrix %%%
    A=LHS(K,NS,NB,M,xDom,yDom,xSource,ySource,xBound,yBound,deltaT);

    %%% Assembling of RHS $$$
    b=RHS(tau(s),K,xDom,yDom,xBound,yBound);

    %%% Coefficients  %%%

    coeff = (A\b);

       
    %%% Solution at current time-step at given point %%%
    Tsol(s,1)=Sol(ksi,eta,deltaT,NS,M,coeff,xSource,ySource, xDom, yDom, psiMONO1, psiMONO2, psiMONO3);
   % Qsol(s,1)=Force(ksi,eta,Tin,coeff,xDom,yDom,NS,M,s,phiMONO1,phiMONO2,phiMONO3);
   % Tin = Tsol(s,1,xSource,ySource, xDom, yDom);  
end
    
    
%TlaplFunc = fit(linspace(1,smax,smax)',Tlapl,'a/(s+b)');
