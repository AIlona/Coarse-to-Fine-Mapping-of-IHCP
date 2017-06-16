function [Qnm1,InLapl]=Force(x,y,InT,coeff,xDom,yDom,NS,M,K,deltaT,lmbda, theta, InLapl,Tsol)
% returns:
% Qnm1 = value of source function at previous time step at given point
% InLapl = value of laplacian of T at current time step at given point
% inputs:
% x,y = coordinates of given point
% InT = value of temperature at given point at previous time step
% coeff = vector with coefficiens W's, alpha's, beta's computed on 
% current time step for given point
% xDom, yDom - vectors with x and y coordinates of points inside of the
% domain resectiely
% K  = # of monomials
% NS = total # of source points
% M  = # of points inside of the domain
% deltaT is the size of the time step
% theta is a parameter, corresponding to discretization in time (0=<theta=<1)
% lmbda = 1/sqrt(theta*deltaT)
% InLapl = laplacian of temperature at given point on previous time step
% Tsol = temperature vlue at given point on currecnt time step

% 08/15/2015

fnm1=0;
for i=1:M
    rm = sqrt((x-xDom(i,1))^2 + (y-yDom(i,1))^2);
    phihat = log(rm)*rm^2;
    fnm1=fnm1+coeff(i+NS,1)*phihat;
end
[psiMONO,phiMONO] = monomials(lmbda,x,y,K);
for j=1:K
    fnm1 = fnm1 + coeff(NS+M+j,1)*phiMONO(j);
end

Qnm1 = -theta*fnm1 - InT/deltaT - (1-theta)*InLapl;
InLapl = fnm1+(1/(theta*deltaT))*Tsol;