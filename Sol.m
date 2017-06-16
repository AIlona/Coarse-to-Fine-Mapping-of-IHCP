function Tsol=Sol(x,y,lmbda,NS,M,coeff,xSource,ySource, xDom, yDom,K)
% computes the temperature on given time step using coefficients W,
% alpha, beta at given point (x,y)
% Inputs:

% x,y = x and y coordinates of given point
% coeff = vector of coefficients, coming from solutionof linear system
% M = total # of points inside of the domain
% xDom, yDom - vectors with x and y coordinates of points inside of the
% domain respectively
% NS = total # of source points
% xSource, ySource - vectors with x and y coordinates of source points
% K = # of monomials
% lmbda = 1/sqrt(theta*deltaT), where
% deltaT is the size of the time step and
% theta is a parameter, corresponding to discretization in time (0=<theta=<1)

% 08/15/2015

[psiMONO,phiMONO] = monomials(lmbda,x,y,K);
Tsol=0;
    for i=1:NS
        Tsol= Tsol + coeff(i,1)*besselk(0,lmbda*sqrt(((x-xSource(i,1))^2+(y-ySource(i,1))^2)));
    end
    for i=NS+1:NS+M
        rm = sqrt((x-xDom(i-NS,1))^2+(y-yDom(i-NS,1))^2);
        psiHat = 0.0;
             if rm>0
                 psiHat = -4/lmbda^4*(besselk(0,rm*lmbda)+log(rm)+1)-(rm^2*log(rm))/lmbda^2;
             else
                 psiHat = 4/lmbda^4*(0.57721+log(lmbda/2)-1);
             end
        Tsol=Tsol+coeff(i,1)*psiHat;
    end
    for j=1:K    
         Tsol = Tsol + coeff(NS+M+j,1)*psiMONO(j);
    end

