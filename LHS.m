function A=LHS(K,NS,NB,M,xDom,yDom,xSource,ySource,xBound,yBound,lmbda)
% assembles the left-hand side matrix A for the problem
% only Dirichlet BC are implemented 

% Inputs:
% K  = # of monomials
% NS = total # of source points
% NB = total # of points on the boundary
% M  = # of points inside of the domain
% xDom, yDom - vectors with x and y coordinates of points inside of the
% domain resectiely
% xSource, ySource - vectors with x and y coordinates of points inside of the
% domain respectively
% xBound, yBound - vectors with x and y coordinates of points on the
% boundary respectively
% lmbda = 1/sqrt(theta*deltaT), where
% deltaT is the size of the time step and
% theta is a parameter, corresponding to discretization in time (0=<theta=<1)

% 08/15/2015

A=zeros(K+NB+M,K+NS+M);

%%% W1&W2 blocks %%%
for j=1:NS
    for i=1:M
        A(i,j) = besselk(0,lmbda*sqrt(((xDom(i,1)-xSource(j,1))^2+(yDom(i,1)-ySource(j,1))^2)));
    end
    for i=(M+1):(NB+M)
        A(i,j) = besselk(0,lmbda*sqrt(((xBound(i-M,1)-xSource(j,1))^2+(yBound(i-M,1)-ySource(j,1))^2)));
    end
end

%%% Beta1&Beta2 blocks %%%
for i=1:M
    [psiMONO,phiMONO] = monomials(lmbda,xDom(i,1),yDom(i,1),K);
    for j=1:K
        A(i,NS+M+j) = psiMONO(j);
    end
end
for i=(M+1):(NB+M)
    [psiMONO,phiMONO] = monomials(lmbda,xBound(i-M,1),yBound(i-M,1),K);
    for j=1:K
        A(i,NS+M+j) = psiMONO(j);
    end
end
    
 %%% Alpha1&Alpha2&Alpha3 blocks
 for j=1:M
     for i=1:M
         rm = sqrt((xDom(i,1)-xDom(j,1))^2+(yDom(i,1)-yDom(j,1))^2);
         phiHat = 0.0;
         if rm>0
             phiHat = -(4/lmbda^4)*(besselk(0,rm*lmbda)+log(rm)+1)-(rm^2*log(rm))/lmbda^2;
         else
             phiHat = (4/lmbda^4)*(0.57721+log(lmbda/2)-1);
         end
         A(i,NS+j)=phiHat;
     end
     for i=(M+1):(M+NB)
         rm = sqrt((xBound(i-M,1)-xDom(j,1))^2+(yBound(i-M,1)-yDom(j,1))^2);
         phiHat = 0.0;
         if rm>0
             phiHat = -(4/lmbda^4)*(besselk(0,rm*lmbda)+log(rm)+1)-(rm^2*log(rm))/lmbda^2;
         else
             phiHat = (4/lmbda^4)*(0.57721+log(lmbda/2)-1);
         end
         A(i,NS+j)=phiHat;
     end
 end
 for i=1:M
     [psiMONO,phiMONO] = monomials(lmbda,xDom(i,1),yDom(i,1),K);
     for j=1:K
        A(M+NB+j,i+NS) = phiMONO(j);
     end
 end
