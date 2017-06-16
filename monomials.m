function [psiMONO,phiMONO] = monomials(lmbda,x,y,K)
% returns values of the first few monomials and their particular
% solutions at given point

%Inputs:
% lmbda = 1/sqrt(theta*deltaT), where
% deltaT is the size of the time step and
% theta is a parameter, corresponding to discretization in time (0=<theta=<1)
% x,y - coordinates of given point
% K - desired # of monomials (K=<10)

% 08/15/2015
 
%%% psi functions %%%
psiMONO1=-1/lmbda^2;
psiMONO2=-x/lmbda^2;
psiMONO3=-y/lmbda^2;
psiMONO4=-x*y/lmbda^2;
psiMONO5=-x^2/lmbda^2-2/lmbda^4;
psiMONO6=-y^2/lmbda^2-2/lmbda^4;
psiMONO7=-x^3/lmbda^2-6*x/lmbda^4;
psiMONO8=-x^2*y/lmbda^2-2*y/lmbda^4;
psiMONO9=-x*y^2/lmbda^2-2*x/lmbda^4;
psiMONO10=-y^3/lmbda^2-6*y/lmbda^4;
psiMONOall = [psiMONO1;psiMONO2;psiMONO3;psiMONO4;psiMONO5;psiMONO6;psiMONO7;psiMONO8;psiMONO9;psiMONO10];
psiMONO = psiMONOall(1:K);

%%% phi functions %%%
phiMONO1= 1;
phiMONO2= x;
phiMONO3= y;
phiMONO4=x*y;
phiMONO5=x^2;
phiMONO6=y^2;
phiMONO7=x^3;
phiMONO8=x^2*y;
phiMONO9=x*y^2;
phiMONO10=y^3;
phiMONOall = [phiMONO1;phiMONO2;phiMONO3;phiMONO4;phiMONO5;phiMONO6;phiMONO7;phiMONO8;phiMONO9;phiMONO10];
phiMONO = phiMONOall(1:K);

end