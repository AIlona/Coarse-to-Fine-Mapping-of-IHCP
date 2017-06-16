function b=RHS(s,K,M,NB,xDom,yDom,xBound,yBound)
% computes the right-hand side for the system on each time step
% using the true solution

% Inputs:
% s is the time index
% K  = # of monomials
% M  = # of points inside of the domain
% NB = total # of points on the boundary
% xDom, yDom - vectors with x and y coordinates of points inside of the
% domain resectiely
% xBound, yBound - vectors with x and y coordinates of points on the
% boundary respectively

% 08/15/2015


b=zeros(M+NB+K,1);

%%% Measurements inside of the domain %%%
% depends on given true RHS function, see  examples below
for i=1:M
     b(i,1) = sin(xDom(i))*sin(yDom(i))*sin(s);  
  %  b(i,1) = s*((xDom(i)-6)^3 + (yDom(i)-6)^3)/6;
  %  b(i,1)=(1-exp(-4*s))*(cos(2*(xDom(i))) + cos(2*(yDom(i))));  
end

%%% Boundary conditions %%%
% depends on given true RHS function, see  examples below
for i=1:NB
    b(i+M,1)=sin(xBound(i))*sin(yBound(i))*sin(s);
 %  b(i+M,1)=s*((xBound(i)-6)^3 + (yBound(i)-6)^3)/6;   
 %  b(i+length(xDom),1)=(1-exp(-4*s))*(cos(2*(xBound(i))) + cos(2*(yBound(i))));

end

%%% Zero part of RHS %%%
for i=1:K
    b(M+NB+i,1)=0;
end

end
