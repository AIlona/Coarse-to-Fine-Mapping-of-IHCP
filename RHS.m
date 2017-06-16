function b=RHS(s,K,xDom,yDom,xBound,yBound)
% computes the right-hand side for the system using the true solution
% each row is an RHS for given time step

syms x;
b=zeros(length(xDom)+length(xBound)+K,1);
%%% Measurements inside of the domain %%%
%meas=zeros(length(xDom),length(tau));
for i=1:length(xDom)
  %  b(i,1)=(1-exp(-4*s))*(cos(2*(xDom(i))) + cos(2*(yDom(i))));  % measurements
    b(i,1)=s*((xDom(i)-6)^3 + (yDom(i)-6)^3)/6;  % measurements
  %  appr = polyfit(tau,meas(i,:)', 5);
  %  f=func(appr);
  % b(i,1) = laplace(f,x,s);
end

%%% Boundary conditions %%%
%bndry=zeros(length(xBound),length(tau));
for i=1:length(xBound)
    b(i+length(xDom),1)=s*((xBound(i)-6)^3 + (yBound(i)-6)^3)/6;
%     b(i+length(xDom),1)=(1-exp(-4*s))*(cos(2*(xBound(i))) + cos(2*(yBound(i))));
%    appr = polyfit(tau,bndry(i,:)', 5);
%    f=func(appr);
%    b(length(xDom)+i,1) = laplace(f,x,s);
end

%%% Zero part of RHS %%%
for i=1:K
    b(length(xDom)+length(xBound)+i,1)=0;
end

end
