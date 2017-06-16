function b=RHS_fromData(s,K,M,NB,BndryMeas,intPntsMeas)
%  on each time step computes the right-hand side for the system 
% using the given data with measurements inside of the domian 
% and on its boundary

% 08/15/2015

b=zeros(M+NB+K,1);

%%% Measurements inside of the domain %%%
for i=1:M
     b(i,1)=intPntsMeas(s+3,i);  
end

%%% Boundary conditions %%%
for i=1:NB
    b(i+M,1)=BndryMeas(s+3,i);
end

%%% Zero part of RHS %%%
for i=1:K
    b(M+NB+i,1)=0;
end

end
