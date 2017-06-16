function A=LHS(K,NS,NB,M,xDom,yDom,xSource,ySource,xBound,yBound,deltaT)
% assembles the left-hand side matrix for a given value of parameter s


[psiMONO1,psiMONO2,psiMONO3,phiMONO1,phiMONO2,phiMONO3] = monomials(deltaT);

% 4NB-4 for boundary points, same for source!
A=zeros(K+3*NB+M,K+3*NB+M);
%%% W1&W2 blocks
for j=1:3*NS
    for i=1:M
        A(i,j) = besselk(0,sqrt((1/deltaT)*((xDom(i,1)-xSource(j,1))^2+(yDom(i,1)-ySource(j,1))^2)));
    end
    for i=(M+1):(3*NB+M)
        A(i,j) = besselk(0,sqrt((1/deltaT)*((xBound(i-M,1)-xSource(j,1))^2+(yBound(i-M,1)-ySource(j,1))^2)));
    end
end

%%% Beta1&Beta2 blocks
for i=1:M
        A(i,3*NS+M+1) = psiMONO1(xDom(i,1),yDom(i,1));
        A(i,3*NS+M+2) = psiMONO2(xDom(i,1),yDom(i,1));
        A(i,3*NS+M+3) = psiMONO3(xDom(i,1),yDom(i,1));
end
for i=(M+1):(3*NB+M)
        A(i,3*NS+M+1) = psiMONO1(xBound(i-M,1),yBound(i-M,1));
        A(i,3*NS+M+2) = psiMONO2(xBound(i-M,1),yBound(i-M,1));
        A(i,3*NS+M+3) = psiMONO3(xBound(i-M,1),yBound(i-M,1));
end
    
 %%% Alpha1&Alpha2&Alpha3 blocks
 for j=1:M
     for i=1:M
         rm = sqrt((xDom(i,1)-xDom(j,1))^2+(yDom(i,1)-yDom(j,1))^2);
         phiHat = 0.0;
         if rm>0
             phiHat = -(4*deltaT^2)*(besselk(0,rm*sqrt(1/deltaT))+log(rm)+1)-(rm^2*log(rm))*deltaT;
         else
             phiHat = (4*deltaT^2)*(0.57721+log(sqrt(1/deltaT)/2)-1);
         end
         A(i,3*NS+j)=phiHat;
     end
     for i=(M+1):(M+3*NB)
         rm = sqrt((xBound(i-M,1)-xDom(j,1))^2+(yBound(i-M,1)-yDom(j,1))^2);
         phiHat = 0.0;
         if rm>0
             phiHat = -(4*deltaT^2)*(besselk(0,rm*sqrt(1/deltaT))+log(rm)+1)-(rm^2*log(rm))*deltaT;
         else
             phiHat = (4*deltaT^2)*(0.57721+log(sqrt(1/deltaT)/2)-1);
         end
         A(i,3*NS+j)=phiHat;
     end
 end
 for i=1:M
        A(M+3*NB+1,i+3*NS) = phiMONO1(xDom(i,1),yDom(i,1));
        A(M+3*NB+2,i+3*NS) = phiMONO2(xDom(i,1),yDom(i,1));
        A(M+3*NB+3,i+3*NS) = phiMONO3(xDom(i,1),yDom(i,1));
 end
