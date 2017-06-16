function Tsol=Sol(x,y,deltaT,NS,M,coeff,xSource,ySource, xDom, yDom, psiMONO1, psiMONO2, psiMONO3)

Tsol=0;
    for i=1:3*NS
        Tsol= Tsol + coeff(i,1)*besselk(0,sqrt((1/deltaT)*((x-xSource(i,1))^2+(y-ySource(i,1))^2)));
    end
    for i=3*NS+1:3*NS+M
        rm = sqrt((x-xDom(i-3*NS,1))^2+(y-yDom(i-3*NS,1))^2);
        phiHat = 0.0;
             if rm>0
                 %phiHat = -(4*deltaT^2)*(besselk(0,rm*sqrt(1/deltaT))+log(rm)+1)-(rm^2*log(rm))*deltaT;
                 phiHat = -(4*deltaT^2)*(besselk(0,rm*sqrt(1/deltaT))+log(rm)+1)-(rm^2*log(rm))*deltaT;
             else
                 phiHat = (4*deltaT^2)*(0.57721+log(sqrt(1/deltaT)/2)-1);
             end
        Tsol=Tsol+coeff(i,1)*phiHat;
    end
    Tsol = Tsol + coeff(3*NS+M+1,1)*psiMONO1(x,y)+coeff(3*NS+M+2,1)*psiMONO2(x,y)+coeff(3*NS+M+3,1)*psiMONO3(x,y);

