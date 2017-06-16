function qnm1=Force(x,y,Tin,coeff,xDom,yDom,NS,M,s,phiMONO1,phiMONO2,phiMONO3)

fnm1=0;
for i=1:M
    phiTild = (sqrt((x-xDom(i,1))^2+(y-yDom(i,1))^2))^2*log(sqrt((x-xDom(i,1))^2+(y-yDom(i,1))^2));
    fnm1=fnm1+coeff(i,s)*phiTild;
end
fnm1 = fnm1 + coeff(3*NS+M+1,s)*phiMONO1(x,y)+coeff(3*NS+M+2,s)*phiMONO2(x,y)+coeff(3*NS+M+3,s)*phiMONO3(x,y);

qnm1=-1*fnm1-Tin;