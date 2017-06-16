% Visualization of the numerical solution VS true solution/measurements

x=ksi;
y=eta; 

tr=intPntsMeas1(3:end,25);  
%tr=sin(tau)*sin(ksi)*sin(eta);  
figure
plot(tau(1:701),Tsol(1:701),'b--',tau(1:701),tr(1:701),'green')
%plot(tau,Tsol,'b--',tau,tr,'green')
relEr=norm(Tsol-tr)/norm(tr)  % relative error
legend('approx solution','true solution')

% Comparison to interpolation (for points on the axis) %
%inter=(intPntsMeas1(3:end,18)+intPntsMeas1(3:end,32))/2;
%figure
%plot(tau(201:501),Tsol(201:501),'b--',tau(201:501),tr(201:501),'green',tau(201:501),inter(201:501),'red')
%legend('approx solution','true solution','interpolation')