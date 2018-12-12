% Using MATLAB built in partial derivative solver
%[T,Y] = ode45(@hmltn,t,a,[],delta,height,width);
[t,y] = ode45('fstirap1',[-4.5e12 9.5e12],[0 0 0 0 0 0 0 -1 -1/(sqrt(3))]);
%calls stirap1 function
% printing pulses
%


%f1 = figure('visible','off');  %---------------------------------
figure(1)
%subplot(2,1,1);      %----------------------------

plot(t,y.*conj(y))
title('Populations','fontsize',14);
xlabel('time','fontsize',14);
ylabel('Populations','fontsize',14)
legend('\rho_{1}','\rho_{2}','\rho_{3}','\rho_{4}','\rho_{5}','\rho_{6}','\rho_{7}','\rho_{8}','\rho_{9}');



%f1 = figure('visible','off');  %---------------------------------
figure(2)
%subplot(2,1,1);      %----------------------------

plot(t,y)
title('Populations','fontsize',14);
xlabel('time','fontsize',14);
ylabel('Populations','fontsize',14)
legend('\rho_{1}','\rho_{2}','\rho_{3}','\rho_{4}','\rho_{5}','\rho_{6}','\rho_{7}','\rho_{8}','\rho_{9}');
