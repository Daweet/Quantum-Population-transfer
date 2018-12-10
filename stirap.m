function dydt= stirap(t,y)
%solve the STIRAP example
dydt =zeros (size(y));
%parameters - and constants
T=0.413e12;
omega0=20/T;  tau=1.2*T;
%omegap=omega0*0.707*exp(-((t-tau)^2)/(2*T^2)); 
%omegas=(omega0*exp(-((t+tau)^2)/(T^2)))+(omega0*0.707*exp(-((t-tau)^2)/(2*T^2))); 
omegap=omega0*exp(-((t-tau)^2)/(2*T^2)); 
omegas=omega0*exp(-((t+tau)^2)/(2*T^2));
A=y(1);
B=y(2);
C=y(3);
%Evaluate the RHS expression
dydt(1)=i*0.5*omegap*B;
dydt(2)=i*0.5*(omegap*A + omegas*C);
dydt(3)=i*0.5*omegas*B;
