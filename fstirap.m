function dydt= fstirap(t,y)
%solve the STIRAP example
dydt =zeros (size(y));
%parameters - and constants
T=0.413e12;delta=1.376e-12;
omega0=0.1*(20/T);  tau=0.7*T;
omegap=omega0*sin((pi/4))*exp(-((t-tau)^2)/(T^2)); 
omegas=(omega0*exp(-((t+tau)^2)/(T^2)))+(omega0*cos((pi/4))*exp(-((t-tau)^2)/(T^2))); 
%omegap=omega0*exp(-((t-tau)^2)/(2*T^2)); 
%omegas=omega0*exp(-((t+tau)^2)/(2*T^2));
A=y(1);
B=y(2);
C=y(3);
%Evaluate the RHS expression
dydt(1)=i*0.5*omegap*B;
dydt(2)=i*0.5*(omegap*A +2*delta*B+ omegas*C);
dydt(3)=i*0.5*omegas*B;
