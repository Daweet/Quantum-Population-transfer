	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%                                    %
	%       FILE NAME : fstirap_prop      %
	%                                    %  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	

function dydt = fstirap1(t,y)

%function stirap -- before the changes for the matlab web server
% programe for analysing population transfer in 3 level lambda STIRAP systems
% 
%
%                                                             ----- |2>
%                                                              / \
%                                                             /   \
%                                                     |1>   ---   --- |3>
% The program propagats the wavefunction in energy eigenstate representation,
% under the unitary Hamiltonian evolution (No dissipation).
% A STIRAP mechanism where the prope proceedes the pump (anti intuitive), is 
% demonstrated.





dydt=zeros(size(y));  % Initial population distribution


%-----parameters and constants -------------

T=0.413e12;
omega0=0.1*(20/T);% Laser intensity
delta=1.376e-12;     % Field  detuning
tau=0.7*T;
%omegap=omega0*sin((pi/180)*45)*exp(-((t-tau)^2)/(T^2));  
%omegas=omega0*exp(-((t+tau)^2)/(T^2))+omega0*cos((pi/180)*45)*exp(-((t-tau)^2)/(T^2)); 
omegap=omega0*sin((pi/4))*exp(-((t-tau)^2)/(T^2));  
omegas=omega0*exp(-((t+tau)^2)/(T^2))+omega0*cos((pi/4))*exp(-((t-tau)^2)/(T^2));

A=y(1);
B=y(2);
C=y(3);
D=y(4);
E=y(5);
F=y(6);
G=y(7);
H=y(8);
I=y(9);

%Evaluating 
% O.D.E 
dydt(1)=0;
dydt(2)=(delta*E+0.5*omegas*G);
dydt(3)=(-delta*F-0.5*omegap*G);
dydt(4)=(0.5*omegas*E-0.5*omegap*F);
dydt(5)=(-delta*B-0.5*omegas*D+omegap*H);
dydt(6)=(delta*C+0.5*omegap*D-omegas*H+ sqrt(3)*omegas*I);
dydt(7)=(-0.5*omegas*B+0.5*omegap*C);
dydt(8)=(-omegap*E+0.5*omegas*F);
dydt(9)=(-sqrt(3)*omegas*F);







