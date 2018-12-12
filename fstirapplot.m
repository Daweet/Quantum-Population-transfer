[t,y]=ode45('fstirap',[-4.5e12 9.5e12],[1 0 0]);
figure(1);
plot(t,y.*conj(y))
xlabel('Time','FontSize',14);
ylabel('Population','FontSize',14);
legend('\rho_{11}(t)','\rho_{22}(t)','\rho_{33}(t)')
axis([-4.5e12 4.5e12 0 1]);
figure(2);
for n=1:3
plot(t,y()*conj(y(n)))
xlabel('Time','FontSize',14);
ylabel('Coherences','FontSize',14);
legend('\rho_{12}(t)','\rho_{23}(t)','\rho_{13}(t)')
axis([-4.5e12 9.5e12 -1 1]);
end
