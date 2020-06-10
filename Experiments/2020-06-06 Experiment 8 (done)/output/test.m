clc;

x0=[pi;0];     % Initial condition (should satisfy the constraints)
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
Tspan=[0 50];

[t,x]=ode45('system',Tspan,x0,options); % x is the state
t=t';
x=x';

figure(1)
hold off
plot(t(1,:),x(1,:))
hold on
yline(5/3*pi)
yline(1/3*pi)