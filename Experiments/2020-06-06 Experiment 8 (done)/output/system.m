function [xdot] = system(t,x)

% paramters
g = 8.8016;
J = 2.4421e-4;
Km=10.5082;
l=0.0411;
M=0.0762;
tau = 0.398;

u=2.75*sign(sin(1*t));

x1 = x(1); % theta
x2 = x(2); % omega
xdot = [x2; (M*g*l/J)*sin(x1)-(1/tau)*x(2)+(Km/tau)*u]; 

end