function [xdot] = BaseSys(t,x)


A = 0.0154;
g = 9.81;
c=15.2;
S=5e-5;

x1=x(1); %tank1
x2=x(2); %tank3

Q = c*S*sign(x1-x2)*sqrt(2*g*abs(x1-x2));

x1dot = -(1/A)*Q;
x2dot = (1/A)*Q;

xdot=[x1dot;x2dot];
end