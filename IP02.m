
xo = [0.2 0.3;0.15 0.4 ;0.5 0.3;0.9 0.7;0.9 0.7] ;
t = [0 100];

A = 0.0154;
S = 5e-5;
S_LD = 0.3;
S_MD = 0.2;
S_RD = 0.1;
So = 0;
g = 9.81;
Q1 = 0;
Q2 = 0;

P = [A S S_LD S_MD S_RD So g Q1 Q2];

x1d = 0.2;
x2d = 0.3;
x3d = 0.35;
xr = [x1d x2d x3d];

sol1 = [];
sol2 = [];
sol3 = [];
sol4 = [];
sol5 = [];
l = size(xo,2);

for i =1:l
    eo = xo(:,i)-[xr'; xo(4:5,i)];
    [t,x] = ode45(@ttse, t,eo,[],P,xr);
    sol1 = [sol1 x(:,1)];
    sol2 = [sol2 x(:,2)];
    sol3 = [sol3 x(:,3)];
    sol4 = [sol4 x(:,4)];
    sol5 = [sol5 x(:,5)];
end

figure(1)
plot(t,sol1(:,1),t,sol2(:,1),t,sol3(:,1));
legend('e_1','e_2','e_3');
xlabel('time (s)');
ylabel('error (m)');

figure(2)
plot(t,sol4(:,1),t,sol5(:,1));
legend('S_LM','S_RM');
xlabel('time (s)');
ylabel('position');
%% 
h= 0.45;
S = 5e-5;
F = 0.8*S*sqrt(2*9.81*h);
Q = 1.02e-4;
r = F/Q;

%% 

% function dxdt = tts(t,x,P)
% 
% A = P(1); S = P(2); S_LM = P(3); S_RM = P(4); S_LD = P(5); S_MD = P(6); S_RD = P(7); So = P(8); g = P(9); Q1 = P(10); Q2 = P(11);
% 
% dxdt = [Q1/A-S_LD*S/A*sqrt(2*g*x(1))-S_LM*S/A*sqrt(2*g/abs(x(1)-x(3)))*(x(1)-x(3));
%         Q2/A-S_RD*S/A*sqrt(2*g*x(2))-So*S/A*sqrt(2*g*x(2))-S_RM*S/A*sqrt(2*g/abs(x(2)-x(3)))*(x(2)-x(3));
%         -S_MD*S/A*sqrt(2*g*x(3))+S_LM*S/A*sqrt(2*g/abs(x(1)-x(3)))*(x(1)-x(3))+S_RM*S/A*sqrt(2*g/abs(x(2)-x(3)))*(x(2)-x(3))];
% 
% end

function dxdt = ttse(t,x,P,xr)

A = P(1); S = P(2); S_LD = P(3); S_MD = P(4); S_RD = P(5); So = P(6); g = P(7); Q1 = P(8); Q2 = P(9);

Q1 = min(1.02e-4,A*(S_LD*S/A*sqrt(2*g*(x(1)+xr(1)))+x(4)*S/A*sqrt(2*g/abs((x(1)+xr(1))-(x(3)+xr(3))))*((x(1)+xr(1))-(x(3)+xr(3)))-x(1)));
Q2 = min(1.02e-4,A*(S_RD*S/A*sqrt(2*g*(x(2)+xr(2)))+So*S/A*sqrt(2*g*(x(2)+xr(2)))+x(5)*S/A*sqrt(2*g/abs((x(2)+xr(2))-(x(3)+xr(3))))*((x(2)+xr(2))-(x(3)+xr(3)))-x(2)));
dxdt = [Q1/A-S_LD*S/A*sqrt(2*g*(x(1)+xr(1)))-x(4)*S/A*sqrt(2*g/abs((x(1)+xr(1))-(x(3)+xr(3))))*((x(1)+xr(1))-(x(3)+xr(3)));
        Q2/A-S_RD*S/A*sqrt(2*g*(x(2)+xr(2)))-So*S/A*sqrt(2*g*(x(2)+xr(2)))-x(5)*S/A*sqrt(2*g/abs((x(2)+xr(2))-(x(3)+xr(3))))*((x(2)+xr(2))-(x(3)+xr(3)));
        -S_MD*S/A*sqrt(2*g*(x(3)+xr(3)))+x(4)*S/A*sqrt(2*g/abs((x(1)+xr(1))-(x(3)+xr(3))))*((x(1)+xr(1))-(x(3)+xr(3)))+x(5)*S/A*sqrt(2*g/abs((x(2)+xr(2))-(x(3)+xr(3))))*((x(2)+xr(2))-(x(3)+xr(3)));
        -(x(1)-xr(2)-0.2);
        -(x(1)-xr(2)-0.2)];
end