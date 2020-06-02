xo = [0 0.2 0.3 0.5 0.3 0.1; 0.005 0.15 0.5 0.55 0.3 0.1; 0.01 0.5 0.3 0.45 0.3 0.1];
t = [0 2000];
global Q1s Q2s ts
x1d = 0.4;
x2d = 0.4;
x3d = 0.4;
xr = [x1d x2d x3d];

A = 0.0154;
S = 5e-5;
S_LM = 0.95;
S_RM = 0.95;
S_LD = 0.1;
S_MD = 0.1;
S_RD = 0.1;
So = 0;
g = 9.81;
Q1 = 0;
Q2 = 0;
Q1s = [];
Q2s = [];
ts = [];

P = [A S S_LM S_RM S_LD S_MD S_RD So g Q1 Q2];

sol1 = [];
sol2 = [];
sol3 = [];

l = size(xo,2);

for i=1:l
    eo = xo(:,i)-xr';
    [t,x] = ode45(@ttse, t,eo,[],P,xr);
    sol1 = [sol1 x(:,1)];
    sol2 = [sol2 x(:,2)];
    sol3 = [sol3 x(:,3)];
end

figure(1)
plot(t,sol1(:,1),t,sol2(:,1),t,sol3(:,1));
legend('e_1','e_2','e_3');
xlabel('time (s)');
ylabel('error (m)');

% figure(2)
% plot(ts,Q1s,ts,Q2s);
% legend('Q_1','Q_2');
% xlabel('time (s)');
% ylabel('flow rate m^3/s');

% figure(3)
% plot(t,sol1(:,1)+x1d,t,sol2(:,1)+x2d,t,sol3(:,1)+x3d);
% legend('T_1','T_2','T_3');
% xlabel('time (s)');
% ylabel('height (m)');
%% 
Q = 1.02e-4;
F=Q;
S = 5e-5;
O = 1.1;
h = (F/(O*S))^2/(2*9.81);


%% 

function dxdt = tts(t,x,P)

A = P(1); S = P(2); S_LM = P(3); S_RM = P(4); S_LD = P(5); S_MD = P(6); S_RD = P(7); So = P(8); g = P(9); Q1 = P(10); Q2 = P(11);

dxdt = [Q1/A-S_LD*S/A*sqrt(2*g*x(1))-S_LM*S/A*sqrt(2*g/abs(x(1)-x(3)))*(x(1)-x(3));
        Q2/A-S_RD*S/A*sqrt(2*g*x(2))-So*S/A*sqrt(2*g*x(2))-S_RM*S/A*sqrt(2*g/abs(x(2)-x(3)))*(x(2)-x(3));
        -S_MD*S/A*sqrt(2*g*x(3))+S_LM*S/A*sqrt(2*g/abs(x(1)-x(3)))*(x(1)-x(3))+S_RM*S/A*sqrt(2*g/abs(x(2)-x(3)))*(x(2)-x(3))];

end

function dxdt = ttse(t,x,P,xr)

A = P(1); S = P(2); S_LM = P(3); S_RM = P(4); S_LD = P(5); S_MD = P(6); S_RD = P(7); So = P(8); g = P(9); Q1 = P(10); Q2 = P(11);

Q1 = min(1.02e-4,A*(S_LD*S/A*sqrt(2*g*(x(1)+xr(1)))+S_LM*S/A*sqrt(2*g/abs((x(1)+xr(1))-(x(3)+xr(3))))*((x(1)+xr(1))-(x(3)+xr(3)))-x(1)));
Q2 = min(1.02e-4,A*(S_RD*S/A*sqrt(2*g*(x(2)+xr(2)))+So*S/A*sqrt(2*g*(x(2)+xr(2)))+S_RM*S/A*sqrt(2*g/abs((x(2)+xr(2))-(x(3)+xr(3))))*((x(2)+xr(2))-(x(3)+xr(3)))-x(2)));
if Q1<0
    Q1 = 0;
end
if Q2<0
    Q2 = 0;
end
y = store(Q1,Q2,t);
% figure(1)
% plot(t,Q1,'o');
% hold on
dxdt = [Q1/A-S_LD*S/A*sqrt(2*g*(x(1)+xr(1)))-S_LM*S/A*sqrt(2*g/abs((x(1)+xr(1))-(x(3)+xr(3))))*((x(1)+xr(1))-(x(3)+xr(3)));
        Q2/A-S_RD*S/A*sqrt(2*g*(x(2)+xr(2)))-So*S/A*sqrt(2*g*(x(2)+xr(2)))-S_RM*S/A*sqrt(2*g/abs((x(2)+xr(2))-(x(3)+xr(3))))*((x(2)+xr(2))-(x(3)+xr(3)));
        -S_MD*S/A*sqrt(2*g*(x(3)+xr(3)))+S_LM*S/A*sqrt(2*g/abs((x(1)+xr(1))-(x(3)+xr(3))))*((x(1)+xr(1))-(x(3)+xr(3)))+S_RM*S/A*sqrt(2*g/abs((x(2)+xr(2))-(x(3)+xr(3))))*((x(2)+xr(2))-(x(3)+xr(3)))];

end

function y = store(Q1,Q2,t)
global Q1s Q2s ts
Q1s = [Q1s Q1];
Q2s = [Q2s Q2];
ts = [ts t];
y = [Q1s;Q2s];
end