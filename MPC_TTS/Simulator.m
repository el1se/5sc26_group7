clear all; close all; clc
%% constants
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S);     % [m]
Atank = 0.0154;         % [m^2]
Hmax = 0.61;            % [m]
Qmax = 0.1;             % [l/s]

L = 0.05;               % [m]
mu = 1.002;             % [Ns/m^2]
rho = 1000;             % [kg/m^3]
g = 9.81;               % [m/s^2]

% system constants (ID or whatever)
v4 = 1;                 % [s]?
v5 = 1;                 % [s]?
% other things
ksim = 250;
Ts = 1;
N = 30;
%% model
BmCT = [0.001/Atank 0 0 0;
        0 0.001/Atank 0 0;
        0 0 0 0 ;
        0 0 0.01*v4*Dvalve 0;
        0 0 0 0.01*v5*Dvalve];
frac = 128*mu*L*Atank;
dfdx = @(x) [   
    -pi*x(4)^4*rho*g/frac   0                       pi*x(4)^4*rho*g/frac                -pi*4*x(4)^3*rho*g*(x(1)-x(3))/frac     0;
    0                       -pi*x(5)^4*rho*g/frac   pi*x(5)^4*rho*g/frac                0                                       pi*4*x(5)^3*rho*g*(x(3)-x(2))/frac;
    pi*x(4)^4*rho*g/frac    pi*x(5)^4*rho*g/frac    -pi*(x(5)^4+x(4)^4)*rho*g/frac      4*pi*x(4)^3*rho*(x(1)-x(3))/frac        -4*pi*x(5)^3*rho*(x(3)-x(2))/frac;
    0                       0                       0                                   -v4                                     0;
    0                       0                       0                                   0                                       -v5
    ];

 Cm = eye(5);

m = size(BmCT,2);
n = size(BmCT,1);
q = size(Cm,1);

Dm = zeros(q,m);
%% constraints
deltaumin = [-Qmax; -Qmax; -100; -100];
deltaumax = -deltaumin;
umin  = zeros(4,1);
umax  = [Qmax; Qmax; 100; 100];
ymin  = zeros(q,1);
ymax  = [zeros(3,1) + Hmax; Dvalve; Dvalve];

%% reference How to design reference for valves?
r = [zeros(q-2,ksim+N)+0.3; zeros(2,ksim+N)];
%% Q & R
Q = [0.001*eye(5) zeros(5,5);
    zeros(5,5) eye(5)];
Qm = eye(5);
R = eye(4);

%% state preparation
X = zeros(n,ksim);
Xaug = zeros(n+q,ksim);
u = zeros(m,ksim);
u0 = zeros(m,1);

X(:,1) = [0; 0; 0;Dvalve;Dvalve];
Xaug(:,1) = [zeros(n,1);X(1:q,1)];
options_qp =  optimoptions('quadprog','Display','off');

for k = 1:ksim
    AmCT = dfdx(X(:,k));
    sysCT = ss(AmCT,BmCT,Cm,Dm);
    sysDT = c2d(sysCT,Ts);
    Am = sysDT.A;
    Bm = sysDT.B;
    
    
    A = [Am zeros(q,n)';
         Cm*Am eye(q)];
    B = [Bm;Cm*Bm];
    C = [zeros(q,n) eye(q)];

    % integral model
    [~,P,~] = dlqr(Am,Bm,Qm,R);
    [ Psi, Omega ] = QRPN2PsiOmega(C*Q*C',R,P,N );    
    [Phi, Gamma] = ABCN2FPhi(A,B,C,N);
    G = 2*(Psi+Gamma'*Omega*Gamma);
    F = 2*Gamma'*Omega;
   

     Rk = [];
    for i = 1:N
        Rk = [Rk; r(:,i)];
    end

    
    if k >1
        [D,M,E,c] = DMEcIntegral(deltaumin,deltaumax,umin,umax,ymin,ymax,u(:,k-1),N);
    else
        [D,M,E,c] = DMEcIntegral(deltaumin,deltaumax,umin,umax,ymin,ymax,zeros(m,1),N);
    end
    
    L = M*Gamma+E;
    W = -D*C-M*Phi;
    [u_qp,fval,exitflag] = quadprog((G+G')/2,(F*Phi*Xaug(:,k)-F*Rk),L,c+W*(Xaug(:,k)),[],[],[],[],[],options_qp);
    ustark = -inv(G)*(F*Phi*Xaug(:,k)-F*Rk);
    deltaustar0k = [eye(m) zeros(m,m*(N-1))]*u_qp;
    if k > 1
        u(:,k) = u(:,k-1) + deltaustar0k;
    else
        u(:,k) = u0 + deltaustar0k;
    end
    
    Xaug(:,k+1) = A*Xaug(:,k) + B*deltaustar0k;
    X(:,k+1) = Xaug(n+1:end,k+1);
end

%%
t = 0:Ts:(size(X,2)-1)*Ts;
figure(1); clf;
subplot(2,2,1)
hold on;
plot(t,X(1,:));plot(t,X(2,:));plot(t,X(3,:));
legend('Water level tank 1','Water level tank 2','Water level tank 3');
axis([-0.1 t(end) -0.01 0.7])
subplot(2,2,2);
hold on;
plot(t,X(4,:));plot(t,X(5,:));
legend('Valve diameter 1','Valve diameter 2')
axis([-0.1 t(end) 0 0.015])
subplot(2,2,3);
hold on;
plot(t(1:end-1),u(1,:));plot(t(1:end-1),u(2,:));
legend('Pump flow 1','Pump flow 2')
axis([-0.1 t(end) 0 0.12])
subplot(2,2,4);
hold on;
plot(t(1:end-1),u(3,:));plot(t(1:end-1),u(4,:));
legend('valve input 1','valve input 2')
axis([-0.1 t(end) -105 105])
 