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
v4 = 3;                 % [s]?
v5 = 3;                 % [s]?
% other things
ksim = 500;
Ts = 0.05;
N = 10;

%% constraints
deltaumin = [-Qmax; -Qmax; -100; -100];
deltaumax = -deltaumin;
umin  = zeros(4,1);
umax  = [Qmax; Qmax; 100; 100];
ymin  = zeros(3,1);
ymax  = zeros(3,1) + Hmax;
%% model
BmCT = [0.001/Atank 0 0 0;
        0 0.001/Atank 0 0;
        0 0 0 0 ;
        0 0 -0.01*v4/Dvalve 0;
        0 0 0 -0.01*v5/Dvalve];
dfdx = @(x) [   -pi*x(4)^4*rho*g/(128*mu*L*Atank) 0 pi*x(4)^4*rho*g/(128*mu*L*Atank) -pi*4*x(4)^3*rho*g*(x(1)-x(3))/(128*mu*L*Atank) 0;
                0 -pi*x(5)^4*rho*g/(128*mu*L*Atank) pi*x(5)^4*rho*g/(128*mu*L*Atank) 0 -pi*4*x(5)^3*rho*g*(x(3)-x(2))/(128*mu*L*Atank);
                pi*x(4)^4*rho*g/(128*mu*L*Atank) pi*x(5)^4*rho*g/(128*mu*L*Atank) -pi*x(5)^4*rho*g/(128*mu*L*Atank) 4*pi*x(4)^3*rho*(x(1)-x(3))/(128*mu*L*Atank) -4*pi*x(5)^3*rho*(x(3)-x(2))/(128*mu*L*Atank);
                0 0 0 v4 0;
                0 0 0 0 v5];
    
    
    Cm = [   1 0 0 0 0;
        0 1 0 0 0;
        0 0 1 0 0];

m = size(BmCT,2);
n = size(BmCT,1);
q = size(Cm,1);

Dm = zeros(q,m);
%% reference
r = zeros(q,ksim)+0.3;
%% Q & R
Q = [0.01*eye(5) zeros(5,3);
    zeros(3,5) eye(3)];
R = eye(4);

%% state preparation
X = zeros(n,ksim);
Xaug = zeros(n+q,ksim);
u = zeros(m,ksim);
u0 = zeros(4,1);

X(:,1) = [0.2; 0.1; 0;0.5*Dvalve;0.2*Dvalve];
Xaug(:,1) = [zeros(5,1);X(1:3,1)];
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


    [~,Pm,~] = dlqr(A,B,Q,R);
    P = pinv(C)'*Pm*pinv(C);
    [ Psi, Omega ] = QRPN2PsiOmega(C*Q*C',R,P,N );    
    [Phi, Gamma] = ABCN2FPhi(A,B,C,N);
%     Phiy = C*Phi*C';
    G = 2*(Psi+Gamma'*Omega*Gamma);
    F = 2*Gamma'*Omega;
    Rk = [];
    for i = 1:q
        Rk = [Rk; r(i,k+1:k+N)'];
    end
    
    if k >1
        [D,M,E,c] = DMEcIntegral(deltaumin,deltaumax,umin,umax,ymin,ymax,u(:,k-1),N);
    else
        [D,M,E,c] = DMEcIntegral(deltaumin,deltaumax,umin,umax,ymin,ymax,zeros(m,1),N);
    end
    
    L = M*Gamma+E;
    W = -D*C-M*Phi;
    [u_qp,fval,exitflag] = quadprog((G+G')/2,F*Phi*Xaug(:,k),L,c+W*Xaug(:,k),[],[],[],[],[],options_qp);
    deltaustar0k = u_qp(1:m,:);
    if k > 1
        u(:,k) = u(:,k-1) + deltaustar0k;
    else
        u(:,k) = u0 + deltaustar0k;
    end
    
    Xaug(:,k+1) = A*Xaug(:,k) + B*deltaustar0k;
    X(:,k+1) = Am* X(:,k)+Bm*u(:,k);
    
    
end
 