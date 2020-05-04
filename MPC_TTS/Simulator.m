clear all; close all; clc
%% constants
S = 5e-5;       % [m^2]
Dvalve = sqrt(4*S);  % [m]
Atank = 0.0154;     % [m^2]
Hmax = 0.61;    % [m]
Qmax = 0.1;     % [l/s]
L = 0.05;       % [m]
mu = 1.002;     % [Ns/m^2]
rho = 1000;     % [kg/m^3]
g = 9.81;       % [m/s^2]
% system constants (ID or whatever)
v4 = 3;         % [s]?
v5 = 3;         % [s]?
% other things
ksim = 500;
Ts = 0.5;
N = 10;

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
    

% Cm = [   1 0 0 0 0;
%         0 1 0 0 0;
%         0 0 1 0 0;
%         0 0 0 1/(5e-5) 0;
%         0 0 0 0 1/(5e-5)];
    
    Cm = [   1 0 0 0 0;
        0 1 0 0 0;
        0 0 1 0 0];

m = size(BmCT,2);
n = size(BmCT,1);
q = size(Cm,1);

Dm = zeros(q,m);
%% reference
r = zeros(q,ksim)+0.3;
%% QPR
Q = [eye(3) zeros(3,2);
    zeros(2,3) eye(2)];
% Q = [eye(3) zeros(3,2);
%     zeros(2,5)];
R = eye(4);
%% constraints
xmin = zeros(5,1);
xmax = [Hmax Hmax Hmax Dvalve Dvalve]';
umin = zeros(4,1);
umax = [Qmax Qmax 100 100];
%% state preparation
X = zeros(n,ksim);
u = zeros(m,ksim);

X(:,1) = [0.2; 0.1; 0.3;0.5*Dvalve;0.2*Dvalve];



for k = 1:ksim
    AmCT = dfdx(X(:,k));
    sysCT = ss(AmCT,BmCT,Cm,Dm);
    sysDT = c2d(sysCT,Ts);
%     sysDTMin = ss(sysDT,'minimal');
    Am = sysDT.A;
    Bm = sysDT.B;
%     [Am,Bm,Cm,Dm]=ssdata(sysDTMin);
    A = [Am zeros(q,n)';
         Cm*Am eye(q)];
    B = [Bm;Cm*Bm];
    C = [zeros(q,n) eye(q)];
    Dtilde = [Cm*Am eye(q)];
    OL = [];
    for i = 1:N
        OL = [OL; Dtilde*A^(i-1)];
    end
    [K,P,~] = dlqr(Am, Bm,Q,R);
    [ Psi, Omega ] = QRPN2PsiOmega(Q,R,P,N )
    sysAug = (ss(A,B,C,zeros(size(C,1),size(B,2)),Ts));
    
%     sysAugMin = ss(sysAug,'minimal');
%     Aa = sysAugMin.A; Ba = sysAugMin.B; Ca = sysAugMin.C; 

    [Aa,Ba,Ca,Da]=ssdata(sysAug);
    
    [Phi, Gamma] = ABCN2FPhi(Aa,Ba,Ca,N);
    
    G = 2*(Psi+Gamma'*Omega*Gamma);
    F = 2*Gamma'*Omega;
    Rk = r(k+1:k+N)';
    deltaUstar = -inv(G)*(F*phi*X(:,k)-F*Rk)
%     deltaU = inv(phi'*phi+Rbar)*(phi'*Rsbar*Rk-phi'*F*X(:,k))
end
 