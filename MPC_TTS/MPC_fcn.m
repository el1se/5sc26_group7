function [u, P] = MPC_fcn(Y, uprev, Pprev, YPrev)
% outputs 4d vector u containing inputs for TTS  (format pump1 pump2 valve LM valve RM)
% outputs 3x3 matrix P, containing solution to DARE, can be used in
% next itration
% Y: 5d input state vector of TTS (format water level tank 1, tank 2, tank 3, connecting pipe diameter LM)
% uprev: previous input sent to plant (4D vector, same format as u), for augmented model necessary
% Pprev: previous solution to DARE, might be used if no stabilizing
% solution can be calculated
% Yprev: previous 5d state vector of TTS, same format as Y
%% constants
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S);     % [m]
Atank = 0.0154;         % [m^2]
g = 9.81;               % [m/s^2]
L = 0.1;                % [m]
mu = 1.002;             % [Ns/m^2]
rho = 1000;             % [kg/m^3]
% system constants (ID or whatever)
v4 = 0.5;                 % [s]?
v5 = 0.5;                 % [s]?

% Constraints 
Hmax = 0.61;            % [m]
Qmax = 0.1;             % [l/s]

% other things
Ts = 2;                 % TODO: check if need be due to .mat file?
N = 25;
%% model
BmCT = [0.001/Atank 0 0 0;
        0 0.001/Atank 0 0;
        0 0 0 0 ;
        0 0 0.01*v4*Dvalve 0;
        0 0 0 0.01*v5*Dvalve];
Cm = eye(5);


m = size(BmCT,2);
n = size(BmCT,1);
q = size(Cm,1);

Dm = zeros(q,m);
%% constraints
deltaumin = [-0.1*Qmax; -0.1*Qmax; -100; -100];
deltaumax = -deltaumin;
umin  = zeros(4,1);
umax  = [Qmax; Qmax; 100; 100];
ymin  = zeros(q,1);
ymax  = [zeros(3,1) + Hmax; Dvalve; Dvalve];

%% reference How to design reference for valves?
r = [zeros(q-2,1)+0.3; zeros(2,1)];
%% Q & R
Q = [eye(3) zeros(3,2);
        zeros(2,5)];
Qm = [  eye(3)      zeros(3,2);
        zeros(2,5)                              ];
R = [eye(2) zeros(2,2);
    zeros(2,2) 0.00001*eye(2)];

%% calculations

% Xaug = [zeros(n,1);Y(1:q,1)];

options_qp =  optimoptions('quadprog','Display','off');

% Pprev = [15.9081147155710,2.96062907928527e-12,0.000675171321078397,-0.874073498351786,4.34237622311551e-05;2.96062907928527e-12,15.9081147155711,0.000675171327952177,4.34237625058848e-05,-0.874073498351258;0.000675171321078397,0.000675171327952177,155486.683787255,9999.99731184461,9999.99731184459;-0.874073498351786,4.34237625058848e-05,9999.99731184461,643.201217531435,643.148089872505;4.34237622311551e-05,-0.874073498351258,9999.99731184459,643.148089872505,643.201217531434];

AmCT = freeFallLinearizationA(Y);
sysCT = ss(AmCT,BmCT,Cm,Dm);
sysDT = c2d(sysCT,Ts);
Am = sysDT.A;
Bm = sysDT.B;


A = [Am zeros(q,n)';
     Cm*Am eye(q)];
B = [Bm;Cm*Bm];
C = [zeros(q,n) eye(q)];

Xaug = [    Y-YPrev;
            Y           ];

% integral model
if rank(ctrb(Am,Bm)) == 5
    [~,P,~] = dlqr(Am,Bm,Qm,R);
else
    P = Pprev;
end

[ Psi, Omega ] = QRPN2PsiOmega(Q,R,P,N);    
[Phi, Gamma] = ABCN2FPhi(A,B,C,N);
G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega;


 Rk = [];
for i = 1:N
    Rk = [Rk; r];
end

[D,M,E,c] = DMEcIntegral(deltaumin,deltaumax,umin,umax,ymin,ymax,uprev,N);


Lmpc = M*Gamma+E;
W = -D*C-M*Phi;
[u_qp,~,~] = quadprog((G+G')/2,(F*Phi*Xaug-F*Rk),Lmpc,c+W*Xaug,[],[],[],[],[],options_qp);

deltaustar0k = [eye(m) zeros(m,m*(N-1))]*u_qp;
u = uprev + deltaustar0k;
%%-------------------------------------------------------------------------
% functions
function A = freeFallLinearizationA(x)
tol = 0;
if ((x(2)-x(3)) <= tol) && ((x(1)-x(3)) > tol)
    A = [   
        -S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank    0    S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank     -sign(x(1)-x(3))*S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank   0;
        zeros(1,5);
        S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank     0    -S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank    sign(x(1)-x(3))*S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank    0;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
elseif ((x(1)-x(3)) <= tol) && ((x(2)-x(3)) > tol)
    A = [   
        zeros(1,5);
        0                       -S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank    S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank     0    -sign(x(2)-x(3))*S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank     -S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank    0    sign(x(2)-x(3))*S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
elseif ((x(2)-x(3)) <= tol) && ((x(1)-x(3)) <= tol)
     A = [   
        zeros(3,5);
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
else
    A = [   
        -S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank     0                                                   sign(x(1)-x(3))*S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank                                                      -S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank   0;
        0                                                   -S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank     sign(x(2)-x(3))*S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank                                                      0                                       -S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank      S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank      -sign(x(1)-x(3))*S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank-sign(x(2)-x(3))*S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank      S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank    S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
end
end

%%-------------------------------------------------------------------------
function [D,M,E,c] = DMEcIntegral(deltaumin,deltaumax,umin,umax,ymin,ymax,ukmin1,N)
m = size(umin,1);
q = size(ymin,1);

bi = [-deltaumin; deltaumax; -umin+ukmin1;umax-ukmin1;-ymin;ymax];
bN = [zeros(4*m,1);-ymin; ymax];
Lb = length(bi);

c = zeros((4*m+2*q)*(N+1),1);
for i = 1:N
    c((i-1)*Lb+1:i*Lb) = bi;
end
c(N*Lb+1:end,:) = bN;


M0 = [zeros(4*m,q); -eye(q); eye(q)];

D = zeros((4*m+2*q)*(N+1),q);
D(1:size(M0,1),:)  = M0;

M = zeros((4*m+2*q)*(N+1),q*N);
for i = 1:N
    M(i*Lb+1:(i+1)*Lb,(i-1)*q+1:i*q) = M0;
end

Ei = [-eye(m); eye(m);-eye(m); eye(m);zeros(2*q,m)];
E0ff = [zeros(2*m,m);-eye(m); eye(m);zeros(2*q,m)];

E = zeros((4*m+2*q)*(N+1),m*N);
for i = 1:N
    E((i-1)*Lb+1:(i)*Lb,(i-1)*m+1:i*m) = Ei;
    if i > 1
        for j = 1:i-1
            E((i-1)*Lb+1:(i)*Lb,(j-1)*m+1:j*m) = E0ff;
        end
    end
end

end
%%-------------------------------------------------------------------------
function [ Psi, Omega ] = QRPN2PsiOmega( Q,R,P,N )
% Omega
if (N==1)
    Omega=P;
else
    Omega=Q;
    for j=2:(N-1)
        Omega=blkdiag(Omega,Q);
    end
    Omega=blkdiag(Omega,P);
%     omega=blkdiag(omega,Q);
end

% Psi
Psi = R;
if N>1
    for j = 2:N
        Psi = blkdiag(Psi,R);
    end
end
end
%%-------------------------------------------------------------------------
function [F, Gamma ] = ABCN2FPhi(A,B,C,N)
%ABN2PhiGamma Summary of this function goes here
%   Detailed explanation goes here

n = size(B,1);
m = size(B,2);
q = size(C,1);


F = [];
for i = 1:N
   F = [F;C*A^i];
end

% phi
Gamma = zeros(q*N,m*N);
for i = 1:N
    for j = 1:i
        Gamma((i-1)*q+1:q*i,(j-1)*m+1:j*m) = C*A^(i-j)*B;
    end
end
end
end

