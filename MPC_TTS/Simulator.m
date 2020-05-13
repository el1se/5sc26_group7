clear all; close all; clc
%% constants
global Atank v4 v5 g Dvalve S rho mu L
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S);     % [m]
Atank = 0.0154;         % [m^2]
g = 9.81;               % [m/s^2]
L = 0.1;                % [m]
mu = 1.002;             % [Ns/m^2]
rho = 1000;             % [kg/m^3]
% system constants (ID or whatever)
v4 = 1;                 % [s]?
v5 = 1;                 % [s]?

% Constraints 
Hmax = 0.61;            % [m]
Qmax = 0.1;             % [l/s]

% other things
ksim = 200;
Ts = 1;
N = 10;
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
deltaumin = [-0.1*Qmax; -0.1*Qmax; -100*10000000; -100*10000000];
deltaumax = -deltaumin;
umin  = zeros(4,1);
umax  = [Qmax; Qmax; 100; 100];
ymin  = zeros(q,1);
ymax  = [zeros(3,1) + Hmax; Dvalve; Dvalve];

%% reference How to design reference for valves?
r = [zeros(q-2,ksim+N)+0.3; zeros(2,ksim+N)];
%% Q & R
Q = [eye(3) zeros(3,2);
        zeros(2,5)];
Qm = [  eye(3)      zeros(3,2);
        zeros(2,5)                              ];
R = eye(4);

%% state preparation
X = zeros(n,ksim);
Xaug = zeros(n+q,ksim);
u = zeros(m,ksim);


X(:,1) = [0.01; 0.01; 0;0*Dvalve;0*Dvalve];
Xaug(:,1) = [zeros(n,1);X(1:q,1)];
options_qp =  optimoptions('quadprog','Display','off');

u0 = [0;0;X(4,1)/Dvalve*100;X(5,1)/Dvalve*100];

for k = 1:ksim
    if ((abs(X(3,k)-X(1,k)) < 1e-5) || (abs(X(3,k)-X(2,k)) < 1e-5))
        error('ERROR: liquid level of tank 1 or 2 to close to liquid level of tank 3, cannot evaluate jacobian')
    else
        AmCT = freeFallLinearizationA(X(:,k));
    end
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
%      P( abs(P) < 1e-9 ) = 0;
    [ Psi, Omega ] = QRPN2PsiOmega(Q,R,P,N );    
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
        [D,M,E,c] = DMEcIntegral(deltaumin,deltaumax,umin,umax,ymin,ymax,u0,N);
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
plot(t,X(1,:),'LineWidth',1.5);plot(t,X(2,:),'LineWidth',1.5);plot(t,X(3,:),'LineWidth',1.5);
legend('Tank 1','Tank 2','Tank 3');
axis([-0.1 t(end) -0.01 0.7])
xlabel('Time [s]');
ylabel('Water level [m]');
subplot(2,2,2);
hold on;
plot(t,X(4,:),'LineWidth',1.5);plot(t,X(5,:),'LineWidth',1.5);
legend('Valve 1','Valve 2')
axis([-0.1 t(end) 0 0.019])
xlabel('Time [s]');
ylabel('Diameter connecting pipe [m]');
subplot(2,2,3);
hold on;
plot(t(1:end-1),u(1,:),'LineWidth',1.5);plot(t(1:end-1),u(2,:),'LineWidth',1.5);
legend('Pump 1','Pump 2')
axis([-0.1 t(end) 0 0.12])
xlabel('Time [s]');
ylabel('Pump volume flow input [l/s]');
subplot(2,2,4);
hold on;
plot(t(1:end-1),u(3,:),'LineWidth',1.5);plot(t(1:end-1),u(4,:),'LineWidth',1.5);
legend('Valve 1','Valve 2')
xlabel('Time [s]');
ylabel('Valve input [%]');
axis([-0.1 t(end) -105 105])
 