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
v4 = 0.5;                 % [s]?
v5 = 0.5;                 % [s]?

% Constraints 
Hmax = 0.61;            % [m]
Qmax = 0.1;             % [l/s]

% other things
ksim = 300;
Ts = 2;
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
r = [zeros(q-2,ksim+N)+0.3; zeros(2,ksim+N)];
%% Q & R
Q = [eye(3) zeros(3,2);
        zeros(2,5)];
Qm = [  eye(3)      zeros(3,2);
        zeros(2,5)                              ];
% R = eye(4);
R = [eye(2) zeros(2,2);
    zeros(2,2) 0.00001*eye(2)];

%% state preparation
Y = zeros(n,ksim);
Xaug = zeros(n+q,ksim);
u = zeros(m,ksim);


Y(:,1) = [0; 0.28; 0;0*Dvalve;0*Dvalve];
X2 = Y;
Xaug(:,1) = [zeros(n,1);Y(1:q,1)];
options_qp =  optimoptions('quadprog','Display','off');

u0 = [0;0;Y(4,1)/Dvalve*100;Y(5,1)/Dvalve*100];
Pprev = [15.9081147155710,2.96062907928527e-12,0.000675171321078397,-0.874073498351786,4.34237622311551e-05;2.96062907928527e-12,15.9081147155711,0.000675171327952177,4.34237625058848e-05,-0.874073498351258;0.000675171321078397,0.000675171327952177,155486.683787255,9999.99731184461,9999.99731184459;-0.874073498351786,4.34237625058848e-05,9999.99731184461,643.201217531435,643.148089872505;4.34237622311551e-05,-0.874073498351258,9999.99731184459,643.148089872505,643.201217531434];
for k = 1:ksim

    AmCT = freeFallLinearizationA(Y(:,k));
    sysCT = ss(AmCT,BmCT,Cm,Dm);
    sysDT = c2d(sysCT,Ts);
    Am = sysDT.A;
    Bm = sysDT.B;
    
    
    A = [Am zeros(q,n)';
         Cm*Am eye(q)];
    B = [Bm;Cm*Bm];
    C = [zeros(q,n) eye(q)];

%     Lobs = place(A',C',[-0.8 0.8 0.9 -0.9 0.7 -0.7 0.6 -0.6 0.5 -0.5])';
    % integral model
    if rank(ctrb(Am,Bm)) == 5
        [~,P,~] = dlqr(Am,Bm,Qm,R);
    else
        P = Pprev;
    end
    Pprev = P;
    
    [ Psi, Omega ] = QRPN2PsiOmega(Q,R,P,N);    
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
    
    Lmpc = M*Gamma+E;
    W = -D*C-M*Phi;
    [u_qp,fval,exitflag] = quadprog((G+G')/2,(F*Phi*Xaug(:,k)-F*Rk),Lmpc,c+W*(Xaug(:,k)),[],[],[],[],[],options_qp);
    deltaustar0k = [eye(m) zeros(m,m*(N-1))]*u_qp;
    if k > 1
        u(:,k) = u(:,k-1) + deltaustar0k;
    else
        u(:,k) = u0 + deltaustar0k;
    end
    
    Xaug(:,k+1) = A*Xaug(:,k) + B*deltaustar0k;

    Y(:,k+1) = C*Xaug(:,k+1);
    X2(:,k+1) = Am*X2(:,k) + Bm*u(:,k);
%     Xhat(:,k+1) = Am*Xhat(:,k) + B*deltaustar0k+Lobs*(Y(:,k)) + Lobs*(Y(:,k)-C*Xhat(:,k));
end

%%
t = 0:Ts:(size(Y,2)-1)*Ts;
figure(1); clf;
subplot(2,2,1)
hold on;
plot(t,Y(1,:),'LineWidth',1.5);plot(t,Y(2,:),'LineWidth',1.5);plot(t,Y(3,:),'LineWidth',1.5);
legend('Tank 1','Tank 2','Tank 3');
axis([-0.1 t(end) -0.01 0.7])
xlabel('Time [s]');
ylabel('Water level [m]');
subplot(2,2,2);
hold on;
plot(t,Y(4,:),'LineWidth',1.5);plot(t,Y(5,:),'LineWidth',1.5);
legend('Valve 1','Valve 2')
axis([-0.1 t(end) 0 Dvalve*1.05])
xlabel('Time [s]');
ylabel('Connecting pipe diameter [m]');
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
figure
plot(X2')
 