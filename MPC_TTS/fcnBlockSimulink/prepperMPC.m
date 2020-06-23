tubeSpeed = 0.0155; clc;
Ts = 0.5;
Ts2 = Ts;
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S/pi);     % [m]
bB = [16.1461609250188,15.5165993538781]./1000;
[Atemp,Btemp] =  PID_fcn([0.3;0.3;0.15;0.5*Dvalve;0.5*Dvalve],bB);
% Q = [1*eye(3) zeros(3,2);
%         zeros(2,3) 0*eye(2)];
Q = [1*eye(2) zeros(2,3);
     0 0 2e0 0 0;
zeros(2,3) 0*eye(2)];
R = [1e1*eye(2) zeros(2,2); %   For delta input 1e1
    zeros(2,2) 0.8e-4*eye(2)]; % 0.8e-4
Ru = [5e4*eye(2) zeros(2,2); % for input (P) % 5e4
    zeros(2,2) 1e0*eye(2)];  %1e0


[~,P,~] = lqrd(Atemp,Btemp,Q,Ru,Ts2);

save('RL_TTS3_Controller_Blank.mat','Ts');
save('nonTunables.mat','Ts2','tubeSpeed','Q','R','P','bB');
load('RL_busses.mat');