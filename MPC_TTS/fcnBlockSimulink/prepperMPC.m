tubeSpeed = 0.0155; clc;
Ts = 0.5;
Ts2 = Ts;
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S/pi);     % [m]
bB = [16.1461609250188,15.5165993538781]./1000;
[Atemp,Btemp] =  PID_fcn([0.2;0.2;0.1;0.5*Dvalve;0.5*Dvalve],bB);
Q = [1*eye(3) zeros(3,2);
        zeros(2,3) 0*eye(2)];
R = [1e-3/0.1*eye(2) zeros(2,2); % For delta input
    zeros(2,2) 1e-4/100*eye(2)]; 
Ru = [1e-3/0.1*eye(2) zeros(2,2); % for input (P)
    zeros(2,2) 1e-5/100*eye(2)]; 


[~,P,~] = lqrd(Atemp,Btemp,Q,Ru,Ts2);

save('RL_TTS3_Controller_Blank.mat','Ts');
save('nonTunables.mat','Ts2','tubeSpeed','Q','R','P','bB');
load('RL_busses.mat');