tubeSpeed = 0.0155; clc;
Ts = 0.01;
Ts2 = 4;
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S/pi);     % [m]
[Atemp,Btemp] =  PID_fcn([0.3;0.3;0.15;0.5*Dvalve;0.5*Dvalve],0.0155);
Q = [eye(3) zeros(3,2);
        zeros(2,5)];
R = 0.75*[0.5*eye(2) zeros(2,2); % 0.5
    zeros(2,2) 0.00004*eye(2)]; % 0.00005
[~,P,~] = lqrd(Atemp,Btemp,Q,R,Ts2);


save('RL_TTS3_Controller_Blank.mat','Ts');
save('nonTunables.mat','Ts2','tubeSpeed','Q','R','P');
