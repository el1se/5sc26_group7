%% constant
addpath('C:\Users\maxva\Google Drive\Documenten\TuE\Master\5SC26 - Integration project\5sc26_group7\GP_valves');
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S/pi);     % [m]
Ts = 0.01;
betaBar = [16.1461609250188,15.5165993538781];
[Act,Bct] = PID_fcn([0.3;0.3;0;0.5*Dvalve;0.5*Dvalve],betaBar);
Q = [eye(3) zeros(3,2);
        zeros(2,5)];
R = [0.00005*eye(2) zeros(2); % 0.005
     zeros(2) 0.00000001*eye(2)]; %0.000001
[K,~,~] = lqrd(Act,Bct,Q,R,Ts); 
save('RL_TTS3_Controller_Blank.mat','Ts');
save('nonTunables.mat','K','Dvalve','S','Dvalve');
load('RL_busses.mat')