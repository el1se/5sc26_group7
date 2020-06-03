%% constant
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S/pi);     % [m]
Ts = 0.01;
tubeSpeed = 0.0155;
[Act,Bct] = PID_fcn([0.3;0.3;0;0.5*Dvalve;0.5*Dvalve],tubeSpeed);
Q = [eye(3) zeros(3,2);
        zeros(2,5)];
R = [0.05*eye(2) zeros(2);
     zeros(2) 0.0001*eye(2)];
[K,~,~] = lqrd(Act,Bct,Q,R,Ts); 
save('RL_TTS3_Controller_Blank.mat','Ts');
save('nonTunables.mat','K','Dvalve');