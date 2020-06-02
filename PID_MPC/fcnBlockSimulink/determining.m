%% constant
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S/pi);     % [m]
Ts = 0.1;
[Act,Bct] = PID_fcn([0.3;0.3;0;0.5*Dvalve;0.5*Dvalve]);
Q = [eye(3) zeros(3,2);
        zeros(2,5)];
R = [eye(2) zeros(2);
     zeros(2) 0.001*eye(2)];
[K,~,~] = lqrd(Act,Bct,Q,R,Ts); 
save('RL_TTS3_Controller_Blank.mat','Ts','K','Dvalve');