clear all;
tubeSpeed = 0.0155; clc; 
Ts = 0.5;
Ts2 = Ts;
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S/pi);     % [m]
bB = [16.1461609250188,15.5165993538781]./1000;
[Atemp,Btemp] =  PID_fcn([0.3;0.3;0.15;0.5*Dvalve;0.5*Dvalve],bB);
% Q = [1*eye(3) zeros(3,2);
%         zeros(2,3) 0*eye(2)];
Q = [1*eye(3) zeros(3,2);
zeros(2,3) 0*eye(2)];
% R = [1e1*eye(2) zeros(2,2); %   For delta input 1e1
%     zeros(2,2) 0.8e-4*eye(2)]; % 0.8e-4
R = [1e1*eye(2) zeros(2);
    zeros(2) 1e-5*eye(2)];
% Ru = [5e4*eye(2) zeros(2,2); % for input (P) % 5e4
%     zeros(2,2) 1e0*eye(2)];  %1e0
Ru = [5e4*eye(2) zeros(2);
    zeros(2) 1e-5*eye(2)];


% [~,P,~] = lqrd(Atemp,Btemp,Q,Ru,Ts2);
% P=1e3*eye(5   );
P = [];
steps = 4;
for i = 1:steps
    for j = 1:steps
        [At,Bt] =  PID_fcn([0.3;0.3;0.15;i*1/steps*Dvalve;j*1/steps*Dvalve],bB);
        [~,Pt(:,(j-1)*5+1:j*5),~] = lqrd(At,Bt,Q,Ru,Ts2);
        indexer(i,j) = (i-1)*steps+j;
    end
    P = [P Pt];
end

save('RL_TTS3_Controller_Blank.mat','Ts');
save('nonTunables.mat','Ts2','tubeSpeed','Q','R','P','bB','Dvalve','steps','indexer');
load('RL_busses.mat');