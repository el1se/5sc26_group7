clc; clear all; close all;

% Filling right tank, opening valve between middle/right tank


%% Set variables

% Parameters
A = 0.0154;     % [m^2]
S = 5*10^-2;    % [m^2]
Hmax = 60;      % [cm]
Qmax = 0.1;     % [L/s]

% Sample time
Ts = 0.01;  

% Time intervals
T0 = 0:Ts:5;                        % initialize
T1 = T0(end)+Ts:Ts:T0(end)+100;     % fill tank 2
T2 = T1(end)+Ts:Ts:T1(end)+120;     % open valve 3-2 (MR)

% Time
t = [T0 T1 T2];

% Pumps
SetPumps = [zeros(1,length(T0))...
            ones(1,length(T1))*Qmax...
            zeros(1,length(T2))];
            
SetPumpR = timeseries(SetPumps,t);
SetPumpL = timeseries(zeros(1,length(t)),t);

% Valve 3-2
Valve_open = [zeros(1,length([T0 T1]))...
              ones(1,length(T2))];

Valve_close = [ones(1,length([T0 T1]))...
               zeros(1,length(T2))];
           
Valve_other = zeros(1,length(t));
               
ValveLD_close = timeseries(Valve_other,t);
ValveLM_close = timeseries(Valve_other,t);
ValveMD_close = timeseries(Valve_other,t);
ValveMR_close = timeseries(Valve_close,t);
ValveRD_close = timeseries(Valve_other,t);

ValveLD_open = timeseries(Valve_other,t);
ValveLM_open = timeseries(Valve_other,t);
ValveMD_open = timeseries(Valve_other,t);
ValveMR_open = timeseries(Valve_open,t);
ValveRD_open = timeseries(Valve_other,t);

%% Check data

% Pumps
figure(1)
subplot(1,2,1);
plot(t,SetPumps)
xlabel('t'); ylabel('pump 2 input [L/s]')
subplot(1,2,2)
plot(t,cumsum(SetPumps)*Ts/A)
xlabel('t'); ylabel('level 2 [mm]')
level = trapz(t,SetPumps/A)

% Valves
figure(2)
plot(t,Valve_open,t,Valve_close)
legend('open','close')
xlabel('t'); ylabel('valve booleans')
title('valveMR')

%% Save data to .mat file

save('RL_TTS3_Controller_Blank.mat',...
     'SetPumpL','SetPumpR','Ts',...
     'ValveLD_close','ValveLM_close','ValveMD_close',...
     'ValveMR_close','ValveRD_close',...
     'ValveLD_open','ValveLM_open','ValveMD_open',...
     'ValveMR_open','ValveRD_open');
 
 