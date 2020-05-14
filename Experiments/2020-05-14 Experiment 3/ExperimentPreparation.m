clear all
%%
% control the valves
% ValvesVecClosed = [LeakX1 x1x3Valve LeakX3 x3x2Valve LeakX2] % 1 means closed
ValvesVecClosed = [1 1 1 1 1];
ValveLD_close = ValvesVecClosed(1);
ValveLM_close = ValvesVecClosed(2);
ValveMD_close = ValvesVecClosed(3);
ValveMR_close = ValvesVecClosed(4);
ValveRD_close = ValvesVecClosed(5);
ValveLD_open = double(~ValvesVecClosed(1));
ValveLM_open = double(~ValvesVecClosed(2));
ValveMD_open = double(~ValvesVecClosed(3));
ValveMR_open = double(~ValvesVecClosed(4));
ValveRD_open = double(~ValvesVecClosed(5));

save('RL_TTS3_Controller_Blank.mat','ValveLD_close','ValveLM_close',...
    'ValveMD_close','ValveMR_close','ValveRD_close','ValveLD_open',...
    'ValveLM_open','ValveMD_open','ValveMR_open','ValveRD_open','-append');

clear
load('RL_TTS3_Controller_Blank.mat')

%% pump input
clc
% assume experiment is 5 min
Ts= 0.01; % sampling time
f = 1/4; % frequency
mag = 0.01:0.01:0.1; % nr of different amplitudes
periods = 8; % nr of periods each magnitude occurs
t = 0:Ts:size(mag,2)*periods/f; % duration time

% creat blog signal with frequency
data1 = max(0,sign(sin(2*pi*f*t-0.5*pi))); 

% add a ramp to the input signal
data1 = data1 .*[0 repelem(mag,periods/(f*Ts))];
  
% plotting
subplot(1,2,1);
plot(t,data1); title('input')
subplot(1,2,2); plot(t,cumsum(data1)*0.01/0.0154); title('level'); ylabel('mm');yline(600); 

% total tank height
trapz(t,data1/0.0154)

%% making inputs for experiment
t2 = 0:Ts:max(t)*3;
data2 = [data1 data1(2:end)*0 fliplr(data1(2:end))];
subplot(1,2,1);
plot(t2,data2); title('input')
subplot(1,2,2); plot(t2,cumsum(data2)*0.01/0.0154); title('level'); ylabel('mm');yline(600); 

% vectors for the valves
ValveVec = [ones(1,size(t,2)) zeros(1,size(t,2)-1) ones(1,size(t,2)-1)];
ValveVec = [ones(1,size(t,2))]>0;

ValveLD_close_Vector = timeseries(ValveVec,t);
ValveLM_close_Vector = timeseries(ValveVec,t);
ValveMD_close_Vector = timeseries(ValveVec,t);
ValveMR_close_Vector = timeseries(ValveVec,t);
ValveRD_close_Vector = timeseries(ValveVec,t);
ValveLD_open_vector = timeseries(double(~ValveVec),t);
ValveLM_open_vector = timeseries(double(~ValveVec),t);
ValveMD_open_vector = timeseries(double(~ValveVec),t);
ValveMR_open_vector = timeseries(double(~ValveVec),t);
ValveRD_open_vector = timeseries(double(~ValveVec),t);

SetPumpR_vector =  timeseries((data1),t);
SetPumpL_vector =  timeseries((data1),t);
save('RL_TTS3_Controller_Blank.mat','SetPumpL_vector','SetPumpR_vector','-append');
save('RL_TTS3_Controller_Blank.mat','ValveLD_close_Vector','ValveLM_close_Vector',...
    'ValveMD_close_Vector','ValveMR_close_Vector','ValveRD_close_Vector',...
    'ValveLD_open_vector','ValveLM_open_vector','ValveMD_open_vector',...
    'ValveMR_open_vector','ValveRD_open_vector','-append');


% clear
% load('RL_TTS3_Controller_Blank.mat')

