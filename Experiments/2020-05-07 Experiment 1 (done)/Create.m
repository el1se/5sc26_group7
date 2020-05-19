clear all

save('x1data.mat')
save('x2data.mat')
save('x3data.mat')

load('RL_TTS3_Controller_Blank.mat')



% if simtime = 200 with ts =0.01 you need 200001 array
%%
StopTime = 30; % in secs

load('')

%%
t = 0:Ts:200;
data = 0.1*max(sign(sin(2*pi*0.1*t)),0);
plot(data)
%%

SetPumpR_vector = [data;0:Ts:200];

