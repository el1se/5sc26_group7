clearvars
%% input rate constraint
% it was not optimally visible las time.
% new approach is to make a slightly different signal.

f = 1/8; % frequency
mag = 0.02:0.02:0.1; % nr of different amplitudes
periods = 6; % nr of periods each magnitude occurs
factor = 1.5; % each period appears 1.5 times
Ts= 0.01; % sampling frequency.
UnitBlock = [zeros(1,1/(f*Ts)) ones(1,1/(f*2*Ts))];

% making input signal
n = periods*size(mag,2);
input1 = zeros(1,n*size(UnitBlock,2)+1);
for i = 1:size(mag,2)
    MagUnit = repmat(UnitBlock,1,periods)*mag(i);
    N = size(MagUnit,2);
    input1(1,N*(i-1)+2:N*i+1) = MagUnit;
end
t1 = 0:Ts:(factor*periods*size(mag,2)/f);
fprintf('duration is %f seconds', max(t1));
x1andx3 = trapz(t1,input1/0.0154)

SetValves_vector1  = boolean([zeros(1,size(input1,2))]);
SetPump_vector = timeseries(input1,t1);
SetValves_vector = timeseries(SetValves_vector1,t1);
delete('ToRemoteLabs1.zip')
delete('RL_TTS3_Controller_Blank.mat')
save('RL_TTS3_Controller_Blank.mat','SetPump_vector','Ts','SetValves_vector')
zip('ToRemoteLabs1',{'RL_TTS3_Controller_Blank.mat','RL_TTS3_Controller_Blank.slx'});
kB=dir('ToRemoteLabs1.zip').bytes/1000

% time of emptying tankes takes about 150 seconds
% if one wants toe extend it with more data
input2 = [input1 zeros(1,150/Ts) fliplr(input1(2:end))];
t2 = 0:Ts:(size(input2,2)-1)*Ts;
fprintf('duration is %f seconds', max(t2));
SetValves_vector2  = boolean([zeros(1,size(input1,2)) ones(1,150/Ts) zeros(1,size(input1,2)-1)]);
SetPump_vector = timeseries(input2,t2);
SetValves_vector = timeseries(SetValves_vector2,t2);

delete('ToRemoteLabs2.zip')
delete('RL_TTS3_Controller_Blank.mat')
save('RL_TTS3_Controller_Blank.mat','SetPump_vector','Ts','SetValves_vector')
zip('ToRemoteLabs2',{'RL_TTS3_Controller_Blank.mat','RL_TTS3_Controller_Blank.slx'});
kB=dir('ToRemoteLabs2.zip').bytes/1000

% if one wants toe extend it with more more more data
input3 = [input1 zeros(1,150/Ts) fliplr(input1(2:end)) zeros(1,150/Ts) (input1(2:end)) zeros(1,150/Ts) fliplr(input1(2:end))];
t3 = 0:Ts:(size(input3,2)-1)*Ts;
fprintf('duration is %f seconds', max(t3));
SetValves_vector3  = boolean([zeros(1,size(input1,2)) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1)]);
SetPump_vector = timeseries(input3,t3);
SetValves_vector = timeseries(SetValves_vector3,t3);

delete('ToRemoteLabs3.zip')
delete('RL_TTS3_Controller_Blank.mat')
save('RL_TTS3_Controller_Blank.mat','SetPump_vector','Ts','SetValves_vector')
zip('ToRemoteLabs3',{'RL_TTS3_Controller_Blank.mat','RL_TTS3_Controller_Blank.slx'});
kB=dir('ToRemoteLabs3.zip').bytes/1000

% if one wants toe extend it with more more more more  more data
input4 = [input1 zeros(1,150/Ts) fliplr(input1(2:end)) zeros(1,150/Ts) (input1(2:end)) zeros(1,150/Ts) fliplr(input1(2:end)) zeros(1,150/Ts) (input1(2:end)) zeros(1,150/Ts) fliplr(input1(2:end))];
t4 = 0:Ts:(size(input4,2)-1)*Ts;
fprintf('duration is %f seconds, %f minutes', max(t4), max(t4)/60);
SetValves_vector4  = boolean([zeros(1,size(input1,2)) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1)]);
SetPump_vector = timeseries(input4,t4);
SetValves_vector = timeseries(SetValves_vector4,t4);

delete('ToRemoteLabs4.zip')
delete('RL_TTS3_Controller_Blank.mat')
save('RL_TTS3_Controller_Blank.mat','SetPump_vector','Ts','SetValves_vector')
zip('ToRemoteLabs4',{'RL_TTS3_Controller_Blank.mat','RL_TTS3_Controller_Blank.slx'});
kB=dir('ToRemoteLabs4.zip').bytes/1000

% save('Output/input.mat','input1','input2','input3','input4','SetValves_vector1','SetValves_vector2','SetValves_vector3','SetValves_vector4','t1','t2','t3','t4','t5','SetValves_vector5','input5')
return
%% FRF EXPERIMENT PREPERATION
% maximum input rate was observed to be 0.3 L/s2 based on 1 point
% Range = [0 0.1];
% N = 180000;
% input5 = idinput([N,1,1],'rbs',[],Range);
% load('InputFRF.mat')
% input5 = input5';
% integral_input5 = cumsum(input5)*0.0cle1/0.0154;
% plot(cumsum(input5)*0.01/0.0154,'linewidth',1)
% input5new = [input5(1:15400) zeros(1,15000) input5(15400+1:30800) zeros(1,15000)...
%     input5(30800+1:46200) zeros(1,15000) input5(46200+1:61600) zeros(1,15000)...
%      input5(61600+1:77000) zeros(1,15000)  input5(77000+1:92400)];
% input5= input5new;
% SetValves_vector5  = boolean([zeros(1,15400) ones(1,150/Ts) zeros(1,15400) ones(1,150/Ts) zeros(1,15400) ones(1,150/Ts) zeros(1,15400) ones(1,150/Ts) zeros(1,15400) ones(1,150/Ts) zeros(1,15400)]);
% t5 = 0:0.01:size(input5,2)*0.01-0.01;
% Ts= 0.01; % sampling frequency.
% SetValves_vector = timeseries(SetValves_vector,t5);
% SetPump_vector = timeseries(input5,t5);
% 
% fprintf('duration is %f seconds, %f minutes', max(t5), max(t5)/60);
% delete('ToRemoteLabs5.zip')
% delete('RL_TTS3_Controller_Blank.mat')
% save('RL_TTS3_Controller_Blank.mat','SetPump_vector','Ts','SetValves_vector')
% zip('ToRemoteLabs5',{'RL_TTS3_Controller_Blank.mat','RL_TTS3_Controller_Blank.slx'});
% kB=dir('ToRemoteLabs5.zip').bytes/1000