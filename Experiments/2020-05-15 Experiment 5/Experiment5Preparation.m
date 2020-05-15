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

SetValves_vector  = boolean([zeros(1,size(input1,2))]);
SetPump_vector = timeseries(input1,t1);
SetValves_vector = timeseries(SetValves_vector,t1);
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
SetValves_vector  = boolean([zeros(1,size(input1,2)) ones(1,150/Ts) zeros(1,size(input1,2)-1)]);
SetPump_vector = timeseries(input2,t2);
SetValves_vector = timeseries(SetValves_vector,t2);

delete('ToRemoteLabs2.zip')
delete('RL_TTS3_Controller_Blank.mat')
save('RL_TTS3_Controller_Blank.mat','SetPump_vector','Ts','SetValves_vector')
zip('ToRemoteLabs2',{'RL_TTS3_Controller_Blank.mat','RL_TTS3_Controller_Blank.slx'});
kB=dir('ToRemoteLabs2.zip').bytes/1000

% if one wants toe extend it with more more more data
input3 = [input1 zeros(1,150/Ts) fliplr(input1(2:end)) zeros(1,150/Ts) (input1(2:end)) zeros(1,150/Ts) fliplr(input1(2:end))];
t3 = 0:Ts:(size(input3,2)-1)*Ts;
fprintf('duration is %f seconds', max(t3));
SetValves_vector  = boolean([zeros(1,size(input1,2)) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1)]);
SetPump_vector = timeseries(input3,t3);
SetValves_vector = timeseries(SetValves_vector,t3);

delete('ToRemoteLabs3.zip')
delete('RL_TTS3_Controller_Blank.mat')
save('RL_TTS3_Controller_Blank.mat','SetPump_vector','Ts','SetValves_vector')
zip('ToRemoteLabs3',{'RL_TTS3_Controller_Blank.mat','RL_TTS3_Controller_Blank.slx'});
kB=dir('ToRemoteLabs3.zip').bytes/1000

% if one wants toe extend it with more more more more  more data
input4 = [input1 zeros(1,150/Ts) fliplr(input1(2:end)) zeros(1,150/Ts) (input1(2:end)) zeros(1,150/Ts) fliplr(input1(2:end)) zeros(1,150/Ts) (input1(2:end)) zeros(1,150/Ts) fliplr(input1(2:end))];
t4 = 0:Ts:(size(input4,2)-1)*Ts;
fprintf('duration is %f seconds, %f minutes', max(t4), max(t4)/60);
SetValves_vector  = boolean([zeros(1,size(input1,2)) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1)]);
SetPump_vector = timeseries(input4,t4);
SetValves_vector = timeseries(SetValves_vector,t4);

delete('ToRemoteLabs4.zip')
delete('RL_TTS3_Controller_Blank.mat')
save('RL_TTS3_Controller_Blank.mat','SetPump_vector','Ts','SetValves_vector')
zip('ToRemoteLabs4',{'RL_TTS3_Controller_Blank.mat','RL_TTS3_Controller_Blank.slx'});
kB=dir('ToRemoteLabs4.zip').bytes/1000

return


