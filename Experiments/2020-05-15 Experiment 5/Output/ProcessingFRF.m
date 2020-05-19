clearvars
clc
%%
load('ToRemoteLabs5FRF_results.mat')

referenceBoolean = boolean([zeros(1,15400) ones(1,150/0.01) zeros(1,15400) ones(1,150/0.01) zeros(1,15400) ones(1,150/0.01) zeros(1,15400) ones(1,150/0.01) zeros(1,15400) ones(1,150/0.01) zeros(1,15400)]);


data.time = seconds(mdfData.Time(1:length(referenceBoolean)));
data.x1 = double(mdfData.x1_mm(1:length(referenceBoolean)));


load('C:\Users\Bart\Documents\5.4 integration project 5SC26\GitHubFolderStuff\5sc26_group7\Experiments\2020-05-15 Experiment 5\InputFRF.mat')
input5 = input5';
integral_input5 = cumsum(input5)*0.01/0.0154;
% plot(cumsum(input5)*0.01/0.0154,'linewidth',1)
input5new = [input5(1:15400) zeros(1,15000) input5(15400+1:30800) zeros(1,15000)...
    input5(30800+1:46200) zeros(1,15000) input5(46200+1:61600) zeros(1,15000)...
     input5(61600+1:77000) zeros(1,15000)  input5(77000+1:92400)];
input5= input5new;

data.time = data.time(~referenceBoolean);
data.x1 = data.x1(~referenceBoolean);
data.u = input5(~referenceBoolean)';

nk=1;
nbrange = 1; % is right one
nfrange = 1; % is right one
ncrange = 1;
ndrange = 1;

c=0;
dataxy = iddata(data.x1,data.u,0.01);
G_bj = bj(dataxy, [1,1,1,1,1], bjOptions('Focus','prediction'));%,'WeightingFilter',[0.1 1.9])); % [nb nc nd nf nk])
figure()
resid(dataxy,G_bj);
plot(data.time,data.x1)
hold on
% plot(data.time,data.cumsum(input5)*0.01/0.0154,'linewidth',1)

input = input5;
tijd  = seconds(mdfData.Time(1:length(input5)))';

lsim(tf(G_bj.B,G_bj.F,0.01),input,tijd)
legend('x1','sim x1')