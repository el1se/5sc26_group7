clc; clear all; close all;

%% Load data

% input
load('2197_input.mat');

% output
load('2197_results.mat');
names = mdfData.Properties.VariableNames;
names{1,1} = 'x1_mm';
names{1,2} = 'x3_mm';
names{1,3} = 'x2_mm';
names{1,4} = 'LD';
names{1,5} = 'LM';
names{1,6} = 'MD';
names{1,7} = 'MR';
names{1,8} = 'RD';
mdfData.Properties.VariableNames = names;

%% Plot data
mdfData.Time(mdfData.LD>0)

% valve positions
figure(1)
subplot(211)
plot(mdfData.Time,mdfData.LD); hold on;
plot(mdfData.Time,mdfData.LM)
plot(mdfData.Time,mdfData.MD)
plot(mdfData.Time,mdfData.MR)
plot(mdfData.Time,mdfData.RD)
legend('LD','LM','MD','MR','RD')
xlabel('Time [s]')
ylabel('Valve positions [%]')
grid
%xlim([0 450])

% input booleans
subplot(212)
plot(ValveLD_open); hold on;
plot(ValveLD_close)
legend('open','close');
xlabel('Time [s]')
ylabel('Valve booleans [-]')
grid
