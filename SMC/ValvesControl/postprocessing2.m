%% set variables and run model

clear all; clc; close all

load RL_busses.mat
load RL_TTS3_Controller_Blank.mat

Dvalve = 2*sqrt(5e-5/pi);

%Tend = 125
Tend = 110

out = sim('RL_TTS3_TOP_SIM.slx');

%% postprocessing

error1 = out.References(:,1)*1000-out.Heights(:,1);
error2 = out.References(:,2)*1000-out.Heights(:,2);
error3 = out.References(:,3)*1000-out.Heights(:,3);

% norms:
L = Tend/0.01;
norm2_1 = 1/L*norm(error1(1:L),2) 
norm2_2 = 1/L*norm(error2(1:L),2) 
norm2_3 = 1/L*norm(error3(1:L),2)
norminf_1 = error1(L) %norm(error1(1:L),Inf) 
norminf_2 = error2(L) %norm(error2(1:L),Inf) 
norminf_3 = error3(L) %norm(error3(1:L),Inf) 

[n,m] = size(out.tout);
pump1 = reshape(out.pump1,[n,1]);
pump2 = reshape(out.pump2,[n,1]);

%% report plots

figure(1)
subplot(2,2,1)
plot(out.tout,out.Heights(:,1),'Linewidth',1.3)
hold on
plot(out.tout,out.Heights(:,2),'Linewidth',1.3)
plot(out.tout,out.Heights(:,3),'Linewidth',1.3)
legend('Tank L','Tank R','Tank M')
grid
xlim([0 Tend])
xlabel('Time [s]')
ylabel('Tank height [mm]')

subplot(2,2,2)
plot(out.tout,out.LMpercentage,'Linewidth',1.3)
hold on
plot(out.tout,out.MRpercentage,'Linewidth',1.3)
legend('Valve LM','Valve MR')
grid
xlim([0 Tend])
xlabel('Time [s]')
ylabel('Valve opening [%]')

subplot(2,2,3)
plot(out.tout,error1,'Linewidth',1.3)
hold on
plot(out.tout,error2,'Linewidth',1.3)
plot(out.tout,error3,'Linewidth',1.3)
legend('Tank L','Tank R','Tank M')
grid
xlim([0 Tend])
xlabel('Time [s]')
ylabel('Error [mm]')

subplot(2,2,4)
plot(out.tout,pump1,'Linewidth',1.3)
hold on
plot(out.tout,pump2,'Linewidth',1.3)
legend('Pump L','Pump R')
grid
xlim([0 Tend])
xlabel('Time [s]')
ylabel('Pump input [L/s]')
