%% set variables and run model

clear all; clc; close all

load RL_busses.mat
load RL_TTS3_Controller_Blank.mat
Dvalve = 2*sqrt(5e-5/pi);

out = sim('RL_TTS3_TOP_SIM.slx');

%% postprocessing

error1 = out.References(:,1)*1000-out.Heights(:,1);
error2 = out.References(:,2)*1000-out.Heights(:,2);
error3 = out.References(:,3)*1000-out.Heights(:,3);

% norms:
L = 125/0.01;
norm2_1 = 1/L*norm(error1(1:L),2) 
norm2_2 = 1/L*norm(error2(1:L),2) 
norm2_3 = 1/L*norm(error3(1:L),2)
norminf_1 = norm(error1(1:L),Inf) 
norminf_2 = norm(error2(1:L),Inf) 
norminf_3 = norm(error3(1:L),Inf) 

[n,m] = size(out.tout);
pump1 = reshape(out.pump1,[n,1]);
pump2 = reshape(out.pump2,[n,1]);

%% plotting results

figure(1)
subplot(2,1,1)
plot(out.tout,out.References(:,1)*1000)
hold on
plot(out.tout,out.References(:,2)*1000)
plot(out.tout,out.References(:,3)*1000)
plot(out.tout,out.Heights(:,1))
plot(out.tout,out.Heights(:,2))
plot(out.tout,out.Heights(:,3))
legend('T1 ref','T2 ref','T3 ref',...
       'T1 sim','T2 sim','T3 sim','NumColumns',2 )
grid
xlabel('Time [s]')
ylabel('Height [mm]')

subplot(2,1,2)
plot(out.tout,error1)
hold on
plot(out.tout,error2)
plot(out.tout,error3)
legend('error1','error2','error3')
grid
xlabel('Time [s]')
ylabel('Error [mm]')

figure(2)
subplot(2,1,1)
plot(out.tout,pump1)
hold on
plot(out.tout,pump2)
legend('pump1','pump2')
grid
xlabel('Time [s]')
ylabel('Pump input [L/s]')

subplot(2,1,2)
plot(out.tout,out.LMpercentage)
hold on
plot(out.tout,out.MRpercentage)
legend('LM','MR')
grid
xlabel('Time [s]')
ylabel('Valve input [%]')

%% report plots

figure(3)
subplot(2,2,1)
plot(out.tout,out.Heights(:,1),'Linewidth',1.3)
hold on
plot(out.tout,out.Heights(:,2),'Linewidth',1.3)
plot(out.tout,out.Heights(:,3),'Linewidth',1.3)
legend('Tank L','Tank R','Tank M')
grid
xlabel('Time [s]')
ylabel('Tank height [mm]')

subplot(2,2,2)
plot(out.tout,out.LMpercentage,'Linewidth',1.3)
hold on
plot(out.tout,out.MRpercentage,'Linewidth',1.3)
legend('Valve LM','Valve MR')
grid
xlabel('Time [s]')
ylabel('Valve opening [%]')

subplot(2,2,3)
plot(out.tout,pump1,'Linewidth',1.3)
hold on
plot(out.tout,pump2,'Linewidth',1.3)
legend('Pump L','Pump R')
grid
xlabel('Time [s]')
ylabel('Pump input [L/s]')

subplot(2,2,4)
plot(out.tout,out.LMpercentage_ip,'Linewidth',2)
hold on
plot(out.tout,out.MRpercentage_ip,'Linewidth',2)
legend('Valve LM','Valve MR')
grid
xlabel('Time [s]')
ylabel('Valve input [%]')
