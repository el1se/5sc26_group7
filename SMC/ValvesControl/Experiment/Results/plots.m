clear all; clc; close all

% results path
base = "C:\Users\edtve\Documents\TUe Systems and Control\Integration Project\Code\GitHub\5sc26_group7\SMC\ValvesControl\Experiment\Results"

% measurement data
%load(fullfile(base,"same_refs1/5217_results.mat"))
%load(fullfile(base,"same_refs2/5940_results.mat"))
%load(fullfile(base,"same_refs3/5942_results.mat"))
%load(fullfile(base,"same_refs4/5948_results.mat"))
%load(fullfile(base,"different_refs1/5218_results.mat"))
load(fullfile(base,"different_refs2/5946_results.mat"))
%load(fullfile(base,"different_refs3/5950_results.mat"))
%load(fullfile(base,"different_refs4/5952_results.mat"))
%load(fullfile(base,"different_refs5/5953_results.mat"))

% reference data
%load(fullfile(base,"same_refs1/sim_out_references.mat"))
load(fullfile(base,"different_refs1/sim_out_references.mat"))

% duration
%Tend = 125;
Tend = 125;
L = Tend/0.01;

% get data from measurement file
time  = mdfData.RunTime;
pump1 = mdfData.ModelRoot_controller_controller_ModelRoot_PumpL_In1;
pump2 = mdfData.ModelRoot_controller_controller_ModelRoot_PumpR_In1;

H1 = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightL_mm__In1;
H3 = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightM_mm__In1;
H2 = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightR_mm__In1;

LMpercentage_ip = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeValveLMIp____In1;
MRpercentage_ip = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeValveMRIp___1_In;
LMpercentage = mdfData.ModelRoot_controller_controller_ModelRoot_Valves_ScopeValveLM__;
MRpercentage = mdfData.ModelRoot_controller_controller_ModelRoot_Valves_ScopeValveMR__;

%% norms

error1 = out.References(1:L,1)*1000-H1(1:L);
error2 = out.References(1:L,2)*1000-H2(1:L);
error3 = out.References(1:L,3)*1000-H3(1:L);

norm2_1 = 1/L*norm(error1(1:L),2) 
norm2_2 = 1/L*norm(error2(1:L),2) 
norm2_3 = 1/L*norm(error3(1:L),2)
norminf_1 = error1(L) %norm(error1(1:L),Inf) 
norminf_2 = error2(L) %norm(error2(1:L),Inf) 
norminf_3 = error3(L) %norm(error3(1:L),Inf) 

%% report plots

figure(1)
subplot(2,2,1)
plot(time,H1,'Linewidth',1.3)
hold on
plot(time,H2,'Linewidth',1.3)
plot(time,H3,'Linewidth',1.3)
legend('Tank L','Tank R','Tank M')
grid
xlim([0 Tend])
xlabel('Time [s]')
ylabel('Tank height [mm]')

subplot(2,2,2)
plot(time,LMpercentage,'Linewidth',1.3)
hold on
plot(time,MRpercentage,'Linewidth',1.3)
legend('Valve LM','Valve MR')
grid
xlim([0 Tend])
xlabel('Time [s]')
ylabel('Valve opening [%]')

subplot(2,2,3)
plot(time(1:L),error1,'Linewidth',1.3)
hold on
plot(time(1:L),error2,'Linewidth',1.3)
plot(time(1:L),error3,'Linewidth',1.3)
legend('Tank L','Tank R','Tank M')
grid
xlim([0 Tend])
xlabel('Time [s]')
ylabel('Error [mm]')

subplot(2,2,4)
plot(time,pump1,'Linewidth',1.3)
hold on
plot(time,pump2,'Linewidth',1.3)
legend('Pump L','Pump R')
grid
xlim([0 Tend])
xlabel('Time [s]')
ylabel('Pump input [L/s]')

% plot(time,LMpercentage_ip,'Linewidth',2)
% hold on
% plot(time,MRpercentage_ip,'Linewidth',2)
% legend('Valve LM','Valve MR')
% grid
% xlim([0 Tend])
% xlabel('Time [s]')
% ylabel('Valve input [%]')
