% clear all; 
close all; clc;
clear X; clear U;
load('LQRresultsDiffTraj.mat');
t = mdfData.RunTime;
X(:,1) = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightL_mm__In1;
X(:,2) = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightR_mm__In1;
X(:,3) = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightM_mm__In1;
X(:,4) = mdfData.ModelRoot_controller_controller_ModelRoot_Valves_ScopeValveLM__;
X(:,5) = mdfData.ModelRoot_controller_controller_ModelRoot_Valves_ScopeValveMR__;
U(:,1) = mdfData.ModelRoot_controller_controller_ModelRoot_PumpL_In1;
U(:,2) = mdfData.ModelRoot_controller_controller_ModelRoot_PumpR_In1;
U(:,3) = mdfData.ModelRoot_controller_controller_ModelRoot_MLValveInput_In1;
U(:,4) = mdfData.ModelRoot_controller_controller_ModelRoot_MRValveInput1_In1;
dat = out.simout.Data;
L = length(t);
error = dat(1:L,10:12)*1000-X(:,1:3);
LW = 1.3; 
I=11001;
figure
subplot(2,2,1)
plot(t,X(:,1),'Linewidth',LW); hold on;
plot(t,X(:,2),'Linewidth',LW);
plot(t,X(:,3),'Linewidth',LW);
% plot(t,zeros(length(t),1)+305,'r--');
% plot(t,zeros(length(t),1)+295,'r--');
xlim([0 t(I)])
legend('Tank L','Tank R','Tank M');
xlabel('Time [s]');
ylabel('Tank heights [mm]')
subplot(2,2,2)
plot(t,X(:,4),'Linewidth',LW); hold on;
plot(t,X(:,5),'Linewidth',LW);
xlim([0 t(I)])
legend('Valve LM','Valve RM');
xlabel('Time [s]');
ylabel('Valve opening [%]')
subplot(2,2,3)
plot(t,U(:,1),'Linewidth',LW); hold on;
plot(t,U(:,2),'Linewidth',LW);
xlim([0 t(I)])
legend('Pump L','Pump R')
xlabel('Time [s]');
ylabel('Pump input [L/s]');
subplot(2,2,4)
plot(t,U(:,3),'Linewidth',LW); hold on;
plot(t,U(:,4),'Linewidth',LW);
xlim([0 t(I)])
legend('Valve ML','Valve MR')
xlabel('Time [s]');
ylabel('Valve input [%]');



%%
1/I*norm(error(1:I,1),2)
1/I*norm(error(1:I,2),2)
1/I*norm(error(1:I,3),2)
error(end,:)