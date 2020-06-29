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
plot(t(1:I),X(1:I,1),'Linewidth',LW); hold on;
plot(t(1:I),X(1:I,2),'Linewidth',LW);
plot(t(1:I),X(1:I,3),'Linewidth',LW);
grid on;
xlim([0 125])
xlabel('Time [s]');
ylabel('Tank heights [mm]')
subplot(2,2,3)
plot(t(1:I),error(1:I,1),'Linewidth',LW); hold on;
plot(t(1:I),error(1:I,2),'Linewidth',LW);
plot(t(1:I),error(1:I,3),'Linewidth',LW);
xlim([0 125])
grid on;
legend('Tank L','Tank R','Tank M');
xlabel('Time [s]');
ylabel('Error[mm]');
subplot(2,2,2)
plot(t(1:I),X(1:I,4),'Linewidth',LW); hold on;
plot(t(1:I),X(1:I,5),'Linewidth',LW);
xlim([0 125])
grid on;
legend('Valve LM','Valve RM');
xlabel('Time [s]');
ylabel('Valve opening [%]')
subplot(2,2,4)
plot(t(1:I),U(1:I,1),'Linewidth',LW); hold on;
plot(t(1:I),U(1:I,2),'Linewidth',LW);grid on;
xlim([0 125])
legend('Pump L','Pump R')
xlabel('Time [s]');
ylabel('Pump input [L/s]');



%%
1/I*norm(error(1:I,1),2)
1/I*norm(error(1:I,2),2)
1/I*norm(error(1:I,3),2)
error(end,:)