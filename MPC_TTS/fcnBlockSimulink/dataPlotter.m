% clear all; 
close all; clc;
clear X; clear U;
load('MPC040302FINALRESULTS.mat');
t = mdfData.RunTime;
X(:,1) = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightL_mm__In1;
X(:,2) = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightR_mm__In1;
X(:,3) = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightM_mm__In1;
X(:,4) = mdfData.ModelRoot_controller_controller_ModelRoot_Valves_ScopeValveLM__;
X(:,5) = mdfData.ModelRoot_controller_controller_ModelRoot_Valves_ScopeValveMR__;
U(:,1) = mdfData.ModelRoot_controller_controller_ModelRoot_PumpL_In1;
U(:,2) = mdfData.ModelRoot_controller_controller_ModelRoot_PumpR_In1;
% U(:,3) = mdfData.ModelRoot_controller_controller_ModelRoot_MLValveInput_In1;
% U(:,4) = mdfData.ModelRoot_controller_controller_ModelRoot_MRValveInput1_In1;
I=221;
LW = 1.3;

ref = out.simout;
error = ref*1000-X(1:I,1:3);
1/I*norm(error(:,1),2)
1/I*norm(error(:,2),2)
1/I*norm(error(:,3),2)
error(end,:)

figure
subplot(2,2,1)
plot(t(1:I),X(1:I,1),'Linewidth',LW); hold on;
plot(t(1:I),X(1:I,2),'Linewidth',LW);
plot(t(1:I),X(1:I,3),'Linewidth',LW);
grid on;
xlim([0 t(I)])
xlabel('Time [s]');
ylabel('Tank heights [mm]')
subplot(2,2,3)
plot(t(1:I),error(1:I,1),'Linewidth',LW); hold on;
plot(t(1:I),error(1:I,2),'Linewidth',LW);
plot(t(1:I),error(1:I,3),'Linewidth',LW);
xlim([0 t(I)])
grid on;
legend('Tank L','Tank R','Tank M');
xlabel('Time [s]');
ylabel('Error[mm]');
subplot(2,2,2)
plot(t(1:I),X(1:I,4),'Linewidth',LW); hold on;
plot(t(1:I),X(1:I,5),'Linewidth',LW);
xlim([0 t(I)])
grid on;
legend('Valve LM','Valve RM');
xlabel('Time [s]');
ylabel('Valve opening [%]')
subplot(2,2,4)
plot(t(1:I),U(1:I,1),'Linewidth',LW); hold on;
plot(t(1:I),U(1:I,2),'Linewidth',LW);grid on;
xlim([0 t(I)])
legend('Pump L','Pump R')
xlabel('Time [s]');
ylabel('Pump input [L/s]');


%%

% subplot(2,2,4)
% plot(t,U(:,3),'Linewidth',LW); hold on;
% plot(t,U(:,4),'Linewidth',LW);
% legend('Valve ML','Valve MR')
% xlabel('Time [s]');
% ylabel('Valve input [%]');
% 
% % 
% figure
% subplot(1,3,1)
% plot(t,X(:,1),'Linewidth',LW); hold on;
% plot(out.tout,out.simout.Data(:,1),'Linewidth',LW);
% xlabel('Time [s]');
% ylabel('Tank heights [mm]')
% subplot(1,3,2)
% plot(t,X(:,2),'Linewidth',LW); hold on;
% plot(out.tout,out.simout.Data(:,2),'Linewidth',LW);
% xlabel('Time [s]');
% ylabel('Tank heights [mm]')
% subplot(1,3,3)
% plot(t,X(:,3),'Linewidth',LW); hold on;
% plot(out.tout,out.simout.Data(:,3),'Linewidth',LW);
% xlabel('Time [s]');
% ylabel('Tank heights [mm]')
% legend('Real system','Simulation (tubeSpeed=0.0155)')
% 
% 
% figure
% subplot(1,3,1)
% plot(t,X(:,1)-out.simout.Data(:,1),'Linewidth',LW); 
% xlabel('Time [s]');
% ylabel('Residual between model and system [mm]')
% subplot(1,3,2)
% plot(t,X(:,2)-out.simout.Data(:,2),'Linewidth',LW); 
% xlabel('Time [s]');
% ylabel('Residual between model and system [mm]')
% subplot(1,3,3)
% plot(t,X(:,3)-out.simout.Data(:,3),'Linewidth',LW); 
% xlabel('Time [s]');
% ylabel('Residual between model and system [mm]')
% %%
% figure
% plot(t,X(:,1)-out.simout.Data(:,4)*1000); hold on
% norm(X(:,1)-out.simout.Data(:,4)*1000,2)
% norm(X(:,2)-out.simout.Data(:,4)*1000,2)
% norm(X(:,3)-out.simout.Data(:,4)*1000,2)