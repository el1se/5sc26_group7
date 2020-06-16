clc; clear all; close all

load LQRresultsNormalTRaj.mat
%load LQRresultsDiffTraj.mat

load RL_busses.mat

time  = mdfData.RunTime;
dataL = mdfData.ModelRoot_controller_controller_ModelRoot_PumpL_In1;
pumpL = timeseries(dataL,time);
dataR = mdfData.ModelRoot_controller_controller_ModelRoot_PumpR_In1;
pumpR = timeseries(dataR,time);

Dvalve = 2*sqrt(5e-5/pi);

expHeightL = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightL_mm__In1;
expHeightM = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightM_mm__In1;
expHeightR = mdfData.ModelRoot_controller_controller_ModelRoot_ScopeHeightR_mm__In1;

% better valves hack
dataVL     = mdfData.ModelRoot_controller_controller_ModelRoot_MLValveInput_In1;
valveMLin  = timeseries(dataVL,time);
dataVR     = mdfData.ModelRoot_controller_controller_ModelRoot_MRValveInput1_In1;
valveMRin  = timeseries(dataVR,time);
dataVL     = mdfData.ModelRoot_controller_controller_ModelRoot_Valves_ScopeValveLM__;
valveMLpos = timeseries(dataVL,time);
dataVR     = mdfData.ModelRoot_controller_controller_ModelRoot_Valves_ScopeValveMR__;
valveMRpos = timeseries(dataVR,time);

tstep = time(2)-time(1);
tend  = time(end) + tstep; %uncomment for NormalTraj
out   = sim('RL_TTS3_System_Blank.slx');

simHeightL = out.simHeightL;
simHeightM = out.simHeightM;
simHeightR = out.simHeightR;

resHeightL = (expHeightL - simHeightL);%./std(expHeightL - simHeightL);
resHeightM = (expHeightM - simHeightM);%./std(expHeightM - simHeightM);
resHeightR = (expHeightR - simHeightR);%./std(expHeightR - simHeightR);

figure(1)
subplot(231)
plot(simHeightL,expHeightL,'.')
hold on
line([0 simHeightL(end)],[0 simHeightL(end)])
xlabel('model height L [mm]')
ylabel('experiment height L [mm]')
subplot(234)
plot(simHeightL,resHeightL,'.')
hold on
line([0 simHeightL(end)],[0 0])
ylim([-6 6])
xlabel('model height L [mm]')
ylabel('Residuals [mm]')

%figure(2)
subplot(232)
plot(simHeightM,expHeightM,'.')
hold on
line([0 simHeightM(end)],[0 simHeightM(end)])
xlabel('model height M [mm]')
ylabel('experiment height M [mm]')
subplot(235)
plot(simHeightM,resHeightM,'.')
hold on
line([0 simHeightM(end)],[0 0])
ylim([-6 6])
xlabel('model height M [mm]')
ylabel('Residuals [mm]')

%figure(3)
subplot(233)
plot(simHeightR,expHeightR,'.')
hold on 
line([0 simHeightR(end)],[0 simHeightR(end)])
xlabel('model height R [mm]')
ylabel('experiment height R [mm]')
subplot(236)
plot(simHeightR,resHeightR,'.')
hold on 
line([0 simHeightR(end)],[0 0])
ylim([-6 6])
xlabel('model height R [mm]')
ylabel('Residuals [mm]')
