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

tstep = time(2)-time(1);
tend  = time(end) + tstep;
out   = sim('RL_TTS3_System_Blank.slx');

simHeightL = out.simHeightL;
simHeightM = out.simHeightM;
simHeightR = out.simHeightR;

resHeightL = (expHeightL - simHeightL)./sqrt(simHeightL);
resHeightM = (expHeightM - simHeightM)./sqrt(simHeightM);
resHeightR = (expHeightR - simHeightR)./sqrt(simHeightR);

figure(1)
subplot(121)
plot(simHeightL,expHeightL,'.')
hold on
line([0 simHeightL(end)],[0 simHeightL(end)])
xlabel('model height L [mm]')
ylabel('experiment height L [mm]')
subplot(122)
plot(simHeightL,resHeightL,'.')
hold on
line([0 simHeightL(end)],[0 0])
xlabel('model height L [mm]')
ylabel('standardized residuals [mm^{1/2}]')

figure(2)
subplot(121)
plot(simHeightM,expHeightM,'.')
hold on
line([0 simHeightM(end)],[0 simHeightM(end)])
xlabel('model height M [mm]')
ylabel('experiment height M [mm]')
subplot(122)
plot(simHeightM,resHeightM,'.')
hold on
line([0 simHeightM(end)],[0 0])
xlabel('model height M [mm]')
ylabel('standardized residuals [mm^{1/2}]')

figure(3)
subplot(121)
plot(simHeightR,expHeightR,'.')
hold on 
line([0 simHeightR(end)],[0 simHeightR(end)])
xlabel('model height R [mm]')
ylabel('experiment height R [mm]')
subplot(122)
plot(simHeightR,resHeightR,'.')
hold on 
line([0 simHeightR(end)],[0 0])
xlabel('model height R [mm]')
ylabel('standardized residuals [mm^{1/2}]')
