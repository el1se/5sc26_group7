clearvars; clc;
load('results.mat')

%% rename data
% valve 13
% x4 is valve position valve 13
[t_1_020, x1_1_020, x3_1_020, x4_1_020] = cutdata(100,Data_1_020.Time,Data_1_020.x1_mm,Data_1_020.x3_mm,Data_1_020.x1x3ValvePos);
[t_1_040, x1_1_040, x3_1_040, x4_1_040] = cutdata(100,Data_1_040.Time,Data_1_040.x1_mm,Data_1_040.x3_mm,Data_1_040.x1x3ValvePos);
[t_1_060, x1_1_060, x3_1_060, x4_1_060] = cutdata(100,Data_1_060.Time,Data_1_060.x1_mm,Data_1_060.x3_mm,Data_1_060.x1x3ValvePos);
[t_1_080, x1_1_080, x3_1_080, x4_1_080] = cutdata(100,Data_1_080.Time,Data_1_080.x1_mm,Data_1_080.x3_mm,Data_1_080.x1x3ValvePos);
[t_1_100, x1_1_100, x3_1_100, x4_1_100] = cutdata(100,Data_1_100.Time,Data_1_100.x1_mm,Data_1_100.x3_mm,Data_1_100.x1x3ValvePos);

% % valve 32
% % x5 is valve position valve 32
[t_2_020, x2_2_020, x3_2_020, x5_2_020] = cutdata(90,Data_2_020.Time,Data_2_020.x2_mm,Data_2_020.x3_mm,Data_2_020.x3x2ValvePos);
[t_2_040, x2_2_040, x3_2_040, x5_2_040] = cutdata(90,Data_2_040.Time,Data_2_040.x2_mm,Data_2_040.x3_mm,Data_2_040.x3x2ValvePos);
[t_2_060, x2_2_060, x3_2_060, x5_2_060] = cutdata(90,Data_2_060.Time,Data_2_060.x2_mm,Data_2_060.x3_mm,Data_2_060.x3x2ValvePos);
[t_2_080, x2_2_080, x3_2_080, x5_2_080] = cutdata(90,Data_2_080.Time,Data_2_080.x2_mm,Data_2_080.x3_mm,Data_2_080.x3x2ValvePos);
[t_2_100, x2_2_100, x3_2_100, x5_2_100] = cutdata(90,Data_2_100.Time,Data_2_100.x2_mm,Data_2_100.x3_mm,Data_2_100.x3x2ValvePos);


valve13Means = [mean(x4_1_020) mean(x4_1_040) mean(x4_1_060) mean(x4_1_080) mean(x4_1_100)];
valve32Means = [mean(x5_2_020) mean(x5_2_040) mean(x5_2_060) mean(x5_2_080) mean(x5_2_100)];





%% plot:
figure(1)
title('Valve13 results')
subplot(1,3,1)
title('levels Tank 1')
hold off; hold on
plot(t_1_020,x1_1_020)
plot(t_1_040,x1_1_040)
plot(t_1_060,x1_1_060)
plot(t_1_080,x1_1_080)
plot(t_1_100,x1_1_100)
xlabel('time [s]')
ylabel('x_1 level[ mm]')

subplot(1,3,2)
title('levels Tank 3')
hold off; hold on
plot(t_1_020,x3_1_020)
plot(t_1_040,x3_1_040)
plot(t_1_060,x3_1_060)
plot(t_1_080,x3_1_080)
plot(t_1_100,x3_1_100)
xlabel('time [s]')
ylabel('x_3 level [mm]')

subplot(1,3,3)
title('valve 13 position')
hold off; hold on
plot(t_1_020,x4_1_020)
plot(t_1_040,x4_1_040)
plot(t_1_060,x4_1_060)
plot(t_1_080,x4_1_080)
plot(t_1_100,x4_1_100)
xlabel('time [s]')
ylabel('x_4 opening percentage [%]')



figure(2)
title('Valve32 results')
subplot(1,3,1)
title('levels Tank 2')
hold off; hold on
plot(t_2_020,x2_2_020)
plot(t_2_040,x2_2_040)
plot(t_2_060,x2_2_060)
plot(t_2_080,x2_2_080)
plot(t_2_100,x2_2_100)
xlabel('time [s]')
ylabel('x_2 level [mm]')

subplot(1,3,2)
title('levels Tank 3')
hold off; hold on
plot(t_2_020,x3_2_020)
plot(t_2_040,x3_2_040)
plot(t_2_060,x3_2_060)
plot(t_2_080,x3_2_080)
plot(t_2_100,x3_2_100)
xlabel('time [s]')
ylabel('x3 level [mm]')

subplot(1,3,3)
title('valve 32 position')
hold off; hold on
plot(t_2_020,x5_2_020)
plot(t_2_040,x5_2_040)
plot(t_2_060,x5_2_060)
plot(t_2_080,x5_2_080)
plot(t_2_100,x5_2_100)
xlabel('time [s]')
ylabel('x_5 opening percentage [%]')

%% processing
% say q13 flow from 1 to 3 is positive
% say q32 flow from 2 to 3 is positive

% for x1 and 20% open 
% dx1/dt = 4.25 for first 105 sec
% 0.0154*4.25 = c(x)*sqrt(2*g*(x1-x3))
c= 0.0154*4.25./sqrt(2*9.81*abs(x1_1_020-x3_1_020));
Q = mean(c).*sqrt(2*9.81*abs(x1_1_020-x3_1_020));
figure()
mean(c)
plot(Q)







%% simulation of ideal system
x0=[x1_1_100(1); 0];     % Initial condition (should satisfy the constraints)

options = odeset('RelTol',1e-6,'AbsTol',1e-8);
Tspan=[0 400];
  
[t,x]=ode45('BaseSys',Tspan,x0,options); % x is the state
t=t';
x=x';

figure(3)
hold off; hold on
% plot(t,x(1,:))
% hold on
% plot(t_1_100-t_1_100(1),x1_1_100)

%1_020 = c = 0.41
%1_040 = c = 6.5
%1_060 = c = 14
%1_080 = c = 15
%1_100 = c = 15.2

%2_020 = c = 0.14
%2_040 = c = 5.8
%2_060 = c = 13
%2_080 = c = 14.5
%2_100 = c = 15

c13 = [0.41 6.5 14 15 15.2];
c32 = [0.14 5.8 13 14.5 15];
plot(valve13Means,c13)
plot(valve32Means,c32)

poly_c13 = 0.1727*[0:0.01:100];
poly_c32 = 0.1657*[0:0.01:100];
plot([0:0.01:100],poly_c13)
plot([0:0.01:100],poly_c32)


xlabel('valve postition x_4 and x_5 [%]')
ylabel('c [-]')
legend('experiment c(x_4)','experiment c(x_5)','linear fit x_4','linear fit x_5')


%%
function [t, x1orx2, x3,valvepos] = cutdata(fromtime,t,x1orx2,x3, valvepos)
t=seconds(t);
x1orx2=double(x1orx2);
x3=double(x3);
valvepos=double(valvepos);

t=t(fromtime/0.01+1:end);
x1orx2=x1orx2(fromtime/0.01+1:end);
x3=x3(fromtime/0.01+1:end);
valvepos=valvepos(fromtime/0.01+1:end);

end

















