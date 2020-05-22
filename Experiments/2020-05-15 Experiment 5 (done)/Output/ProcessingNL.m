clearvars; close all; clc;
%% load input and reference and plot those
load('ToRemoteLabs4_results.mat');
load('input.mat');
Ts = 0.01;

% experiment:
ex_x1 = double(Data4.x1_mm);
ex_x2 = double(Data4.x2_mm);
ex_t = seconds(Data4.Time);% time
% simulation
sim_t = t4;
sim_x1 = cumsum(input4*0.01/0.0154);
for i = 2:length(input4)
    if(~SetValves_vector4(i-1) && (SetValves_vector4(i)))
        sim_x1(i:end)=sim_x1(i:end)-sim_x1(i);
    end
end
sim_x2 = sim_x1;
% plots
figure(1)
hold on
plot(ex_t,ex_x1);
plot(ex_t,ex_x2);
plot(sim_t,sim_x1);
plot(sim_t,sim_x2);
grid on
legend('x1','x2','x1_{ref}','x2_{ref}')

% now the data will be cut into sections and a correct offset will be added
% for training and validation
N_sec = 36000; % length of each section
sim_x1 = sim_x1(~SetValves_vector4(2:end));
sim_x1_train = sim_x1(1:length(sim_x1)*2/3);
sim_x1_val = sim_x1(length(sim_x1)*2/3+1:end);
sim_x2_train = sim_x1_train; sim_x2_val = sim_x1_val;
% make the experiment and simulation of the same length
ex_x1 = ex_x1(~[1 SetValves_vector4(2:end) ones(1,length(ex_x1)-length(SetValves_vector4))]);
ex_x2 = ex_x2(~[1 SetValves_vector4(2:end) ones(1,length(ex_x2)-length(SetValves_vector4))]);
% offset the data to 0
ex_x2(1*N_sec+1:2*N_sec) = ex_x2(1*N_sec+1:2*N_sec)-ex_x2(1*N_sec+1);
ex_x2(2*N_sec+1:3*N_sec) = ex_x2(2*N_sec+1:3*N_sec)-ex_x2(2*N_sec+1);
ex_x2(3*N_sec+1:4*N_sec) = ex_x2(3*N_sec+1:4*N_sec)-ex_x2(3*N_sec+1);
ex_x2(4*N_sec+1:5*N_sec) = ex_x2(4*N_sec+1:5*N_sec)-ex_x2(4*N_sec+1);
ex_x2(5*N_sec+1:6*N_sec) = ex_x2(5*N_sec+1:6*N_sec)-ex_x2(5*N_sec+1);
ex_x1_train = ex_x1(1:length(ex_x1)*2/3)';
ex_x1_val = ex_x1(length(ex_x1)*2/3+1:end)';
ex_x2_train = ex_x2(1:length(ex_x2)*2/3)';
ex_x2_val = ex_x2(length(ex_x2)*2/3+1:end)';

input = input4(~SetValves_vector4(2:end));
input_val = input(length(input)*2/3+1:end);
input_train = input(1:length(input)*2/3);
t_train = Ts:Ts:length(input_train)*Ts;
t_val = Ts:Ts:length(input_val)*Ts;

figure()
hold on
plot(Ts:Ts:length(sim_x1)*Ts,sim_x1)
hold on
plot(Ts:Ts:length(ex_x1)*Ts,ex_x1)
hold on
plot(Ts:Ts:length(ex_x1)*Ts,input*1000)
legend('simx1','experimentx1','input*1000')
% looks like a 1 second delay everywhere

%% model identification and validation for entire training set
c_x1_InitGues = [1 1 1 1 1 1];
c_x2_InitGues = [1 1 1 1 1 1];

% x1
options = optimoptions('fmincon','Display','off');
fun = @(c) CostFunction(t_train,input_train,ex_x1_train,c);
c_x1_optimum = fmincon(fun,c_x1_InitGues,[],[],[],[],[],[],[],options);

% plot training to see if makes sense correct;
N = (size(input_train,2))/20;
C = [repelem(c_x1_optimum(2:6),1,N) fliplr(repelem(c_x1_optimum(2:6),1,N)) repelem(c_x1_optimum(2:6),1,N) fliplr(repelem(c_x1_optimum(2:6),1,N))];
fitted_train_input = input_train.*C;
fitted_train_tx1 = t_train+C(1);
N = (size(input_train,2))/4;
fitted_train_x1 = [cumsum(fitted_train_input(1:N)) cumsum(fitted_train_input(N+1:N*2))...
    cumsum(fitted_train_input(1+N*2:N*3)) cumsum(fitted_train_input(N*3+1:4*N))]*0.01/0.0154;

% now test on valication set
C=C(1:length(ex_x1_val));
fitted_val_input = input_val.*C;
N = (size(fitted_val_input,2))/2;
fitted_val_x1 = [cumsum(fitted_val_input(1:N)) cumsum(fitted_val_input(N+1:N*2))]*0.01/0.0154;
fitted_val_tx1 = t_val+C(1);

%x2
fun = @(c) CostFunction(t_train,input_train,ex_x2_train,c);
c_x2_optimum = fmincon(fun,c_x2_InitGues,[],[],[],[],[],[],[],options);

% plot training to see if makes sense correct;
N = (size(input_train,2))/20;
C = [repelem(c_x2_optimum(2:6),1,N) fliplr(repelem(c_x2_optimum(2:6),1,N)) repelem(c_x2_optimum(2:6),1,N) fliplr(repelem(c_x2_optimum(2:6),1,N))];
fitted_train_input = input_train.*C;
fitted_train_tx2 = t_train+C(1);
N = (size(input_train,2))/4;
fitted_train_x2 = [cumsum(fitted_train_input(1:N)) cumsum(fitted_train_input(N+1:N*2))...
    cumsum(fitted_train_input(1+N*2:N*3)) cumsum(fitted_train_input(N*3+1:4*N))]*0.01/0.0154;

% now test on valication set
C=C(1:length(ex_x2_val));
fitted_val_input = input_val.*C;
N = (size(fitted_val_input,2))/2;
fitted_val_x2 = [cumsum(fitted_val_input(1:N)) cumsum(fitted_val_input(N+1:N*2))]*0.01/0.0154;
fitted_val_tx2 = t_val+C(1);

%% plot results
figure()
subplot(2,3,[1 2])
hold on
plot(t_train,ex_x1_train)
plot(fitted_train_tx1,fitted_train_x1)
plot(t_train,sim_x1_train)
legend('expeirment','fitted','unfitted')
title('pump1-x1 training'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,3)
hold on
plot(t_val,ex_x1_val)
plot(fitted_val_tx1,fitted_val_x1)
plot(t_val,sim_x1_val)
legend('expeirment','fitted','unfitted')
title('pump1-x1 validation'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,[4 5])
hold on
plot(t_train,ex_x2_train)
plot(fitted_train_tx2,fitted_train_x2)
plot(t_train,sim_x2_train)
legend('expeirment','fitted','unfitted')
title('pump1-x2 training'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,6)
hold on
plot(t_val,ex_x2_val)
plot(fitted_val_tx2,fitted_val_x2)
plot(t_val,sim_x2_val)
legend('expeirment','fitted','unfitted')
title('pump1-x2 validation'); ylabel('height [mm]');xlabel('time [s]')
grid on

%% generate function based on fit
figure()
u_values = 0.02:0.02:0.1;
subplot(1,2,1)
plot(u_values,c_x1_optimum(2:6))
hold on
p = polyfit(0.02:0.02:0.1,c_x1_optimum(2:6),3);
u_values = 0:0.001:0.1;
yfit = polyval(p,u_values);
plot(u_values,yfit)
title('pump 1 fit')
subplot(1,2,2)
u_values = 0.02:0.02:0.1;
plot(u_values,c_x2_optimum(2:6))
hold on
p = polyfit(0.02:0.02:0.1,c_x2_optimum(2:6),3);
u_values = 0:0.001:0.1;
yfit = polyval(p,u_values);
plot(u_values,yfit)
title('pump 2 fit')




%% functions
function [V] = CostFunction(t,u,x1,c)
N = (size(u,2))/20;

% optimize input
% C = repelem(repelem(c(2:6),1,N),1,4);
C = [repelem(c(2:6),1,N) fliplr(repelem(c(2:6),1,N)) repelem(c(2:6),1,N) fliplr(repelem(c(2:6),1,N))];
u = u.*C;
% figure()
% plot(C)
% offset with a delay time
t_fitted = t+c(1);
N = (size(u,2))/4;

%make fit without time adjustment
x1_fitted = [cumsum(u(1:N)) cumsum(u(N+1:N*2)) cumsum(u(1+N*2:N*3)) cumsum(u(N*3+1:4*N))]*0.01/0.0154;

% interpolate between new timeoffset to put in cost function
x1_fitted = interp1(t_fitted,x1_fitted,t);
% the beginnin contains NaN. replace those with zeros
x1_fitted(isnan(x1_fitted))=0;

% any(isnan(ref_fitted))

V = norm(x1_fitted-x1,2);
end






