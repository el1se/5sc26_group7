clearvars; close all; clc;
%% load input and reference and plot those
load('ToRemoteLabs4_results.mat');
load('ToRemoteLabs6_results.mat');
load('input.mat');
Ts = 0.01;

% experiment:
ex4_x1 = double(Data4.x1_mm);
ex4_x2 = double(Data4.x2_mm);
ex4_t = seconds(Data4.Time);% time
ex6_x1 = double(Data6.x1_mm);
ex6_x2 = double(Data6.x2_mm);
ex6_t = seconds(Data6.Time);% time

% simulation
sim4_t = t4;
sim4_x1 = cumsum(input4*0.01/0.0154);
sim6_t = t6;
sim6_x1 = cumsum(input6*0.01/0.0154);
for i = 2:length(input4)
    if(~SetValves_vector4(i-1) && (SetValves_vector4(i)))
        sim4_x1(i:end)=sim4_x1(i:end)-sim4_x1(i);
    end
end
for i = 2:length(input6)
    if(~SetValves_vector6(i-1) && (SetValves_vector6(i)))
        sim6_x1(i:end)=sim6_x1(i:end)-sim6_x1(i);
    end
end
sim4_x2 = sim4_x1; sim6_x2 = sim6_x1;

% plots
figure(1)
subplot(1,2,1)
hold on
plot(ex4_t,ex4_x1);
plot(ex4_t,ex4_x2);
plot(sim4_t,sim4_x1);
plot(sim4_t,sim4_x2);
title('experiment 4')
grid on
legend('x1','x2','x1_{ref}','x2_{ref}')
subplot(1,2,2)
hold on
plot(ex6_t,ex6_x1);
plot(ex6_t,ex6_x2);
plot(sim6_t,sim6_x1);
plot(sim6_t,sim6_x2);
title('experiment 6')
grid on
legend('x1','x2','x1_{ref}','x2_{ref}')


% now the data will be cut into sections and a correct offset will be added
% for training and validation
N_sec = 36000; % length of each section
sim4_x1 = sim4_x1(~SetValves_vector4(2:end));
sim4_x1_train = sim4_x1(1:length(sim4_x1)*2/3);
sim4_x1_val = sim4_x1(length(sim4_x1)*2/3+1:end);
sim4_x2_train = sim4_x1_train; sim4_x2_val = sim4_x1_val;
% make the experiment and simulation of the same length
ex4_x1 = ex4_x1(~[1 SetValves_vector4(2:end) ones(1,length(ex4_x1)-length(SetValves_vector4))]);
ex4_x2 = ex4_x2(~[1 SetValves_vector4(2:end) ones(1,length(ex4_x2)-length(SetValves_vector4))]);
% offset the data to 0
ex4_x2(1*N_sec+1:2*N_sec) = ex4_x2(1*N_sec+1:2*N_sec)-ex4_x2(1*N_sec+1);
ex4_x2(2*N_sec+1:3*N_sec) = ex4_x2(2*N_sec+1:3*N_sec)-ex4_x2(2*N_sec+1);
ex4_x2(3*N_sec+1:4*N_sec) = ex4_x2(3*N_sec+1:4*N_sec)-ex4_x2(3*N_sec+1);
ex4_x2(4*N_sec+1:5*N_sec) = ex4_x2(4*N_sec+1:5*N_sec)-ex4_x2(4*N_sec+1);
ex4_x2(5*N_sec+1:6*N_sec) = ex4_x2(5*N_sec+1:6*N_sec)-ex4_x2(5*N_sec+1);
ex4_x1_train = ex4_x1(1:length(ex4_x1)*2/3)';
ex4_x1_val = ex4_x1(length(ex4_x1)*2/3+1:end)';
ex4_x2_train = ex4_x2(1:length(ex4_x2)*2/3)';
ex4_x2_val = ex4_x2(length(ex4_x2)*2/3+1:end)';

input4 = input4(~SetValves_vector4(2:end));
input4_val = input4(length(input4)*2/3+1:end);
input4_train = input4(1:length(input4)*2/3);
t4_train = Ts:Ts:length(input4_train)*Ts;
t4_val = Ts:Ts:length(input4_val)*Ts;

sim6_x1 = sim6_x1(~SetValves_vector6(2:end));
sim6_x1_train = sim6_x1(1:length(sim6_x1)*2/3);
sim6_x1_val = sim6_x1(length(sim6_x1)*2/3+1:end);
sim6_x2_train = sim6_x1_train; sim6_x2_val = sim6_x1_val;
% make the experiment and simulation of the same length
ex6_x1 = ex6_x1(~[1 SetValves_vector6(2:end) ones(1,length(ex6_x1)-length(SetValves_vector6))]);
ex6_x2 = ex6_x2(~[1 SetValves_vector6(2:end) ones(1,length(ex6_x2)-length(SetValves_vector6))]);
% offset the data to 0
ex6_x2(1*N_sec+1:2*N_sec) = ex6_x2(1*N_sec+1:2*N_sec)-ex6_x2(1*N_sec+1);
ex6_x2(2*N_sec+1:3*N_sec) = ex6_x2(2*N_sec+1:3*N_sec)-ex6_x2(2*N_sec+1);
ex6_x2(3*N_sec+1:4*N_sec) = ex6_x2(3*N_sec+1:4*N_sec)-ex6_x2(3*N_sec+1);
ex6_x2(4*N_sec+1:5*N_sec) = ex6_x2(4*N_sec+1:5*N_sec)-ex6_x2(4*N_sec+1);
ex6_x2(5*N_sec+1:6*N_sec) = ex6_x2(5*N_sec+1:6*N_sec)-ex6_x2(5*N_sec+1);
ex6_x1_train = ex6_x1(1:length(ex6_x1)*2/3)';
ex6_x1_val = ex6_x1(length(ex6_x1)*2/3+1:end)';
ex6_x2_train = ex6_x2(1:length(ex6_x2)*2/3)';
ex6_x2_val = ex6_x2(length(ex6_x2)*2/3+1:end)';

input6 = input6(~SetValves_vector6(2:end));
input6_val = input6(length(input6)*2/3+1:end);
input6_train = input6(1:length(input6)*2/3);
t6_train = Ts:Ts:length(input6_train)*Ts;
t6_val = Ts:Ts:length(input6_val)*Ts;

% plot
figure(2)
subplot(1,2,1)
hold on
plot(Ts:Ts:length(sim4_x1)*Ts,sim4_x1)
hold on
plot(Ts:Ts:length(ex4_x1)*Ts,ex4_x1)
hold on
plot(Ts:Ts:length(ex4_x1)*Ts,input4*1000)
legend('simx1','experimentx1','input*1000')
title('experiment 4')
% looks like a 1 second delay everywhere
subplot(1,2,2)
hold on
plot(Ts:Ts:length(sim6_x1)*Ts,sim6_x1)
hold on
plot(Ts:Ts:length(ex6_x1)*Ts,ex6_x1)
hold on
plot(Ts:Ts:length(ex6_x1)*Ts,input6*1000)
legend('simx1','experimentx1','input*1000')
title('experiment 6')

%% model identification and validation for entire training set4
c4_x1_InitGues = [1 1 1 1 1 1];
c4_x2_InitGues = [1 1 1 1 1 1];

% x1
options = optimoptions('fmincon','Display','off');
fun = @(c) CostFunction4(t4_train,input4_train,ex4_x1_train,c);
c4_x1_optimum = fmincon(fun,c4_x1_InitGues,[],[],[],[],[],[],[],options);

% plot training to see if makes sense correct;
N = (size(input4_train,2))/20;
C4 = [repelem(c4_x1_optimum(2:6),1,N) fliplr(repelem(c4_x1_optimum(2:6),1,N)) repelem(c4_x1_optimum(2:6),1,N) fliplr(repelem(c4_x1_optimum(2:6),1,N))];
fitted4_train_input = input4_train.*C4;
fitted4_train_tx1 = t4_train+C4(1);
N = (size(input4_train,2))/4;
fitted4_train_x1 = [cumsum(fitted4_train_input(1:N)) cumsum(fitted4_train_input(N+1:N*2))...
    cumsum(fitted4_train_input(1+N*2:N*3)) cumsum(fitted4_train_input(N*3+1:4*N))]*0.01/0.0154;

% now test on valication set
C4=C4(1:length(ex4_x1_val));
fitted4_val_input = input4_val.*C4;
N = (size(fitted4_val_input,2))/2;
fitted4_val_x1 = [cumsum(fitted4_val_input(1:N)) cumsum(fitted4_val_input(N+1:N*2))]*0.01/0.0154;
fitted4_val_tx1 = t4_val+C4(1);

%x2
fun = @(c) CostFunction4(t4_train,input4_train,ex4_x2_train,c);
c4_x2_optimum = fmincon(fun,c4_x2_InitGues,[],[],[],[],[],[],[],options);

% plot training to see if makes sense correct;
N = (size(input4_train,2))/20;
C4 = [repelem(c4_x2_optimum(2:6),1,N) fliplr(repelem(c4_x2_optimum(2:6),1,N)) repelem(c4_x2_optimum(2:6),1,N) fliplr(repelem(c4_x2_optimum(2:6),1,N))];
fitted4_train_input = input4_train.*C4;
fitted4_train_tx2 = t4_train+C4(1);
N = (size(input4_train,2))/4;
fitted4_train_x2 = [cumsum(fitted4_train_input(1:N)) cumsum(fitted4_train_input(N+1:N*2))...
    cumsum(fitted4_train_input(1+N*2:N*3)) cumsum(fitted4_train_input(N*3+1:4*N))]*0.01/0.0154;

% now test on valication set
C4=C4(1:length(ex4_x2_val));
fitted4_val_input = input4_val.*C4;
N = (size(fitted4_val_input,2))/2;
fitted4_val_x2 = [cumsum(fitted4_val_input(1:N)) cumsum(fitted4_val_input(N+1:N*2))]*0.01/0.0154;
fitted4_val_tx2 = t4_val+C4(1);

%% model identification and validation for entire training set6
c6_x1_InitGues = [1 1 1 1 1 1 1 1 1 1 1];
c6_x2_InitGues = [1 1 1 1 1 1 1 1 1 1 1];

% x1
options = optimoptions('fmincon','Display','off');
fun = @(c) CostFunction6(t6_train,input6_train,ex6_x1_train,c);
c6_x1_optimum = fmincon(fun,c6_x1_InitGues,[],[],[],[],[],[],[],options);

% plot training to see if makes sense correct;
N = (size(input6_train,2))/40;
C6 = [repelem(c6_x1_optimum(2:11),1,N) fliplr(repelem(c6_x1_optimum(2:11),1,N)) repelem(c6_x1_optimum(2:11),1,N) fliplr(repelem(c6_x1_optimum(2:11),1,N))];
fitted6_train_input = input6_train.*C6;
fitted6_train_tx1 = t6_train+C6(1);
N = (size(input6_train,2))/4;
fitted6_train_x1 = [cumsum(fitted6_train_input(1:N)) cumsum(fitted6_train_input(N+1:N*2))...
    cumsum(fitted6_train_input(1+N*2:N*3)) cumsum(fitted6_train_input(N*3+1:4*N))]*0.01/0.0154;

% now test on valication set
C6=C6(1:length(ex6_x1_val));
fitted6_val_input = input6_val.*C6;
N = (size(fitted6_val_input,2))/2;
fitted6_val_x1 = [cumsum(fitted6_val_input(1:N)) cumsum(fitted6_val_input(N+1:N*2))]*0.01/0.0154;
fitted6_val_tx1 = t6_val+C6(1);

%x2
fun = @(c) CostFunction6(t6_train,input6_train,ex6_x2_train,c);
c6_x2_optimum = fmincon(fun,c6_x2_InitGues,[],[],[],[],[],[],[],options);

% plot training to see if makes sense correct;
N = (size(input6_train,2))/40;
C6 = [repelem(c6_x2_optimum(2:11),1,N) fliplr(repelem(c6_x2_optimum(2:11),1,N)) repelem(c6_x2_optimum(2:11),1,N) fliplr(repelem(c6_x2_optimum(2:11),1,N))];
fitted6_train_input = input6_train.*C6;
fitted6_train_tx2 = t6_train+C6(1);
N = (size(input6_train,2))/4;
fitted6_train_x2 = [cumsum(fitted6_train_input(1:N)) cumsum(fitted6_train_input(N+1:N*2))...
    cumsum(fitted6_train_input(1+N*2:N*3)) cumsum(fitted6_train_input(N*3+1:4*N))]*0.01/0.0154;

% now test on valication set
C6=C6(1:length(ex6_x2_val));
fitted6_val_input = input6_val.*C6;
N = (size(fitted6_val_input,2))/2;
fitted6_val_x2 = [cumsum(fitted6_val_input(1:N)) cumsum(fitted6_val_input(N+1:N*2))]*0.01/0.0154;
fitted6_val_tx2 = t6_val+C6(1);

%% model identification and validation for bot sets 4 and 6
c46_x1_InitGues = [1 1 1 1 1 1 1 1 1 1 1];
c46_x2_InitGues = [1 1 1 1 1 1 1 1 1 1 1];

% x1
options = optimoptions('fmincon','Display','off');
t46_train = Ts:Ts:(length(input4_train)+length(input6_train))*Ts;
t46_val = Ts:Ts:(length(input4_val)+length(input6_val))*Ts;
input46_train = [input4_train input6_train];
ex46_x1_train = [ex4_x1_train ex6_x1_train];
ex46_x1_val = [ex4_x1_val ex6_x1_val];
input46_val =[input4_val input6_val];
sim46_x1_train = [sim4_x1_train sim6_x1_train];
sim46_x1_val = [sim4_x1_val sim6_x1_val];

fun = @(c) CostFunction46(t46_train,input46_train,ex46_x1_train,c);
c46_x1_optimum = fmincon(fun,c46_x1_InitGues,[],[],[],[],[],[],[],options);

% plot training to see if makes sense correct;
N = (size(input46_train,2)/2)/20;
C46 = [repelem(c46_x1_optimum(3:2:11),1,N) fliplr(repelem(c46_x1_optimum(3:2:11),1,N))...
    repelem(c46_x1_optimum(3:2:11),1,N) fliplr(repelem(c46_x1_optimum(3:2:11),1,N))];
N = (size(input46_train,2)/2)/40;
C46 = [C46 repelem(c46_x1_optimum(2:11),1,N) fliplr(repelem(c46_x1_optimum(2:11),1,N))...
    repelem(c46_x1_optimum(2:11),1,N) fliplr(repelem(c46_x1_optimum(2:11),1,N))];
fitted46_train_input = input46_train.*C46;
fitted46_train_tx1 = t46_train+C46(1);
N = (size(input46_train,2))/8;
fitted46_train_x1 = [cumsum(fitted46_train_input(1:N)) cumsum(fitted46_train_input(N+1:N*2))...
    cumsum(fitted46_train_input(1+N*2:N*3)) cumsum(fitted46_train_input(N*3+1:4*N))...
    cumsum(fitted46_train_input(N*4+1:5*N)) cumsum(fitted46_train_input(N*5+1:6*N))...
    cumsum(fitted46_train_input(N*6+1:7*N)) cumsum(fitted46_train_input(N*7+1:8*N))]*0.01/0.0154;

% now test on valication set
C46=C46(1:length(ex46_x1_val));
fitted46_val_input = input46_val.*C46;
N = (size(fitted46_val_input,2))/4;
fitted46_val_x1 = [cumsum(fitted46_val_input(1:N)) cumsum(fitted46_val_input(N+1:N*2))...
    cumsum(fitted46_val_input(N*2+1:N*3)) cumsum(fitted46_val_input(N*3+1:N*4))]*0.01/0.0154;
fitted46_val_tx1 = t46_val+C46(1);

% x2
ex46_x2_train = [ex4_x2_train ex6_x2_train];
ex46_x2_val = [ex4_x2_val ex6_x2_val];
sim46_x2_train = [sim4_x2_train sim6_x2_train];
sim46_x2_val = [sim4_x2_val sim6_x2_val];

fun = @(c) CostFunction46(t46_train,input46_train,ex46_x2_train,c);
c46_x2_optimum = fmincon(fun,c46_x2_InitGues,[],[],[],[],[],[],[],options);

% plot training to see if makes sense correct;
N = (size(input46_train,2)/2)/20;
C46 = [repelem(c46_x2_optimum(3:2:11),1,N) fliplr(repelem(c46_x2_optimum(3:2:11),1,N))...
    repelem(c46_x2_optimum(3:2:11),1,N) fliplr(repelem(c46_x2_optimum(3:2:11),1,N))];
N = (size(input46_train,2)/2)/40;
C46 = [C46 repelem(c46_x2_optimum(2:11),1,N) fliplr(repelem(c46_x2_optimum(2:11),1,N))...
    repelem(c46_x2_optimum(2:11),1,N) fliplr(repelem(c46_x2_optimum(2:11),1,N))];
fitted46_train_input = input46_train.*C46;
fitted46_train_tx2 = t46_train+C46(1);
N = (size(input46_train,2))/8;
fitted46_train_x2 = [cumsum(fitted46_train_input(1:N)) cumsum(fitted46_train_input(N+1:N*2))...
    cumsum(fitted46_train_input(1+N*2:N*3)) cumsum(fitted46_train_input(N*3+1:4*N))...
    cumsum(fitted46_train_input(N*4+1:5*N)) cumsum(fitted46_train_input(N*5+1:6*N))...
    cumsum(fitted46_train_input(N*6+1:7*N)) cumsum(fitted46_train_input(N*7+1:8*N))]*0.01/0.0154;

% now test on valication set
C46=C46(1:length(ex46_x2_val));
fitted46_val_input = input46_val.*C46;
N = (size(fitted46_val_input,2))/4;
fitted46_val_x2 = [cumsum(fitted46_val_input(1:N)) cumsum(fitted46_val_input(N+1:N*2))...
    cumsum(fitted46_val_input(N*2+1:N*3)) cumsum(fitted46_val_input(N*3+1:N*4))]*0.01/0.0154;
fitted46_val_tx2 = t46_val+C46(1);


%% plot results for set4

figure(3)
subplot(2,3,[1 2])
hold on
plot(t4_train,ex4_x1_train)
plot(fitted4_train_tx1,fitted4_train_x1)
plot(t4_train,sim4_x1_train)
legend('expeirment','fitted','unfitted')
title('experiment 4 pump1-x1 training'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,3)
hold on
plot(t4_val,ex4_x1_val)
plot(fitted4_val_tx1,fitted4_val_x1)
plot(t4_val,sim4_x1_val)
legend('expeirment','fitted','unfitted')
title('experiment 4 pump1-x1 validation'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,[4 5])
hold on
plot(t4_train,ex4_x2_train)
plot(fitted4_train_tx2,fitted4_train_x2)
plot(t4_train,sim4_x2_train)
legend('expeirment','fitted','unfitted')
title('experiment 4 pump1-x2 training'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,6)
hold on
plot(t4_val,ex4_x2_val)
plot(fitted4_val_tx2,fitted4_val_x2)
plot(t4_val,sim4_x2_val)
legend('expeirment','fitted','unfitted')
title('expeirment 4 pump1-x2 validation'); ylabel('height [mm]');xlabel('time [s]')
grid on

%% plot results for set6
figure(4)
subplot(2,3,[1 2])
hold on
plot(t6_train,ex6_x1_train)
plot(fitted6_train_tx1,fitted6_train_x1)
plot(t6_train,sim6_x1_train)
legend('expeirment','fitted','unfitted')
title('experiment 6 pump1-x1 training'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,3)
hold on
plot(t6_val,ex6_x1_val)
plot(fitted6_val_tx1,fitted6_val_x1)
plot(t6_val,sim6_x1_val)
legend('expeirment','fitted','unfitted')
title('experiment 6 pump1-x1 validation'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,[4 5])
hold on
plot(t6_train,ex6_x2_train)
plot(fitted6_train_tx2,fitted6_train_x2)
plot(t6_train,sim6_x2_train)
legend('expeirment','fitted','unfitted')
title('experiment 6 pump1-x2 training'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,6)
hold on
plot(t6_val,ex6_x2_val)
plot(fitted6_val_tx2,fitted6_val_x2)
plot(t6_val,sim6_x2_val)
legend('expeirment','fitted','unfitted')
title('expeirment 6 pump1-x2 validation'); ylabel('height [mm]');xlabel('time [s]')
grid on

%% plot results for set46
figure(5)
subplot(2,3,[1 2])
hold on
plot(t46_train,ex46_x1_train)
plot(fitted46_train_tx1,fitted46_train_x1)
plot(t46_train,sim46_x1_train)
legend('expeirment','fitted','unfitted')
title('experiment 46 pump1-x1 training'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,3)
hold on
plot(t46_val,ex46_x1_val)
plot(fitted46_val_tx1,fitted46_val_x1)
plot(t46_val,sim46_x1_val)
legend('expeirment','fitted','unfitted')
title('experiment 46 pump1-x1 validation'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,[4 5])
hold on
plot(t46_train,ex46_x2_train)
plot(fitted46_train_tx2,fitted46_train_x2)
plot(t46_train,sim46_x2_train)
legend('expeirment','fitted','unfitted')
title('experiment 46 pump1-x2 training'); ylabel('height [mm]');xlabel('time [s]')
grid on
subplot(2,3,6)
hold on
plot(t46_val,ex46_x2_val)
plot(fitted46_val_tx2,fitted46_val_x2)
plot(t46_val,sim46_x2_val)
legend('expeirment','fitted','unfitted')
title('expeirment 46 pump1-x2 validation'); ylabel('height [mm]');xlabel('time [s]')
grid on


%% generate function based on fit
figure(6)
u_values = 0.02:0.02:0.1;
subplot(1,2,1)
plot(u_values,c4_x1_optimum(2:6))
hold on
p = polyfit(u_values,c4_x1_optimum(2:6),2);
u_values = 0:0.001:0.1;
yfit = polyval(p,u_values);
plot(u_values,yfit)
u_values = 0.02:0.02:0.1;
net = newrb(u_values,c4_x1_optimum(2:6),1,4);
u_values = 0:0.001:0.1;
Y = net(u_values);
plot(u_values,Y)
title('experiment 4 pump 1 fit')
subplot(1,2,2)
u_values = 0.02:0.02:0.1;
plot(u_values,c4_x2_optimum(2:6))
hold on
p = polyfit(u_values,c4_x2_optimum(2:6),2);
u_values = 0:0.001:0.1;
yfit = polyval(p,u_values);
plot(u_values,yfit)
u_values = 0.02:0.02:0.1;
net = newrb(u_values,c4_x2_optimum(2:6),1,4);
u_values = 0:0.001:0.1;
Y = net(u_values);
plot(u_values,Y)

title('experiment 4 pump 2 fit')

figure(7)
u_values = 0.01:0.01:0.1;
subplot(1,2,1)
plot(u_values,c6_x1_optimum(2:11))
hold on
p = polyfit(u_values,c6_x1_optimum(2:11),2);
u_values = 0:0.001:0.1;
yfit = polyval(p,u_values);
plot(u_values,yfit)
u_values = 0.01:0.01:0.1;
net = newrb(u_values,c6_x1_optimum(2:11),1,4);
u_values = 0:0.001:0.1;
Y = net(u_values);
plot(u_values,Y)
title('experiment 6 pump 1 fit')
subplot(1,2,2)
u_values = 0.01:0.01:0.1;
plot(u_values,c6_x2_optimum(2:11))
hold on
p = polyfit(u_values,c6_x2_optimum(2:11),2);
u_values = 0:0.001:0.1;
yfit = polyval(p,u_values);
plot(u_values,yfit)
u_values = 0.01:0.01:0.1;
net = newrb(u_values,c6_x2_optimum(2:11),1,4);
u_values = 0:0.001:0.1;
Y = net(u_values);
plot(u_values,Y)

title('experiment 6 pump 2 fit')

figure(8)
u_values = 0.01:0.01:0.1;
subplot(1,2,1)
plot(u_values,c46_x1_optimum(2:11),'+')
hold on
p = polyfit(u_values,c46_x1_optimum(2:11),2);
u_values = 0:0.001:0.1;
yfit = polyval(p,u_values);
plot(u_values,yfit)
u_values = 0.01:0.01:0.1;
net = newrb(u_values,c46_x1_optimum(2:11),1,4);
u_values = 0:0.001:0.1;
Y = net(u_values);
plot(u_values,Y)
title('experiment 46 pump 1 fit')
legend('fit fmincon','fitted 2nd poly','fitted GP')
xlabel('u [l/s]'); ylabel('correction factor c')
subplot(1,2,2)
u_values = 0.01:0.01:0.1;
plot(u_values,c46_x2_optimum(2:11),'+')
hold on
p = polyfit(u_values,c46_x2_optimum(2:11),2);
u_values = 0:0.001:0.1;
yfit = polyval(p,u_values);
plot(u_values,yfit)
title('experiment 46 pump 2 fit')
u_values = 0.01:0.01:0.1;
net = newrb(u_values,c46_x2_optimum(2:11),1,4);
u_values = 0:0.001:0.1;
Y = net(u_values);
plot(u_values,Y)
legend('fit fmincon','fitted 2nd poly','fitted GP')
xlabel('u [l/s]'); ylabel('correction factor c')

% add the new found stuff to existing graphs
figure(5)
subplot(2,3,[1 2])
hold on
% p46_x1 = polyfit(0.01:0.01:0.1,c46_x1_optimum(2:11),2);
% p46_x2 = polyfit(0.01:0.01:0.1,c46_x2_optimum(2:11),2);
fittedfit46_train_input_x1 = Pump1Function(input46_train);
% t46_train
N = (size(input46_train,2))/8;
fittedfit46_train_x1 = [cumsum(fittedfit46_train_input_x1(1:N)) cumsum(fittedfit46_train_input_x1(N+1:N*2))...
    cumsum(fittedfit46_train_input_x1(1+N*2:N*3)) cumsum(fittedfit46_train_input_x1(N*3+1:4*N))...
    cumsum(fittedfit46_train_input_x1(N*4+1:5*N)) cumsum(fittedfit46_train_input_x1(N*5+1:6*N))...
    cumsum(fittedfit46_train_input_x1(N*6+1:7*N)) cumsum(fittedfit46_train_input_x1(N*7+1:8*N))]*0.01/0.0154;
plot(t46_train,fittedfit46_train_x1)
legend('expeirment','fitted','unfitted','fitted fit')

subplot(2,3,[4 5])
hold on
fittedfit46_train_input_x2 = Pump2Function(input46_train);
fittedfit46_train_x2 = [cumsum(fittedfit46_train_input_x2(1:N)) cumsum(fittedfit46_train_input_x2(N+1:N*2))...
    cumsum(fittedfit46_train_input_x2(1+N*2:N*3)) cumsum(fittedfit46_train_input_x2(N*3+1:4*N))...
    cumsum(fittedfit46_train_input_x2(N*4+1:5*N)) cumsum(fittedfit46_train_input_x2(N*5+1:6*N))...
    cumsum(fittedfit46_train_input_x2(N*6+1:7*N)) cumsum(fittedfit46_train_input_x2(N*7+1:8*N))]*0.01/0.0154;
plot(t46_train,fittedfit46_train_x2)
legend('expeirment','fitted','unfitted','fitted fit')

subplot(2,3,3)
hold on
fittedfit46_val_input_x1 = Pump1Function(input46_val);
N = (size(input46_val,2))/4;
fittedfit46_val_x1 = [cumsum(fittedfit46_val_input_x1(1:N)) cumsum(fittedfit46_val_input_x1(N+1:N*2))...
    cumsum(fittedfit46_val_input_x1(N*2+1:N*3)) cumsum(fittedfit46_val_input_x1(N*3+1:N*4))]*0.01/0.0154;
plot(t46_val,fittedfit46_val_x1)
legend('expeirment','fitted','unfitted','fitted fit')

subplot(2,3,6)
hold on
fittedfit46_val_input_x2 = Pump2Function(input46_val);
N = (size(input46_val,2))/4;
fittedfit46_val_x2 = [cumsum(fittedfit46_val_input_x2(1:N)) cumsum(fittedfit46_val_input_x2(N+1:N*2))...
    cumsum(fittedfit46_val_input_x2(N*2+1:N*3)) cumsum(fittedfit46_val_input_x2(N*3+1:N*4))]*0.01/0.0154;
plot(t46_val,fittedfit46_val_x2)
legend('expeirment','fitted','unfitted','fitted fit')

%% figures for report and SSA
% training pump 1
figure(9)
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
plot(t46_train,ex46_x1_train)
plot(t46_train,sim46_x1_train)
plot(fitted46_train_tx1,fitted46_train_x1)
plot(t46_train,fittedfit46_train_x1)
ylabel('x_1 [mm]');xlabel('time [s]')
legend('expirement','theoretic','fitted','fitted polynomial')
grid on

% validation pump 1
figure(10)
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
plot(t46_val,ex46_x1_val)
plot(t46_val,sim46_x1_val)
plot(fitted46_val_tx1,fitted46_val_x1)
plot(t46_val,fittedfit46_val_x1)
ylabel('x_1 [mm]');xlabel('time [s]')
legend('expirement','theoretic','fitted','fitted polynomial')
grid on

% training pump 2
figure(11)
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
plot(t46_train,ex46_x2_train)
plot(t46_train,sim46_x2_train)
plot(fitted46_train_tx2,fitted46_train_x2)
plot(t46_train,fittedfit46_train_x2)
ylabel('x_2 [mm]');xlabel('time [s]')
legend('expirement','theoretic','fitted','fitted polynomial')
grid on

% validation pump 2
figure(12)
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
plot(t46_val,ex46_x2_val)
plot(t46_val,sim46_x2_val)
plot(fitted46_val_tx2,fitted46_val_x2)
plot(t46_val,fittedfit46_val_x2)
ylabel('x_2 [mm]');xlabel('time [s]')
legend('expirement','theoretic','fitted','fitted polynomial')
grid on




%% functions
function [V] = CostFunction4(t,u,x1,c)
N = (size(u,2))/20;

% optimize input
% C = repelem(repelem(c(2:6),1,N),1,4);
C = [repelem(c(2:6),1,N) fliplr(repelem(c(2:6),1,N)) repelem(c(2:6),1,N) fliplr(repelem(c(2:6),1,N))];
u = u.*C;

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

function [V] = CostFunction6(t,u,x1,c)
N = (size(u,2))/40;

% optimize input
% C = repelem(repelem(c(2:6),1,N),1,4);
C = [repelem(c(2:11),1,N) fliplr(repelem(c(2:11),1,N)) repelem(c(2:11),1,N) fliplr(repelem(c(2:11),1,N))];
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

function [V] = CostFunction46(t,u,x1,c)

% optimize input
N = (size(u,2)/2)/20;
C = [repelem(c(3:2:11),1,N) fliplr(repelem(c(3:2:11),1,N)) repelem(c(3:2:11),1,N) fliplr(repelem(c(3:2:11),1,N))];
N = (size(u,2)/2)/40;
C = [C repelem(c(2:11),1,N) fliplr(repelem(c(2:11),1,N)) repelem(c(2:11),1,N) fliplr(repelem(c(2:11),1,N))];
u = u.*C;

% figure()
% plot(C)
% offset with a delay time
t_fitted = t+c(1);
N = (size(u,2))/8;

%make fit without time adjustment
x1_fitted = [cumsum(u(1:N)) cumsum(u(N+1:N*2)) cumsum(u(1+N*2:N*3))...
    cumsum(u(N*3+1:4*N)) cumsum(u(N*4+1:5*N)) cumsum(u(N*5+1:6*N))...
    cumsum(u(N*6+1:7*N)) cumsum(u(N*7+1:8*N))]*0.01/0.0154;

% interpolate between new timeoffset to put in cost function
x1_fitted = interp1(t_fitted,x1_fitted,t);
% the beginnin contains NaN. replace those with zeros
x1_fitted(isnan(x1_fitted))=0;

% any(isnan(ref_fitted))

V = norm(x1_fitted-x1,2);
end





