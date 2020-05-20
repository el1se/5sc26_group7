clearvars;
close all;
clc;
load('ToRemoteLabs4_results.mat');
%% load input and reference 
%input and reference
f = 1/8; % frequency
mag = 0.02:0.02:0.1; % nr of different amplitudes
periods = 6; % nr of periods each magnitude occurs
factor = 1.5; % each period appears 1.5 times
Ts= 0.01; % sampling frequency.
UnitBlock = [zeros(1,1/(f*Ts)) ones(1,1/(f*2*Ts))];

% making input signal
n = periods*size(mag,2);
input1 = zeros(1,n*size(UnitBlock,2)+1);
for i = 1:size(mag,2)
    MagUnit = repmat(UnitBlock,1,periods)*mag(i);
    N = size(MagUnit,2);
    input1(1,N*(i-1)+2:N*i+1) = MagUnit;
end

input4 = [input1 zeros(1,150/Ts) fliplr(input1(2:end)) zeros(1,150/Ts)...
    (input1(2:end)) zeros(1,150/Ts) fliplr(input1(2:end)) zeros(1,150/Ts)...
    (input1(2:end)) zeros(1,150/Ts) fliplr(input1(2:end))];
Valves_vector  = boolean([zeros(1,size(input1,2)) ones(1,150/Ts)...
    zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1)...
    ones(1,150/Ts) zeros(1,size(input1,2)-1) ones(1,150/Ts)...
    zeros(1,size(input1,2)-1) ones(1,150/Ts) zeros(1,size(input1,2)-1)]);

t4 = 0:Ts:(size(input4,2)-1)*Ts;
reference = cumsum(input4*0.01/0.0154); % factor 100 because x1 in cm
for i = 2:length(Valves_vector)
   if(~Valves_vector(i-1) && (Valves_vector(i)))
       reference(i:end)=reference(i:end)-reference(i);
   end
end




%% plot refernce and input and result
%result
t = seconds(mdfData.Time);
x1 = double(mdfData.x1_mm);
figure()
plot(t,x1)
hold on
plot(t4, reference)

%% get first section of interest
N         = size(input1,2);
sec_t     = 0:Ts:(length(input1)-1)*Ts;
sec_u     = input1;
sec_x1    = double(mdfData.x1_mm(1:length(input1)))';
sec_ref   = reference(1:length(input1));

% calculate estimated delay
delayed_refsec = zeros(1,N);
t_delay = sqrt(2*0.001.*(620-sec_ref)/9.81);


figure()
plot(t_delay)
ylabel('max delay time [s]')
text(1e4,0.28,'delay is likely to be shorter due to inital velocity>0')

figure()
% plot(sec_t,sec_x1)
% hold on
plot(sec_t,sec_ref-sec_x1)

% close all

%% calculation of optimal values

c_InitGues = [0 1 1 1 1 1];
fun = @(t0) CostFunction(sec_t,sec_u,sec_x1,t0);
optimalSol = fmincon(fun,c_InitGues,[],[]);

t_fitted = sec_t+optimalSol(1);
N = (size(sec_u,2)-1)/5;
sec_u(2:N+1) = sec_u(2:N+1)*optimalSol(2);
sec_u(N*1+2:N*2+1) = sec_u(N*1+2:N*2+1)*optimalSol(3);
sec_u(N*2+2:N*3+1) = sec_u(N*2+2:N*3+1)*optimalSol(4);
sec_u(N*3+2:N*4+1) = sec_u(N*3+2:N*4+1)*optimalSol(5);
sec_u(N*4+2:N*5+1) = sec_u(N*4+2:N*5+1)*optimalSol(6);
ref_fitted = cumsum(sec_u*0.01/0.0154);
figure()
subplot(1,2,1)
plot(t_fitted,ref_fitted)
hold on
plot(sec_t,sec_x1)

%% validation

secval_u = fliplr(input1(2:end));
secval_t = Ts:Ts:size(secval_u,2)*Ts;
secval_x1 = double(mdfData.x1_mm(36001+1+15000:72001+15000))';

valt_fitted = secval_t+optimalSol(1);
N = (size(secval_u,2))/5;
secval_u(1:N) = secval_u(1:N)*optimalSol(6);
secval_u(N*1+1:N*2) = secval_u(N*1+1:N*2)*optimalSol(5);
secval_u(N*2+1:N*3) = secval_u(N*2+1:N*3)*optimalSol(4);
secval_u(N*3+1:N*4) = secval_u(N*3+1:N*4)*optimalSol(3);
secval_u(N*4+1:N*5) = secval_u(N*4+1:N*5)*optimalSol(2);
refval_fitted = cumsum(secval_u*0.01/0.0154);
% figure()
subplot(1,2,2)
plot(valt_fitted,refval_fitted)
hold on
plot(secval_t,secval_x1)
%% cost fun
function [V] = CostFunction(sec_t,sec_u,sec_x1,c)
N = (size(sec_u,2)-1)/5;
% optimize input
sec_u(2:N+1) = sec_u(2:N+1)*c(2);
sec_u(N*1+2:N*2+1) = sec_u(N*1+2:N*2+1)*c(3);
sec_u(N*2+2:N*3+1) = sec_u(N*2+2:N*3+1)*c(4);
sec_u(N*3+2:N*4+1) = sec_u(N*3+2:N*4+1)*c(5);
sec_u(N*4+2:N*5+1) = sec_u(N*4+2:N*5+1)*c(6);


t_fitted = sec_t+c(1);
ref_fitted = cumsum(sec_u*0.01/0.0154);
ref_fitted = interp1(t_fitted,ref_fitted,sec_t);
ref_fitted(isnan(ref_fitted))=0;
any(isnan(ref_fitted))
V = norm(ref_fitted-sec_x1,2)
end















