%% reset add folder
clc; clear all; close all
addpath('Machine_Learning_Basic_Scripts');
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S/pi);     % [m]
%% training data
xTraining = [20.8214436313299,40.6830579434519,60.8568621315079,80.8895468214911,100]';
y(:,1) = [0.41 6.5 14 15 15.2]'; % ydata 1
y(:,2) = [0.14 5.8 13 14.5 15]'; % ydata 2

n =250;                                     % number of test points
N = length(xTraining);                            % number of training points
xPrior = linspace(0,100,n)';                % test points

h = @(x) exp(x*12/100-6)./(exp(x*12/100-6)+1);
mh = size(h(1),2);
%% opt vars
x01 = linspace(1,1000,10);
x02 = logspace(-4,4,10);
x03 = linspace(1,1000,10);
ub = [10000,1,10000];         % lower and upper bounds for hyper parameters
lb= [0.05^2 1e-14 0.05^2];
options = optimoptions('fmincon','Display','off',...
    'Algorithm','interior-point',...          % interior point does not work correctly with specifyobjectivegradient on
    'SpecifyObjectiveGradient',false,...
    'CheckGradients',false,...
    'StepTolerance',1e-10);
[X01,X02,X03] = ndgrid(x01,x02,x03);
%% loop over ydatasets
for i = 1:2
    H = h(xTraining)';
    Hs = h(xPrior)';
    betaBar(:,i) = inv(H*H')*H*y(:,i);
    for ii = 1:(size(X01,1))
        for jj = 1:(size(X02,2))
            for iii = 1:size(X03,3)
                fval(ii,jj,iii) = marLikelihood4hyp(xTraining,y(:,i),h,[X01(ii,jj,iii); X02(ii,jj,iii); X03(ii,jj,iii)],0);
            end
        end
    end
    
    mini = min(min(min(fval)));
    [I]=find(fval==mini);
    xres0 = [X01(I); X02(I);X03(I)];
    [xres,~] = fmincon(@(x) marLikelihood4hyp(xTraining,y(:,i),h,x,0),xres0,[],[],[],[],lb,ub,[],options);
    k = GPSEKernel(xTraining,xTraining,xres(1));
    k_s = xres(3)*GPSEKernel(xTraining,xPrior,xres(1));
    %L and Lk
    Ky = xres(3)*k+xres(2)*eye(N);
    L = chol(Ky,'lower');
    Lk = L \ k_s;
    
    R = Hs-H*inv(Ky)*k_s;
        
    % kernel of prediction
    k_ss = xres(3)*GPSEKernel(xPrior',xPrior',xres(1)) + R'*inv(H*inv(Ky)*H')*R;
    
    % mu etc
    mu(:,i) = ((Lk') * (L \ y(:,i))+R'*betaBar(:,i));
    var(:,i) = (diag(k_ss)' - sum(Lk.^2,1))';
    
end
mu(mu<0) = 0;
figure
subplot(1,2,1)
inBetween = [(mu(:,1)+3*sqrt(var(:,1)))' fliplr((mu(:,1)-3*sqrt(var(:,1)))')];
x2 = [xPrior', fliplr(xPrior')];
fill(x2,inBetween, [7 7 7]/8);hold on;
plot(xTraining,y(:,1),'s','MarkerSize',12);
plot(xPrior,mu(:,1),'LineWidth',1.3);
plot(xPrior,h(xPrior)*betaBar(:,1),'LineWidth',1.3);
xlabel('Valve opening percentage [%]');
ylabel('c');
legend('$\mu \pm 3\sigma$ for GP+LS','Data valve 1-3','Fit valve 1-3 GP+LS','Fit valve 1-3 LS only','Interpreter','Latex');
subplot(1,2,2)
inBetween = [(mu(:,2)+3*sqrt(var(:,2)))' fliplr((mu(:,2)-3*sqrt(var(:,2)))')];
x2 = [xPrior', fliplr(xPrior')];
fill(x2,inBetween, [7 7 7]/8);hold on;
plot(xTraining,y(:,2),'d','MarkerSize',12); 
plot(xPrior,mu(:,2),'LineWidth',1.3);
plot(xPrior,h(xPrior)*betaBar(:,2),'LineWidth',1.3);
xlabel('Valve opening percentage [%]');
ylabel('c');
legend('$\mu \pm 3\sigma$ for GP+LS','Data valve 2-3','Fit valve 2-3 GP+LS','Fit valve 1-3 LS only','Interpreter','Latex');

