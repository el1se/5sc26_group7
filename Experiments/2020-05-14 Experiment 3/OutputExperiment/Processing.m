load('results.mat')
%% changing names
% names = mdfData.Properties.VariableNames;
% names{1,1} = 'x1_mm';
% names{1,2} = 'x3_mm';
% names{1,3} = 'x2_mm';
% names{1,4} = 'x1LeakPos';
% names{1,5} = 'x1x3ValvePos';
% names{1,6} = 'x2LeakPos';
% names{1,7} = 'x3x2ValvePos';
% names{1,8} = 'x3LeakPos';
% mdfData.Properties.VariableNames = names;
% % save('2080_results.mat','mdfData')

%% plotting part1
% reference
Ts= 0.01; % sampling time
f = 1/4; % frequency
mag = 0.01:0.01:0.1; % nr of different amplitudes
periods = 8; % nr of periods each magnitude occurs
t = 0:Ts:size(mag,2)*periods/f; % duration time

% creat blog signal with frequency
data1 = max(0,sign(sin(2*pi*f*t-0.5*pi))); 

% add a ramp to the input signal
data1 = data1 .*[0 repelem(mag,periods/(f*Ts))];
  
% total tank height
trapz(t,data1/0.0154)

% level in x2 was in experiments mean 

% plots
% state over time
figure(1)
set(gca,'LooseInset',get(gca,'TightInset'));

hold off
plot(seconds(mdfData_Part1.Time),double(mdfData_Part1.x1_mm),'linewidth',1)
hold on
plot(seconds(mdfData_Part1.Time),double(mdfData_Part1.x2_mm),'linewidth',1)
plot(t,cumsum(data1)*0.01/0.0154,'linewidth',1)
% plot(seconds(mdfData_Part1.Time),movmean(double(mdfData_Part1.x1_mm),5),'linewidth',1)
legend('x_1','x_2','reference','location','south')
grid on

% valve postions
figure(2)
hold off
plot(seconds(mdfData_Part1.Time),double(mdfData_Part1.x1LeakPos),'linewidth',1)
hold on
plot(seconds(mdfData_Part1.Time),double(mdfData_Part1.x2LeakPos),'linewidth',1)
plot(seconds(mdfData_Part1.Time),double(mdfData_Part1.x3LeakPos),'linewidth',1)
plot(seconds(mdfData_Part1.Time),double(mdfData_Part1.x1x3ValvePos),'linewidth',1)
plot(seconds(mdfData_Part1.Time),double(mdfData_Part1.x3x2ValvePos),'linewidth',1)
legend('x1LeakPos','x2LeakPos','x3LeakPos','x1x3ValvePos','x3x2ValvePos')
% % derived from x1
% timearray = seconds(mdfData_Part1.Time); x1 = movmean(double(mdfData_Part1.x1_mm),20);
% size(x1,1)
% x1dot=zeros(size(x1,1),1);
% x2dot=zeros(size(x1,1),1);
% 
% 
% for i=2:size(x1,1)-1
%     
%     x1dot(i) = (x1(i+1)-x1(i-1))/(2*Ts);
%     x2dot(i) = (x1(i+1)+x1(i-1)-2*x1(i))/(Ts*Ts);
% end
% 
% 
% figure(3)
% hold off
% plot(timearray,x2dot,'linewidth',1)



%% plotting part1
figure(4)
t2 = 0:Ts:max(t)*3;
data2 = [data1 data1(2:end)*0 fliplr(data1(2:end))];
data3 = cumsum(data2)*0.01/0.0154; data3(size(t,2)+1:end) = data3(size(t,2)+1:end)- data3(size(t,2));


set(gca,'LooseInset',get(gca,'TightInset'));

hold off
plot(seconds(mdfData_Part2.Time),double(mdfData_Part2.x1_mm),'linewidth',1)
hold on
plot(seconds(mdfData_Part2.Time),double(mdfData_Part2.x2_mm),'linewidth',1)
plot(t2,data3,'linewidth',1)
% plot(seconds(mdfData_Part1.Time),movmean(double(mdfData_Part1.x1_mm),5),'linewidth',1)
legend('x_1','x_2','reference','location','south')
grid on
figure(2)
hold off
plot(seconds(mdfData_Part2.Time),double(mdfData_Part2.x1LeakPos),'linewidth',1)
hold on
plot(seconds(mdfData_Part2.Time),double(mdfData_Part2.x2LeakPos),'linewidth',1)
plot(seconds(mdfData_Part2.Time),double(mdfData_Part2.x3LeakPos),'linewidth',1)
plot(seconds(mdfData_Part2.Time),double(mdfData_Part2.x1x3ValvePos),'linewidth',1)
plot(seconds(mdfData_Part2.Time),double(mdfData_Part2.x3x2ValvePos),'linewidth',1)
legend('x1LeakPos','x2LeakPos','x3LeakPos','x1x3ValvePos','x3x2ValvePos')