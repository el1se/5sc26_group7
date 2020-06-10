X = out.simout.Data(:,1:5);
U = out.simout.Data(:,6:9);
t = out.tout;
LW = 1.3;
figure
subplot(2,2,1)
plot(t,X(:,1),'Linewidth',LW); hold on;
plot(t,X(:,2),'Linewidth',LW);
plot(t,X(:,3),'Linewidth',LW);
% plot(t,zeros(length(t),1)+305,'r--');
% plot(t,zeros(length(t),1)+295,'r--');
legend('Tank L','Tank R','Tank M');
xlabel('Time [s]');
ylabel('Tank heights [mm]')
subplot(2,2,2)
plot(t,X(:,4),'Linewidth',LW); hold on;
plot(t,X(:,5),'Linewidth',LW);
legend('Valve LM','Valve RM');
xlabel('Time [s]');
ylabel('Valve opening [%]')
subplot(2,2,3)
plot(t,U(:,1),'Linewidth',LW); hold on;
plot(t,U(:,2),'Linewidth',LW);
legend('Pump L','Pump R')
xlabel('Time [s]');
ylabel('Pump input [L/s]');
subplot(2,2,4)
plot(t,U(:,3),'Linewidth',LW); hold on;
plot(t,U(:,4),'Linewidth',LW);
legend('valve ML','Valve MR')
xlabel('Time [s]');
ylabel('Valve Input [%]');