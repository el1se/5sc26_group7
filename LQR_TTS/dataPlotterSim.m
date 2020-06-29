dat = out.simout.Data;
X = dat(:,1:5);
U = dat(:,6:9);
error = dat(:,10:12)*1000-X(:,1:3);
t = out.tout;
L = length(t);
LW = 1.3;
I = 11001;
%%
figure
subplot(2,2,1)
plot(t(1:I),X(1:I,1),'Linewidth',LW); hold on;
plot(t(1:I),X(1:I,2),'Linewidth',LW);
plot(t(1:I),X(1:I,3),'Linewidth',LW);
grid on;
xlim([0 125])
xlabel('Time [s]');
ylabel('Tank heights [mm]')
subplot(2,2,3)
plot(t(1:I),error(1:I,1),'Linewidth',LW); hold on;
plot(t(1:I),error(1:I,2),'Linewidth',LW);
plot(t(1:I),error(1:I,3),'Linewidth',LW);
xlim([0 125])
grid on;
legend('Tank L','Tank R','Tank M');
xlabel('Time [s]');
ylabel('Error[mm]');
subplot(2,2,2)
plot(t(1:I),X(1:I,4),'Linewidth',LW); hold on;
plot(t(1:I),X(1:I,5),'Linewidth',LW);
xlim([0 125])
grid on;
legend('Valve LM','Valve RM');
xlabel('Time [s]');
ylabel('Valve opening [%]')
subplot(2,2,4)
plot(t(1:I),U(1:I,1),'Linewidth',LW); hold on;
plot(t(1:I),U(1:I,2),'Linewidth',LW);grid on;
xlim([0 125])
legend('Pump L','Pump R')
xlabel('Time [s]');
ylabel('Pump input [L/s]');
%%
1/L*norm(error(:,1),2)
1/L*norm(error(:,2),2)
1/L*norm(error(:,3),2)
error(end,:)