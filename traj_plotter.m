figure
subplot(1,3,1)
plot(out.simout.Time,out.simout.data(:,1),'LineWidth',1.3); hold on;
plot(out.simout.Time,out.simout.data(:,2),'LineWidth',1.3);
plot(out.simout.Time,out.simout.data(:,3),'LineWidth',1.3);
xlabel('Time [s]');
ylabel('Position reference [m]');
xlim([0 120])
subplot(1,3,2)
plot(out.simout.Time,out.simout.data(:,4),'LineWidth',1.3); hold on;
plot(out.simout.Time,out.simout.data(:,5),'LineWidth',1.3);
plot(out.simout.Time,out.simout.data(:,6),'LineWidth',1.3);
xlabel('Time [s]');
ylabel('Velocity reference reference [m/s]');
xlim([0 120])
subplot(1,3,3)
plot(out.simout.Time,out.simout.data(:,7),'LineWidth',1.3); hold on;
plot(out.simout.Time,out.simout.data(:,8),'LineWidth',1.3);
plot(out.simout.Time,out.simout.data(:,9),'LineWidth',1.3);
xlabel('Time [s]');
ylabel('Acceleration reference [m/s2]');
legend('Tank 1','Tank 2','Tank 3');
xlim([0 120])