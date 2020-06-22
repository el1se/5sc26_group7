close all

% Tank heights referencs vs simulation
figure(1)
subplot(2,1,1)
plot(out.tout,out.References(:,1)*1000)
hold on
plot(out.tout,out.References(:,2)*1000)
%plot(out.tout,out.References(:,3)*1000)
plot(out.tout,out.Heights(:,1))
plot(out.tout,out.Heights(:,2))
%plot(out.tout,out.Heights(:,3))
legend('T1 ref','T2 ref',...%'T3 ref',...
       'T1 sim','T2 sim')%,'T3 sim')
grid
xlabel('Time [s]')
ylabel('Height [mm]')

%subplot(3,1,2)
%plot(out.tout,out.error1)
%hold on
%plot(out.tout,out.error2)
%legend('error1','error2')
%grid
%xlabel('Time [s]')
%ylabel('Error [mm]')

subplot(2,1,2)
plot(out.tout,out.pump1)
hold on
plot(out.tout,out.pump2)
legend('pump1','pump2')
grid
xlabel('Time [s]')
ylabel('Input [L/s]')
