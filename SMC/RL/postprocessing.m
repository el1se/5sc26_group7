close all

% Tank heights referencs vs simulation
figure(1)
plot(out.References(:,1)*1000)
hold on
plot(out.References(:,2)*1000)
plot(out.References(:,3)*1000)
plot(out.Heights(:,1))
plot(out.Heights(:,2))
plot(out.Heights(:,3))
legend('T1 ref','T2 ref','T3 ref',...
       'T1 sim','T2 sim','T3 sim')
grid