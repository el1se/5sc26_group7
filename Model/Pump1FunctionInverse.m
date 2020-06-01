function [ActualOutput] = Pump1FunctionInverse(TheoreticInput)
% max = 0.0935
U=0:0.001:0.1;

% coeficients
c = [-106.218740175917,14.5693164637851,0.543210017066342];

% output
ActualOutputUncompensated = (c(1)*U.*U + c(2)*U + c(3)).*U;
figure()
plot(U,ActualOutputUncompensated)

% now compensate
ActualOutput=interp1(ActualOutputUncompensated,U,TheoreticInput);

end
