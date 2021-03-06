function [ActualOutput] = Pump2FunctionInverse(TheoreticInput)
% max input = 0.105
if(TheoreticInput>0.105)
    return
end

U=0:0.001:0.1;

% coeficients
c = [-77.3280037134416,10.5485815604300,0.778096847837038];

% output
ActualOutput2Uncompensated = (c(1)*U.*U + c(2)*U + c(3)).*U;
figure()
plot(U,ActualOutput2Uncompensated)

% now compensate
ActualOutput=interp1(ActualOutput2Uncompensated,U,TheoreticInput);

end
