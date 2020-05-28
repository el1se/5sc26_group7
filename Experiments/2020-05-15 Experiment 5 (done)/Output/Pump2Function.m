function [ActualOutput] = Pump2Function(TheoreticInput)
%{
it is a second order polynomial with vector input possible
u_actual = c(u_theory)*u_theory
%}

% coeficients
c = [-77.3280037134416,10.5485815604300,0.778096847837038];

% output
ActualOutput = (c(1)*TheoreticInput.*TheoreticInput + c(2)*TheoreticInput + c(3)).*TheoreticInput;

end


