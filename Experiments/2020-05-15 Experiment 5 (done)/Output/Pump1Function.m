function [ActualOutput] = Pump1Function(TheoreticInput)
%{
so when you apply this there will be delay of a second for this pump 1,
however the fit results in a delay which is almost neglegible.

it is a second order polynomial with vector input possible
u_actual = c(u_theory)*u_theory
%}

% coeficients
c = [-106.218740175917,14.5693164637851,0.543210017066342];

% output
ActualOutput = (c(1)*TheoreticInput.*TheoreticInput + c(2)*TheoreticInput + c(3)).*TheoreticInput;

end