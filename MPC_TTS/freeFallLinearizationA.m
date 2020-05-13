function A = freeFallLinearizationA(x)
global Atank v4 v5 g Dvalve S
A = [   
    -S*x(4)/Dvalve*(2^(1/2)*g)/(2*(g*(x(1) - x(3)))^(1/2))/Atank   0  S*x(4)/Dvalve*(2^(1/2)*g)/(2*(g*(x(1) - x(3)))^(1/2))/Atank -S/Dvalve*sqrt(2*g*(x(1)-x(3)))/Atank 0;
    0 -S*x(5)/Dvalve*(2^(1/2)*g)/(2*(g*(x(2) - x(3)))^(1/2))/Atank S*x(5)/Dvalve*(2^(1/2)*g)/(2*(g*(x(2) - x(3)))^(1/2))/Atank 0 -S/Dvalve*sqrt(2*g*(x(2)-x(3)))/Atank;
    S*x(4)/Dvalve*(2^(1/2)*g)/(2*(g*(x(1) - x(3)))^(1/2))/Atank S*x(5)/Dvalve*(2^(1/2)*g)/(2*(g*(x(2) - x(3)))^(1/2))/Atank -S*x(4)/Dvalve*(2^(1/2)*g)/(2*(g*(x(1) - x(3)))^(1/2))/Atank-S*x(5)/Dvalve*(2^(1/2)*g)/(2*(g*(x(2) - x(3)))^(1/2))/Atank S/Dvalve*sqrt(2*g*(x(1)-x(3)))/Atank S/Dvalve*sqrt(2*g*(x(2)-x(3)))/Atank;
    0                       0                       0                                   -v4                                     0;
    0                       0                       0                                   0                                       -v5
    ];

end

