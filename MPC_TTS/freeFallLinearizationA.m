function A = freeFallLinearizationA(x)
global Atank v4 v5 g Dvalve S
tol = 0;
if ((x(2)-x(3)) < tol) && ((x(1)-x(3)) >= tol)
    A = [   
        -S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank    0    S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank     -sign(x(1)-x(3))*S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank   0;
        zeros(1,5);
        S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank     0    -S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank    sign(x(1)-x(3))*S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank    0;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
elseif ((x(1)-x(3)) < tol) && ((x(2)-x(3)) >= tol)
    A = [   
        zeros(1,5);
        0                       -S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank    S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank     0    -sign(x(2)-x(3))*S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank     -S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank    0    sign(x(2)-x(3))*S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
elseif ((x(2)-x(3)) < tol) && ((x(1)-x(3)) < tol)
     A = [   
        zeros(3,5);
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
else
    A = [   
        -S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank     0                                                   sign(x(1)-x(3))*S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank                                                      -S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank   0;
        0                                                   -S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank     sign(x(2)-x(3))*S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank                                                      0                                       -S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank      S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank      -sign(x(1)-x(3))*S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank-sign(x(2)-x(3))*S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank      S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank    S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
end
end

