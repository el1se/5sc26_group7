function A = laminarPoisseuilleFlowLinearizatonA(x)
global Atank v4 v5 g rho mu L
frac = 128*mu*L*Atank;
A = [   
    -pi*x(4)^4*rho*g/frac   0                       pi*x(4)^4*rho*g/frac                -pi*4*x(4)^3*rho*g*(x(1)-x(3))/frac     0;
    0                       -pi*x(5)^4*rho*g/frac   pi*x(5)^4*rho*g/frac                0                                       pi*4*x(5)^3*rho*g*(x(3)-x(2))/frac;
    pi*x(4)^4*rho*g/frac    pi*x(5)^4*rho*g/frac    -pi*(x(5)^4+x(4)^4)*rho*g/frac      4*pi*x(4)^3*rho*(x(1)-x(3))/frac        -4*pi*x(5)^3*rho*(x(3)-x(2))/frac;
    0                       0                       0                                   -v4                                     0;
    0                       0                       0                                   0                                       -v5
    ];
end

