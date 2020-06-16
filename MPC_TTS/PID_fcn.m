function [A,B]  = PID_fcn(Y,bB)
% outputs 4d vector u containing inputs for TTS  (format pump1 pump2 valve LM valve RM)
% outputs 3x3 matrix P, containing solution to DARE, can be used in
% next itration
% Y: 5d input state vector of TTS (format water level tank 1, tank 2, tank 3, connecting pipe diameter LM)
% uprev: previous input sent to plant (4D vector, same format as u), for augmented model necessary
% Pprev: previous solution to DARE, might be used if no stabilizing
% solution can be calculated
% Yprev: previous 5d state vector of TTS, same format as Y
%% ini
u = zeros(4,1);
%% constant
S = 5e-5;               % [m^2]
Dvalve = sqrt(4*S/pi);     % [m]
Atank = 0.0154;         % [m^2]
g = 9.81;               % [m/s^2]
L = 0.1;                % [m]
mu = 1.002;             % [Ns/m^2]
rho = 1000;             % [kg/m^3]
% system constants (ID or whatever)
v4 = 0.5;                 % [s]?
v5 = 0.5;                 % [s]?

% Constraints 
Hmax = 0.61;            % [m]
Qmax = 0.1;             % [l/s]



%% model
B = [0.001/Atank 0 0 0;
        0 0.001/Atank 0 0;
        0 0 0 0 ;
        0 0 0.01*v4*Dvalve 0;
        0 0 0 0.01*v5*Dvalve];
C = eye(5);


m = size(B,2);
n = size(B,1);
q = size(C,1);

Dm = zeros(q,m);

A = freeFallLinearizationA(Y);
%% ------------------------------------------------------------------------
% functions
function Afcn = freeFallLinearizationA(x)
tol = 0;
f= @(xi,beta) beta * (exp(xi/Dvalve*12-6))./(exp(xi/Dvalve*12-6)+1);
dfdx = @(xi,beta) Dvalve*12*exp(xi/Dvalve*12+6)./(exp(xi/Dvalve*12)+exp(6))^2;
if ((x(2)-x(3)) <= tol) && ((x(1)-x(3)) > tol)
    Afcn = [   
        -S*f(x(4),bB(1))*g/sqrt(2*g*abs(x(1)-x(3)))/Atank           0      S*f(x(4),bB(1))*g/sqrt(2*g*abs(x(1)-x(3)))/Atank         -dfdx(x(4),bB(1))*S*sqrt(2*g*abs(x(1)-x(3)))/Atank   0;
        zeros(1,5);
        S*f(x(4),bB(1))*g/sqrt(2*g*abs(x(1)-x(3)))/Atank            0      -S*f(x(4),bB(1))*g/sqrt(2*g*abs(x(1)-x(3)))/Atank        dfdx(x(4),bB(1))*S*sqrt(2*g*abs(x(1)-x(3)))/Atank    0;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
elseif ((x(1)-x(3)) <= tol) && ((x(2)-x(3)) > tol)
    Afcn = [   
        zeros(1,5);
        0    -S*f(x(5),bB(2))*g/sqrt(2*g*abs(x(2)-x(3)))/Atank      S*f(x(5),bB(2))*g/sqrt(2*g*abs(x(2)-x(3)))/Atank   0             -dfdx(x(5),bB(2))*S*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0    S*f(x(5),bB(2))*g/sqrt(2*g*abs(x(2)-x(3)))/Atank       -S*f(x(5),bB(2))*g/sqrt(2*g*abs(x(2)-x(3)))/Atank  0             dfdx(x(5),bB(2))*S*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
elseif ((x(2)-x(3)) <= tol) && ((x(1)-x(3)) <= tol)
     Afcn = [   
        zeros(3,5);
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
else % x1 > x3 & x2 > x3
    Afcn = [   
        -S*f(x(4),bB(1))*g/sqrt(2*g*abs(x(1)-x(3)))/Atank           0      S*f(x(4),bB(1))*g/sqrt(2*g*abs(x(1)-x(3)))/Atank         -dfdx(x(4),bB(1))*S*sqrt(2*g*abs(x(1)-x(3)))/Atank   0;
        0    -S*f(x(5),bB(2))*g/sqrt(2*g*abs(x(2)-x(3)))/Atank      S*f(x(5),bB(2))*g/sqrt(2*g*abs(x(2)-x(3)))/Atank   0             -dfdx(x(5),bB(2))*S*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        S*f(x(4),bB(1))*g/sqrt(2*g*abs(x(1)-x(3)))/Atank            S*f(x(5),bB(2))*g/sqrt(2*g*abs(x(2)-x(3)))/Atank      -S*f(x(4),bB(1))*g/sqrt(2*g*abs(x(1)-x(3)))/Atank-S*f(x(5),bB(2))*g/sqrt(2*g*abs(x(2)-x(3)))/Atank      dfdx(x(4),bB(1))*S*sqrt(2*g*abs(x(1)-x(3)))/Atank    dfdx(x(5),bB(2))*S*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
end
end
end


