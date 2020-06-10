function [A,B]  = PID_fcn(Y,val)
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

tubeSpeed = val;


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



%%-------------------------------------------------------------------------
% functions
function Afcn = freeFallLinearizationA(x)
tol = 0;
if ((x(2)-x(3)) <= tol) && ((x(1)-x(3)) > tol)
    Afcn = [   
        -tubeSpeed*S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank    0    tubeSpeed*S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank     -tubeSpeed*S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank   0;
        zeros(1,5);
        tubeSpeed*S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank     0    -tubeSpeed*S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank    tubeSpeed*S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank    0;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
elseif ((x(1)-x(3)) <= tol) && ((x(2)-x(3)) > tol)
    Afcn = [   
        zeros(1,5);
        0                       -tubeSpeed*S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank    tubeSpeed*S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank     0    -tubeSpeed*S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       tubeSpeed*S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank     -tubeSpeed*S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank    0    tubeSpeed*S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
elseif ((x(2)-x(3)) <= tol) && ((x(1)-x(3)) <= tol)
     Afcn = [   
        zeros(3,5);
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
else
    Afcn = [   
        -S*tubeSpeed*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank     0                                                   tubeSpeed*S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank                                                      -tubeSpeed*S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank   0;
        0                                                   -S*tubeSpeed*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank     tubeSpeed*S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank                                                      0                                       -tubeSpeed*S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        S*tubeSpeed*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank      S*tubeSpeed*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank      -tubeSpeed*S*x(4)/Dvalve*g/sqrt(2*g*abs(x(1)-x(3)))/Atank-tubeSpeed*S*x(5)/Dvalve*g/sqrt(2*g*abs(x(2)-x(3)))/Atank      tubeSpeed*S/Dvalve*sqrt(2*g*abs(x(1)-x(3)))/Atank    tubeSpeed*S/Dvalve*sqrt(2*g*abs(x(2)-x(3)))/Atank;
        0                       0                       0                                   -v4                                     0;
        0                       0                       0                                   0                                       -v5
        ];
end
end
end


