function [c] = GPGibValuePlz(x)
% input x 2d vector with openings
% output c 2d vector with c values
x = x(:);
c = zeros(2,1);
xTraining = [20.8214436313299,40.6830579434519,60.8568621315079,80.8895468214911,100]';
xres = [[3009.84105647172;0.999996056412786;118.502573365194] [3217.69210643946;0.885620951056251;112.696852154946]];
y(:,1) = [0.41 6.5 14 15 15.2]'; % ydata 1
y(:,2) = [0.14 5.8 13 14.5 15]'; % ydata 2
N = length(xTraining);
h = @(x) exp(x*12/100-6)./(exp(x*12/100-6)+1);
betaBar = [16.1461609250188,15.5165993538781];
for i = 1:2
    H = h(xTraining)';
    Hs = h(x(i))';
    k = GPSEKernelInline(xTraining,xTraining,xres(1,i));
    k_s = xres(3,i)*GPSEKernelInline(xTraining,x(i),xres(1,i));
    %L and Lk
    Ky = xres(3,i)*k+xres(2,i)*eye(N);
    L = chol(Ky,'lower');
    Lk = L \ k_s;
    
    R = Hs-H*inv(Ky)*k_s;
        
    % kernel of prediction
    k_ss = xres(3,i)*GPSEKernelInline(x(i)',x(i)',xres(1,i)) + R'*inv(H*inv(Ky)*H')*R;
    
    % mu etc
    c(i) = (Lk') * (L \ y(:,i))+R'*betaBar(i);
end
function val = GPSEKernelInline(a,b,l)
    %% calculation
    a = a(:);
    b = b(:);
    sqdist = (repmat(a,size(b'))-repmat(b',size(a))).^2;
    val = exp(-0.5/l*sqdist);
end


end

