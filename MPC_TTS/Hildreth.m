function [u,lambda] = Hildreth(E,invE,F,M,gamma,umin,umax,maxiter,ErrorBound)
%% 1 input (u) for now
% invE = inv(E);
u_unc = -invE*F; % calculate unconstrained solution

% construct matrices H and K used for calculation of w
H = M*invE*M';
K = gamma + M*invE*F;

ni = size(H,1);
disp(ni);

w = zeros(ni,1);
lambda_v = zeros(ni,maxiter);
lambda_mplus1 = zeros(ni,1);

for m = 1:maxiter
    for i = 1:ni
            term1 = 0;
            term2 = 0;
            if i > 1
                for j =1:i-1
                    term1 = term1+H(i,j)*lambda_mplus1(j);
                end
            end
            if i < ni
                for j = i+1:ni
                    term2 = term2+H(i,j)*lambda_v(j,m);
                end
            end
            w(i) = - 1/H(i,i)*(K(i)+term1+term2);
            lambda_mplus1(i,:) = max(0,w(i));
    end
    lambda_v(:,m+1) = lambda_mplus1;
    if abs(lambda_v(:,m+1)-lambda_v(:,m)) < ErrorBound*ones(ni,1)
        u = -invE*(F+M'*lambda_v(:,m+1));
        lambda = lambda_v(:,m+1);
        return;
    end
end
    disp('Warning: max iterations, convergence might not be achieved.')
    for i = 1:length(u_unc)
        if u_unc(i) < umin
            u(i) = umin;
        elseif u_unc(i) > umax
            u(i) = umax;
        else
            u(i)=u_unc(i);
        end
    end
end

