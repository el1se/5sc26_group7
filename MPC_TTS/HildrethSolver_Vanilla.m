function [U_k] = HildrethSolver_Vanilla(G,F,L,cplusWxk,MaxIter,ErrorBound)
% this is our vanilla Hildreth solver and this script has some differnt varianames, but is for the rest identical
% The equation (1/2)u'Gu+u'F, s.t. Lu <= cplusWxk is optimized

% calculate unconstraint solution first
U_k_unconstraint = -G\(F);

if(any(L*U_k_unconstraint>cplusWxk))
    % constraints violated
    H = L*(G\transpose(L));
    K = cplusWxk + L*(G\(F));
    % now ideally -H*lambdastar=K, but often not invertible
    n = size(H,1); % nr of lambdas to solve
    lambdas = zeros(n,1);
    m=1; % iteration counter
    
    while(1)
        % update for next interation
        lambdas_old = lambdas;
        
        % update each scalar value
        for i = 1:n
            lambdas(i,1) = max(0,-(K(i) + H(i,:)*lambdas - H(i,i)*lambdas(i,1))/H(i,i));
        end     
        % if maximum is reached without finding a feasible solution
        if(m>MaxIter)
           fprintf('No feasible solution found within nr. of iterations.\n');
           U_k = -G\(F+transpose(L)*lambdas);
           return
%            break
        end
        % check for convergence or timelimit
        if(transpose(lambdas-lambdas_old)*(lambdas-lambdas_old)<=ErrorBound)
            break
        end
        
        % update nr of iterations
        m = m+1;
    end
    
    % calculate optimal input
    U_k = -G\(F+transpose(L)*lambdas);
else
    % constraints not violated
    U_k = U_k_unconstraint;
%     lambdas = zeros(size(L,1),1);
end
end