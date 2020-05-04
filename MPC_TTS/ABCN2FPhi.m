function [F, Phi ] = ABCN2FPhi(A,B,C,N)
%ABN2PhiGamma Summary of this function goes here
%   Detailed explanation goes here

n = size(B,1);
m = size(B,2);
q = size(C,1);


F = [];
for i = 1:N
   F = [F;C*A^i] 
end

% phi
Phi = zeros(q*N,m*N);
for i = 1:N
    for j = 1:i
        Phi((i-1)*q+1:q*i,(j-1)*m+1:j*m) = C*A^(i-j)*B;
    end
end
end
