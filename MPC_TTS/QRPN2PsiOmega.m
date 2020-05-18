function [ Psi, Omega ] = QRPN2PsiOmega( Q,R,P,N )
% Omega
if (N==1)
    Omega=P;
else
    Omega=Q;
    for j=2:(N-1)
        Omega=blkdiag(Omega,Q);
    end
    Omega=blkdiag(Omega,P);
%     omega=blkdiag(omega,Q);
enda

% Psi
Psi = R;
if N>1
    for j = 2:N
        Psi = blkdiag(Psi,R);
    end
end
end

