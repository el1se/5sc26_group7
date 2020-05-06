function [D,M,E,c] = DMEcIntegral(deltaumin,deltaumax,umin,umax,ymin,ymax,ukmin1,N)
m = size(umin,1);
q = size(ymin,1);

bi = [-deltaumin; deltaumax; -umin+ukmin1;umax-ukmin1;-ymin;ymax];
bN = [zeros(4*m,1);-ymin; ymax];
Lb = length(bi);

c = zeros((4*m+2*q)*(N+1),1);
for i = 1:N
    c((i-1)*Lb+1:i*Lb) = bi;
end
c(N*Lb+1:end,:) = bN;


M0 = [zeros(4*m,q); -eye(q); eye(q)];

D = zeros((4*m+2*q)*(N+1),q);
D(1:size(M0,1),:)  = M0;

M = zeros((4*m+2*q)*(N+1),q*N);
for i = 1:N
    M(i*Lb+1:(i+1)*Lb,(i-1)*q+1:i*q) = M0;
end

Ei = [-eye(m); eye(m);-eye(m); eye(m);zeros(2*q,m)];
E0ff = [zeros(2*m,m);-eye(m); eye(m);zeros(2*q,m)];

E = zeros((4*m+2*q)*(N+1),m*N);
for i = 1:N
    E((i-1)*Lb+1:(i)*Lb,(i-1)*m+1:i*m) = Ei;
    if i > 1
        for j = 1:i-1
            E((i-1)*Lb+1:(i)*Lb,(j-1)*m+1:j*m) = E0ff;
        end
    end
end

end

