syms S x1 x2 x3 x4 x5 g c1 c2 cc A cd  x1r x2r x3r x3dot real

% k = sqrt(2*g/abs(x(1)-x(3)))*(x(1)-x(3));
% Q13 = S*c1*(exp(x4*cc))/(exp(x4*cc)+1)*sqrt(2*g*(x1-x3));
% Q23 = S*c2*(exp(x5*cc))/(exp(x5*cc)+1)*sqrt(2*g*(x2-x3));
Q13 = S*c1*(exp(x4*cc))/(exp(x4*cc)+1)*sqrt(2*g/abs(x1-x3))*(x1-x3);
Q23 = S*c2*(exp(x5*cc))/(exp(x5*cc)+1)*sqrt(2*g/abs(x2-x3))*(x2-x3);
% Q13 = (x4*S*sqrt(2*g/abs(x1-x3))*(x1-x3));
% Q23 = (x5*S*sqrt(2*g/abs(x2-x3))*(x2-x3));


f = 1/A*[-Q13 -Q23 Q13+Q23 -cd*x4 -cd*x5]';
Lfh3 = [0 0 1 0 0]*f;
dLfh3dx = [diff(Lfh3,x1) diff(Lfh3,x2) diff(Lfh3,x3) diff(Lfh3,x4) diff(Lfh3,x5)];
Lf2h3 = dLfh3dx*f;
% disp(simplify(Lf2h3));

alpha = [-1/A*Q13; -1/A*Q23; Lf2h3];

g1 = [1/A; 0 ;0;0;0];
g2 = [0;1/A;0;0;0];
g3 = [0;0;0;cd;0];
g4 = [0;0;0;0;cd];
Lg1Lfh3 = dLfh3dx*g1;
Lg2Lfh3 = dLfh3dx*g2;
Lg3Lfh3 = dLfh3dx*g3;
Lg4Lfh3 = dLfh3dx*g4;
beta = [1/A 0 0 0;
        0 1/A 0 0;
        Lg1Lfh3 Lg2Lfh3 Lg3Lfh3 Lg4Lfh3];
    
betaP = pinv(beta);

v = [-(x1-x1r); -(x2-x2r); -(x2-x3r)-2*(x3dot)];

u = simplify(betaP*(v-alpha));

%% 

% xr = [0.2 0.3 0.25];
% x3dot = 0.5;
% % xr = [ref(1) ref(2) ref(3)];
% A = 0.0154; S = 5e-5; g = 9.81;
% x1 = 0.1; x2 = 0.15; x3 = 0.12; x4 = 50; x5 = 70; g = 9.81; 
% c1 = 16.1462; c2 = 15.5166; cc = 12/100; cd = 100/99;
% x1r = 0.4; x2r = 0.3; x3r = 0.2;
% syms S x1 x2 x3 x4 x5 g c1 c2 cc A cd  x1r x2r x3r x3dot real
u = subs(u,[S x1 x2 x3 x4 x5 g c1 c2 cc A cd x1r x2r x3r x3dot],[5e-5 0.1 0.2 0.15 90 70 9.81 16 15 12/100 0.0154 100/99 0.4 0.3 0.35 0.5]);
