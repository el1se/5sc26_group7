syms S x1 x2 x3 x4 x5 g c1 c2 cc A cd real
Q13 = S*c1*(exp(x4*cc))/(exp(x4*cc)+1)*sqrt(2*g*(x1-x3));
Q23 = S*c2*(exp(x5*cc))/(exp(x5*cc)+1)*sqrt(2*g*(x2-x3));

f = 1/A*[-Q13 -Q23 Q13+Q23 -cd*x4 -cd*x5]';
Lfh3 = [0 0 1 0 0]*f;
dLfh3dx = [diff(Lfh3,x1) diff(Lfh3,x2) diff(Lfh3,x3) diff(Lfh3,x4) diff(Lfh3,x5)];
Lf2h3 = dLfh3dx*f;
disp(simplify(Lf2h3));

alpha = [-1/A*Q13; -1/A*Q23; Lf2h3];
