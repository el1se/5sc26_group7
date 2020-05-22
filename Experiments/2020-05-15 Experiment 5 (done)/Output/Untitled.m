x = 0.01:0.01:0.1;
y = c46_x1_optimum(2:11);

eg = 1; % sum-squared error goal
sc = 4;    % spread constant
net = newrb(x,y,1,4);

X = 0.01:.0001:0.1;
X= x;
Y = net(X);

plot(x,y,'+');

title('Training Vectors');
xlabel('Input Vector P');
ylabel('Target Vector T');
hold on;
plot(X,Y);
hold off;
legend({'Target','Output'})

% variance = norm(y-Y);
% variance0 = 10;
% variance = zeros(length(0.001:0.001:1),length(0.01:0.01:4));
% egvec = 0.01:0.01:1;
% scvec = 0.1:0.1:4;
% variance = zeros(length(egvec),length(scvec));
% 
% for i = 1:length(egvec)
%     for ii = 1:length(scvec)
%         eg = egvec(i); % sum-squared error goal
%         sc = scvec(ii);    % spread constant
%         net = newrb(x,y,eg,sc);
%         Y = net(x); 
%         
%         variance(i,ii) = norm(y-Y);
%         
%     end
% end