% Number of paths computed in parallel.
q=50; 
% Initial conditions.
x = rand(q,1)-1/2;
niter = 1000;
E = [];
for i=1:niter-1
    u = rand(q,1)>.5;
    E(:,i) = F(x(:,end));
    tau = 1/(10+i);
    x(:,end+1) = x(:,end) - tau * ( u.*df{1}(x(:,end)) + (1-u).*df{2}(x(:,end)) );
end
%
clf;
plot(x', 'r');
axis tight;
set(gca, 'FontSize', 20);
