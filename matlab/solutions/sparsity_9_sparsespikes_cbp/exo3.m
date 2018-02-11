% list of regul param
Lmax = max(PhiS(y));
Llist = linspace(Lmax,1e-6,30);
% warmup stage
u = zeros(N,2); % 
niter = 3000;
lambda = Llist(1);
for j=1:niter
	u = ProxJ( u - tau * gradF(u), lambda*tau );
end
% outer loop
niter = 500;
A = []; B = [];
for i=1:length(Llist)
    % progressbar(i, length(Llist));
    lambda = Llist(i);
    % inner loop
    for j=1:niter
        u = ProxJ( u - tau * gradF(u), lambda*tau );
    end
    a = u(:,1);  b = u(:,2);
    A(:,end+1) = a; B(:,end+1) = b;
end
% display
Ic= setdiff(1:N,I);
clf;
subplot(2,1,1); hold on;
plot(Llist, A(I,:)', 'LineWidth', lw);
plot(Llist, A(Ic,:)', 'k-', 'LineWidth', lw); axis tight; box on;
title('a');
subplot(2,1,2); hold on;
plot(Llist, B(I,:)', 'LineWidth', lw);
plot(Llist, B(Ic,:)', 'k-', 'LineWidth', lw); axis tight; box on;
axis([0 max(Llist) -1 1]);
title('b');
