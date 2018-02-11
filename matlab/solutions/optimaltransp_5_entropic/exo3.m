b = ones(N,1);
niter = 2000;
Err_p = []; Err_q = []; 
for i=1:niter
    a = p ./ (xi*b);
    Err_q(end+1) = norm( b .* (xi*a) - q )/norm(q);
    b = q ./ (xi'*a);    
    Err_p(end+1) = norm( a .* (xi'*b) - p )/norm(p);
end
% Display the violation of constraint error in log-plot. 
clf;
subplot(2,1,1);
plot(log10(Err_p)); axis tight; title('log|\pi 1 - p|');
subplot(2,1,2);
plot(log10(Err_q)); axis tight; title('log|\pi^T 1 - q|');
