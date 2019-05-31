niter = 1000;
f = zeros(n,1);
Err = [];
for it=1:niter
    g = mina(C-f,epsilon);
    f = minb(C-g,epsilon);
    % generate the coupling
    P = a .* exp((f+g-C)/epsilon) .* b;
    % check conservation of mass
    Err(it) = norm(sum(P,1)-b,1);    
end
clf;
plot(log10(Err), 'LineWidth', 2);
