mynorm = @(x)norm(x(:));
ndisp = 5; disprate = 10; qdisp = 1;
E1 = []; E2 = [];
clf;
for i=1:niter
    Pi = ProjC1(Pi);
    E2(i) = mynorm( sum(Pi,1)-permute(P,[2 1 3]) ) / mynorm(P);
    %
    Pi = ProjC2(Pi,P);
    p = prod(sum(Pi,2), 3).^(1/K); % average
    pp = repmat(p, [1 1 K]);
    E1(i) = mynorm( sum(Pi,2) - pp ) / mynorm(pp);
    if mod(i,disprate)==1 && qdisp<=ndisp
        subplot(ndisp,1,qdisp);
        bar(x, p, 'k');
        axis tight;
        qdisp = qdisp+1;
    end
end
