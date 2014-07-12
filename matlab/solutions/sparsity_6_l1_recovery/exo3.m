N = 600;
P = N/2;
Phi = PhiRand(P,N);
klist = round(linspace( 1,P/7,20 ));
ntrials = 60;
proba = [];
for i=1:length(klist)
    proba(i,1:3) = 0;
    k = klist(i);
    for j=1:ntrials
        s = zeros(N,1); I = randperm(N); I = I(1:k);
        s(I) = sign(randn(k,1));
        proba(i,1) = proba(i,1) + (F(Phi,s)<1);
        proba(i,2) = proba(i,2) + (erc(Phi,I)<1);
        proba(i,3) = proba(i,3) + (werc(Phi,I)>0 & werc(Phi,I)<1);        
    end
end
clf;
h = plot(klist, proba/ntrials); set(h, 'LineWidth', 2);
xlabel('k');  legend('F<1', 'ERC<1', 'w-ERC<1');
title(sprintf('N=%d, P=%d', N, P));
axis tight;
