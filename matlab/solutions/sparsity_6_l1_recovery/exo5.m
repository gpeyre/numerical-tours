N = 600;
P = N/2;
Phi = PhiRand(P,N);
klist = round(linspace( 1,P/15,10 ));
ntrials = 2000;
rip_val = [];
for i=1:length(klist)
    rip_val(i,1:2) = 0;
    k = klist(i);
    for j=1:ntrials
        I = randperm(N); I = I(1:k);
        [a,b] = ric(Phi(:,I));
        rip_val(i,1:2) = max(rip_val(i,1:2), [a b]);
    end
end
clf;
hold on;
h = plot(klist, rip_val);  set(h, 'LineWidth', 2); 
h = plot(klist, klist*0+sqrt(2)-1, 'r:');  set(h, 'LineWidth', 2); 
legend('\delta^2_k', '\delta^2_k', '0.41'); xlabel('k');
title(sprintf('N=%d, P=%d', N, P));
axis tight;
