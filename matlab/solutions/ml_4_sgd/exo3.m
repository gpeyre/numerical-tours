ElistA = [];
tau0 = .05;
ell0 = 100;
nsamples = 10;
err_rate = 50;
for is=1:nsamples
    w = zeros(p,1); w1 = w;
    for it=1:niter
        if mod(it,err_rate)==1
            ElistA( 1+(it-1)/err_rate,is) = E(w,X,y);
        end
        tau = tau0 / (1+sqrt(it/ell0));
        i = 1+floor(rand*n); % draw uniformly
        w1 = w1 - tau * nablaEi(w1,i);
        w = 1/it*w1 + (1-1/it)*w;
    end
end
clf;
hold on;
plot(1,Inf, 'b'); plot(1,Inf, 'r');
plot(1:err_rate:niter, log10(ElistS-min(Elist)), 'b');
plot(1:err_rate:niter, log10(ElistA-min(Elist)), 'r');
axis tight; box on;
SetAR(1/2);
title('log(E(w_l) - min E)'); set(gca, 'FontSize', fs);
legend('SGD', 'SGA');
