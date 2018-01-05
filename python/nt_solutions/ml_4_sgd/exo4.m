ElistG = [];
tau = .002/n;
nsamples = 10;
err_rate = 50;
for is=1:nsamples
    w = zeros(p,1);
    G = zeros(p,n); % keep track of gradients
    g = zeros(p,1);
    for it=1:niter
        if mod(it,err_rate)==1
            ElistG( 1+(it-1)/err_rate,is) = E(w,X,y);
        end
        i = 1+floor(rand*n); % draw uniformly
        g1 = nablaEi(w,i);
        % update grad
        g = g - G(:,i) + g1;
        G(:,i) = g1;
        %
        w = w - tau * g;
    end
end
clf;
hold on;
plot(1,Inf, 'b'); plot(1,Inf, 'r'); plot(1,Inf, 'g');
plot(1:err_rate:niter, log10(ElistS-min(Elist)), 'b');
plot(1:err_rate:niter, log10(ElistA-min(Elist)), 'r');
plot(1:err_rate:niter, log10(ElistG-min(Elist)), 'g');
axis tight; box on;
SetAR(1/2);
title('log(E(w_l) - min E)'); set(gca, 'FontSize', fs);
legend('SGD', 'SGA', 'SAG');
