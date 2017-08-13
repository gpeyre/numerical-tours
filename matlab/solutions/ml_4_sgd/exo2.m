niter = 100000;
nsamples = 10;
err_rate = 50;
ElistS = [];
for is=1:nsamples
    w = zeros(p,1);
    for it=1:niter
        if mod(it,err_rate)==1
            ElistS( 1+(it-1)/err_rate,is) = E(w,X,y);
        end
        tau = tau0 / (1+it/l0);
        i = 1+floor(rand*n); % draw uniformly
        w = w - tau * nablaEi(w,i);
    end
end
clf;
subplot(2,1,1);
%
hold on;
plot(1:err_rate:niter, ElistS, 'b');
plot(1+(0:length(Elist)-1)*n, Elist, 'k--');
axis([1 niter min(ElistS(:)) max(ElistS(:))]); box on;
title('E(w_l)'); set(gca, 'FontSize', fs);
%
subplot(2,1,2);
hold on;
u = log10(ElistS-min(Elist));
v = log10(Elist -min(Elist));
plot(1:err_rate:niter, u, 'b');
plot(1+(0:length(Elist)-1)*n, v, 'k--');
axis([1 niter min(u(:)) max(u(:))]); box on;
title('log(E(w_l) - min E)'); set(gca, 'FontSize', fs);
