klist = [10 30 50];
P = 200;
ntrials = 200;
tmin = 0; tmax = 2.5;
q = 50;
t = linspace(tmin,tmax,q);
t1 = linspace(tmin,tmax,1000);
dt = (tmax-tmin)/q;
clf;
for j=1:length(klist)
    k = klist(j);
    % simulation    
    v = [];
    for i=1:ntrials
        v = [v; svd(randn(P,k)/sqrt(P)).^2];
    end
    subplot(length(klist),1,j);
    hold on;
    h = hist(v,t);
    h = h/sum(h)/dt;
    bar(t,h);
    % theoritical law
    beta = k/P;
    a = (1-sqrt(beta))^2;
    b = (1+sqrt(beta))^2;
    z = sqrt( max(t1-a,0).*max(b-t1,0) )./(2*pi*beta*t1);
    hh = plot(t1, z, 'r' ); set(hh, 'LineWidth', 2);
    
    axis([tmin tmax 0 max(h)*1.05]); % max([z(:);h(:)])]);
    
    title(sprintf('P=%d, k=%d',P,k ));
end
