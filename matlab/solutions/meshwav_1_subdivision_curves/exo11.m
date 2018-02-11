q = 15;
plist = round(linspace(8,80,q));
d = [];
for k=1:length(plist)
    p = plist(k);
    % sample
    t0 = (0:1/n:1-1/n)';
    t = (0:1/p:1-1/p)';
    f0 = interp1(t0,F,t);
    % interpolate
    f = f0;
    Jmax = ceil(log2(n/p));
    for j=0:Jmax
        f = subdivide(f, h);
    end
    % record distance    
    d(end+1) = hausdorff(F,f);
end
clf;
plot(plist, d, 'LineWidth', 2);
xlabel('N_0'); ylabel('d');
axis tight;
