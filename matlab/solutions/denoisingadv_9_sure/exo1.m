ntrials = 100; nlaunch = 20;
E0 = []; E = [];
for i=1:nlaunch
	f = repmat(f0, [1 1 ntrials]) + sigma*randn(n,n,ntrials);
	hf = h(f,lambda);
    %
    e = sum(sum( (hf - repmat(f0, [1 1 ntrials])).^2 )); 
    E0 = [E0; e(:)];
    %
    e = -N*sigma^2 + ...
        sum(sum( (hf - f).^2 )) + ...
        2 * sigma^2 * df(lambda); E = E(:);
    E = [E; e(:)];
end
% mean
v_true = mean(E0);
v_sure = mean(E);
% display
a = mean(E0) - 8*std(E0);
b = mean(E0) + 8*std(E0);
t = linspace(a,b,31);
mybar = @(e)hist(e(e>a & e<b), t);
clf;
subplot(2,1,1); hold on;
s = mybar(E0); bar(t,s); axis tight;
u = plot(v_true*[1 1], [0 max(s)], 'r'); set(u, 'LineWidth', 3);
subplot(2,1,2); hold on;
s = mybar(E); bar(t,s); axis tight;
u = plot(v_sure*[1 1], [0 max(s)], 'r'); set(u, 'LineWidth', 3);
