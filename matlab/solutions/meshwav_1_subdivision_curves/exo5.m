wlist = [1 2 3 6];
subd = @(f,w)f;
for i=1:5
    subd = @(f,w)subdivide(subd(f,w),hcc(w));
end
lgd = {};
F = [];
for i=1:length(wlist);
    w = wlist(i);
    F(:,end+1) = subd(f0,w);
    lgd{i} = ['w=' num2str(w*16)];
end
clf;
hold on;
plot(F([1:end 1],:), 'LineWidth', 2);
myplot(f0, 'r.--');
legend(lgd);
myaxis(.03);
