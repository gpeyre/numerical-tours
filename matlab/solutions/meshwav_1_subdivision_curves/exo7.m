wlist = [.5 1 1.5 2]/16;
subd = @(f,w)f;
for i=1:5
    subd = @(f,w)subdivide(subd(f,w),h4pt(w));
end
lgd = {};
F = [];
for i=1:length(wlist);
    w = wlist(i);
    F(:,end+1) = subd(f0,w);
    lgd{i} = ['w=' num2str(w*16) '/16'];
end
clf;
hold on;
plot(F([1:end 1],:), 'LineWidth', 2);
myplot(f0, 'r.--');
legend(lgd);
axis tight; axis off; axis equal;
