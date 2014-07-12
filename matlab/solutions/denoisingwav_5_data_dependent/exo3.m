lmax = 10;
ntest = 100000;
[V,U] = meshgrid(1:lmax,ones(ntest,1));
W = poissrnd(V);
Mstab1 = 2*sqrt(W+3/8);
Mstab2 = sqrt(W+1)+sqrt(W);
vstab = [];
vstab(1,:) = std(Mstab1,1).^2;
vstab(2,:) = std(Mstab2,1).^2;
v = std(W,1).^2;
% plot
clf;
%subplot(2,1,1);
%h = plot(v); axis('tight');
%if using_matlab()
%    set(h, 'LineWidth', 2);
%end
%title('Poisson variance');
%set_label('\lambda', 'Variance');
%subplot(2,1,2);
h = plot(vstab'); axis('tight');
if using_matlab()
    set(h, 'LineWidth', 2);
    legend('Anscombe', 'Freeman & Tukey');
end
set_label('\lambda', 'Variance');
