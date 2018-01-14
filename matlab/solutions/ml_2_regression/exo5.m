lmax =  max(X0'*y0); % rescaling factor
clf; hold on;
axis tight;
set(gca, 'FontSize', 15);
xlabel('\lambda/|X^* y|_\infty'); ylabel('w_i');
for i=1:p
    lgd{i} = class_names{i};
end
% plot(lambda_list, W', 'LineWidth', 2);
plot(lambda_list/lmax, W', 'LineWidth', 2);
plot( lambda0/lmax*[1 1], [min(W(:)) max(W(:))], 'r--', 'LineWidth', 2);
legend(lgd);
box on;
