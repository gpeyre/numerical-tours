clf; hold on;
plot(lambda_list, W', 'LineWidth', 2);
plot( lambda0*[1 1], [min(W(:)) max(W(:))], 'r--', 'LineWidth', 2);
axis tight;
set(gca, 'FontSize', 15);
xlabel('\lambda'); ylabel('w_i');
for i=1:p
    lgd{i} = class_names{i};
end
legend(lgd); box on;
