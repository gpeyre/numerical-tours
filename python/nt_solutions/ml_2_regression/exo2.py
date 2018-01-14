clf
for i in arange(0,p):
    plot(lambda_list/lmax, W[i,:], label=class_names[0][i])
axis('tight')
xlabel('$\lambda/|X|^2$')
ylabel('$w_i$')
legend()
