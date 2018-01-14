clf
for i in arange(0,p):
    plot(lambda_list/lmax, W[i,:], label=class_names[0][i])
plot( [lambda0/lmax,lambda0/lmax], [W.flatten().min(),W.flatten().max()], 'r--')
axis('tight')
xlabel('$\lambda/|X^* y|_\infty$')
ylabel('$w_i$')
legend()
