R = zeros(2, n)
R[:,B] = Z
Y = (L1 \ R')'
plot_mesh([Y;zeros(1,n)],F, lwdt=1, c="lightgrey");
