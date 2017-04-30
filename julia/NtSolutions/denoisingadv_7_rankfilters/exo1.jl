figure(figsize = (10, 7))

beta_list = linspace(0, 1, 6)
index = 1

for i in 1 : length(beta_list)
        beta_c = beta_list[i]
        imageplot(phi(f, beta_c)[:, :], @sprintf("Beta = %.1f", beta_c), [2, 3, i])
        index += 1
end
