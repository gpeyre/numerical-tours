figure(figsize = (10, 5))
Y, X = meshgrid(collect(1 : n), collect(1 : n))
k = 4 #number of circles
r = .3*n/k
phi0 = zeros(n, n) + Inf

for i in 1 : k
    for j in 1 : k
        c = ([i,j] - 1).*(n/k) + (n/k)*.5
        phi0 = min(phi0, sqrt(abs(X - c[1]).^2 + abs(Y - c[2]).^2) - r)
    end
end
        
subplot(1, 2, 1)
plot_levelset(phi0, 0)
subplot(1, 2, 2)
plot_levelset(phi0, 0, f0)