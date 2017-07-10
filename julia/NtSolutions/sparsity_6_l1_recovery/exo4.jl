klist = [10, 30, 50]
P = 200
ntrials = 200
tmin = 0
tmax = 2.5
q = 50
t = linspace(tmin, tmax, q)
t1 = linspace(tmin, tmax, 1000)
dt = (tmax - tmin)/q

for j in 1 : length(klist)
    k = klist[j]
    
    # simulation    
    v = []
    for i in 1 : ntrials
        v = [v; svd(randn(P, k)./sqrt(P))[2].^2]
    end
    
    figure(figsize = (10, 10))
    subplot(length(klist), 1, j)
    h = hist(v, t)[2]
    h = h/sum(h)/dt
    h = [h; 0]
    bar(t[1 : end], h, width = 1/20, color = "darkblue", edgecolor = "black")
    
    # theoritical law
    beta = k/P
    a = (1 - sqrt(beta))^2
    b = (1 + sqrt(beta))^2
    z = sqrt(max(t1 - a, zeros(length(t1))).*max(b - t1, zeros(length(t1))))./(2*pi.*beta.*t1)
    
    plot(t1, z, "r", linewidth = 3)
    xlim(tmin, tmax)
    ylim(0, maximum(h)*1.05)
    title(@sprintf("P = %d, k = %d", P, k))
    
    show()
end