f0w = W(f0)
ntrials = 100
nlaunch = 20
E0 = []
E = []
dfw = (fw, lambd) -> sum(abs(fw) .> lambd, (1,2))

for i in 1 : nlaunch
    
    fw = repeat(f0w, inner = [1, 1, ntrials]) + sigma.*rand(Normal(), n, n, ntrials)
    hfw = S(fw, lambd)
    
    #quadratic error
    e = sum((hfw - repeat(f0w, inner = [1, 1, ntrials])).^2, (1, 2))
    E0 = [E0; e[:]]
    
    #sure error
    e = -N*sigma^2 + sum((hfw - fw).^2, (1, 2)) + 2*sigma^2*dfw(fw, lambd)
    E = [E; e[:]]
    
end
    
v_true = mean(E0)
v_sure = mean(E)
a = v_true - 8*stdm(E0, mean(E0))
b = v_true + 8*stdm(E0, mean(E0))
t = linspace(a, b, 31)
mybar = e -> hist(e[collect((i > a) & (i < b) for i in E0)], t)

figure(figsize = (10, 7))


subplot(2,1,1)
s = mybar(E0)[2]
s = [s; 0]
bar(t[1 : end], s, width = (b-a)/31, color = "darkblue", edgecolor = "white")
axvline(v_true, color = "red", linewidth = 3)

subplot(2,1,2)
s = mybar(E)[2]
s = [s; 0]
bar(t[1 : end], s, width = (b-a)/31, color = "darkblue",edgecolor = "white")
axvline(v_sure, color = "red", linewidth = 3)

show()