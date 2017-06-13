q = floor(p/4)
t = collect(0:(q-1))'/q
t1 = collect(0:(p-3*q-1))'/(p-3*q)
Z = [[t t*0+1 1-t t1*0]; [t*0 t t*0+1 1-t1]]

figure(figsize=(10,10))
axis("off")
xlim(-.1,1.1)
ylim(-.1,1.1)
plot(Z[1, [1:p; 1]], Z[2, [1:p; 1]], ".-", c="blue", linewidth=2, markerfacecolor="red", markeredgecolor="red", ms=10)
