g = (C, I) -> sum(C[:, I], 2)

figure(figsize = (8, 5))
dlist = Array{Int64,1}(1 : N/20 - 1)
criter = zeros(length(dlist), 3)

for i in 1 : length(dlist)
    s = twosparse(dlist[i])
    I = (supp(s))
    criter[i, :] = [F(Phi, s), erc(Phi,I), werc(Phi,I)]
end
    
criter[criter .< 0] = float("inf")

plot(dlist, criter, linewidth = 2)
plot(dlist, dlist.*0 + 1, "k--", linewidth = 2)
xlim(1, maximum(dlist))
ylim(minimum(criter), maximum(criter))
xlabel("d")
legend(["F", "ERC", "w-ERC"])
           
show()
