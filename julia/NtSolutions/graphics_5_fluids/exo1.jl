figure(figsize = (10,10))
rho_list = [2,4,8,16]

for i in 1:length(rho_list)
    rho = rho_list[i]
    imageplot(W(f, rho*U), "rho = $rho", [2,2,i])
end
