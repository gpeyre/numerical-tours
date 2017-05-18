figure(figsize = (8, 5))

names = ["regular3", "phantom", "lena", "mandrill"]

for i in 1 : size(fList)[3]
    fW = perform_wavortho_transf(fList[: , : , i], Jmin, + 1, h)
    cR = sort(abs(fW)[:])[end : -1 : 1]
    err = collect(e for e in vecnorm(fList[:, :, i]).^2 - cumsum(cR.^2))
    Err = err[10 : Int(round(n*n/10))]
    plot(log10(10 : Base.div(n*n, 10)), log10(Err/Err[1]), linewidth = 2, label = names[i])
end
    
    
title(L"$\log_{10}(\epsilon^2[M]/ ||f||^2)$")
xlim(1, log10(Base.div(n*n, 10)))
ylim(-7, 0)
legend(loc = 3)
show()
