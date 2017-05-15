MW = perform_haar_transf(Mnoisy, 1, +1);
Tlist = linspace(1,4,20).*sigma;
err_hard = []; err_soft = [];
for i in 1:length(Tlist)
    MWT = perform_thresholding(MW, Tlist[i], "hard");
    M1 = perform_haar_transf(MWT, 1, -1);
    err_hard[i] = snr(M, M1);
    MWT = perform_thresholding(MW, Tlist[i], "soft");
    M1 = perform_haar_transf(MWT, 1, -1);
    err_soft[i] = snr(M, M1);
    if i>1
      if err_soft[i] > maximum(err_soft[1:i-1])
        Mwav = M1;
      end
    end
end
figure(figsize = (7, 5))

plot(Tlist./sigma, err_hard, ".-", label = "hard", color = "b")
plot(Tlist./sigma, err_soft, ".-", label = "soft", color = "r")
xlabel(L"$T/\sigma$")
ylabel("SNR")
ylim(minimum(err_hard), maximum(err_soft))
legend(loc = "upper right")

show()
