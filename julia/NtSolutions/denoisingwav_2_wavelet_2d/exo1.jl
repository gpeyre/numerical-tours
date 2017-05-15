Tlist = linspace(.8, 4.5, 25).*sigma;
err_soft = []; err_hard = [];
for i in 1:length(Tlist)
    aT = NtToolBox.perform_thresholding(a, Tlist[i], "hard");
	  fWav = perform_wavortho_transf(aT, Jmin, -1, h);
    err_hard = [err_hard; snr(f0,fWav)];
    aT = NtToolBox.perform_thresholding(a, Tlist[i], "soft");
    aT[1:2^Jmin, 1:2^Jmin] = a[1:2^Jmin, 1:2^Jmin];
	  fWav = perform_wavortho_transf(aT, Jmin, -1, h);
    err_soft = [err_soft; snr(f0, fWav)];
end
h1, = plot(Tlist/sigma, err_hard)
h2, = plot(Tlist/sigma, err_soft)
axis("tight")
legend([h1, h2], ["Hard", "Soft"]);
