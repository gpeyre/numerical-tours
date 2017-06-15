ntests = 20;
slist = linspace(.01, 1.5, ntests);
err = [];
Mblur = copy(Mh)
for i in 1:ntests
    h = exp( -(X.^2 + Y.^2 + Z.^2)./(2*slist[i]^2) );
    h = h./sum(h[:]);
    Mh = real( plan_ifft( (plan_fft(Mnoisy)*Mnoisy) .* (plan_fft(fftshift(h))*fftshift(h)) )*( (plan_fft(Mnoisy)*Mnoisy) .* (plan_fft(fftshift(h))*fftshift(h)) ));
    err = [err; NtToolBox.snr(M, Mh)]

    if (i > 1)
        if err[i] > maximum(err[1:i-1])
            Mblur = copy(Mh)
        end
    end
end

figure(figsize = (7, 5))
plot(slist, err, ".-")
xlabel("s")
ylabel("SNR")
show()
