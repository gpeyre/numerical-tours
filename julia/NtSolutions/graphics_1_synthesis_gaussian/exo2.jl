k = 200
t = collect(linspace(0,.5,k))
subplot(2,1,1)
#h,temp = plt[:hist](vec(f0),k)
bar(t,plt[:hist](vec(f0),k)[1]*k/n^2)
xlim(0,.5)
title("Input")
subplot(2,1,2)
bar(t, plt[:hist](vec(f),k)[1]*k/n^2);
xlim(0,.5)
title("Synthesized");
