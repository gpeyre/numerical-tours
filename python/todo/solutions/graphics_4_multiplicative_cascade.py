def exo1():
    """
    Perform the cascade. Display intermediate steps.
    """
    f = ones(n, 1) * exp(H1) / rmin^H1
    k = 0
    for i in 1: N:
       I = find(abs(t-ti(i)) <= ri(i)/ 2)
       f(I) = f(I) * Wi(i)
       if i = =10 || i = =80 || i = =round(.1*N) || i = =round(.5*N)
           k = k + 1
           subplot(2, 2, k)
           plot(t, f)
    		axis([0 T 0 1.1*max(f)])


def exo2():
    """
    Compute several realization for the same log-normal parameters.
    """
    for k in 1: 4:
       Wi = exp(randn(N, 1)*sqrt(sigma2) + mu)
       ti = -1/ 2 + rand(1, N) * (T + 1)
       ti_1 = -1/ 2 + rand(1, round(T + 1))*(T + 1)
       ti = [ti_1 ti];   % for exact scale invariance
       umax = 1/ rmin-1
       ui = [zeros(1, length(ti_1)) rand(1, N) * umax]; % ui = 1/ ri-1
       ri = (1 + ui).^(-1)
       f = ones(n, 1)* exp(H1) / rmin^H1
    for i in 1: N:
           I = find(abs(t-ti(i)) <= ri(i)/ 2)
           f(I) = f(I) * Wi(i)
       subplot(2, 2, k)
       plot(t, f)
       axis([0 T 0 1.1*max(f)])


def exo3():
    """
    Compute realizations for different log-normal parameters |mu| and |sigma2|.
    Use the same distribution of points.
    """
    sigma2list = [0.1 0.1 0.5 0.5]
    mult = [1 10 1 10]
    mulist = -mult.*sigma2list
    ti = -1/ 2 + rand(1, N) * (T + 1)
    ti_1 = -1/ 2 + rand(1, round(T + 1))*(T + 1)
    ti = [ti_1 ti];   % for exact scale invariance
    umax = 1/ rmin-1
    ui = [zeros(1, length(ti_1)) rand(1, N) * umax]; % ui = 1-1/ ri^{1 + beta}
    ri = (1 + ui).^(-1)
    for k in 1: 4:
       mu = mulist(k)
       sigma2 = sigma2list(k)
       if -(exp(2*(sigma2 + mu))-1) <-1  % Condition of non-degeneracy
           disp('Be careful ! This cascade will degenerate as rmin - > 0 !')
       H1 = 1 - exp(mu + sigma2/ 2)
       Wi = exp(randn(N, 1)*sqrt(sigma2) + mu)
       f = ones(n, 1)* exp(H1) / rmin^H1
    for i in 1: N:
           I = find(abs(t-ti(i)) <= ri(i)/ 2)
           f(I) = f(I) * Wi(i)
       subplot(2, 2, k)
       plot(t, f)
       axis([0 T 0 1.1*max(f)]) 
       title(['mu = ' num2str(mu) ', sigma^2 = ' num2str(sigma2)])


def exo4():
    """
    Perform the full cascade, display intermediate steps.
    """
    f = ones(n)/ rmin^H1
    k = 0
    for i in 1: N:
       %progressbar(i, N)
       I = find((X-xi(i)).^2 + (Y-yi(i)).^2 <= ri(i)^2/ 4)
       f(I) = f(I) * Wi(i)
       if i = =30 || i = =200 || i = =round(N*.5) || i = =round(N*.8)  
           k = k + 1
           subplot(2, 2, k)
           imageplot(f)


def exo5():
    """
    Compute the fractional integration for several values of alpha.
    """
    alphalist = [.1 .3 .6 1]
    for i in 1: 4:
       alpha = alphalist(i)
       F = real(ifft2(fft2(f)./ S.^alpha))
       subplot(2, 2, i)
       imageplot(F, ['alpha = ' num2str(alpha)])


def exo6():
    """
    Perform the  cascade for several log-normal parameters |mu| and |sigma2|.
    """
    sigma2list = [0.02 0.02 0.2 0.2]
    mult = [1 10 1 1]
    mulist = -mult.*sigma2list/ 2
    k = 0
    for k in 1: 4:
       mu = mulist(k)
       sigma2 = sigma2list(k)
       if -(exp(2*(sigma2 + mu))-1) <-1  % Condition of non-degeneracy
           disp('Be careful ! This cascade will degenerate as rmin - > 0 !')
       H1 = 1 - exp(mu + sigma2/ 2)
       Wi = exp(randn(N, 1)*sqrt(sigma2) + mu)
       f = ones(n)/ rmin^H1
    for i in 1: N:
           %progressbar(i, N)
           I = find((X-xi(i)).^2 + (Y-yi(i)).^2 <= ri(i)^2/ 4)
           f(I) = f(I) * Wi(i)
       subplot(2, 2, k)
       imageplot(f, ['mu = ' num2str(mu) ', sigma^2 = ' num2str(sigma2)])


