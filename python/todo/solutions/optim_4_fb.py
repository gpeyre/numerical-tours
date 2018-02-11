def exo1():
    """
    Compute the solution of L1 deconvolution.
    Keep track of the degay of the energy $E=f+g$.
    isplay energy decay
    """
    niter = 4000
    x = y
    E = []
    for i in 1: niter:
        x = proxg(x - gamma*gradf(x), gamma)
        % record the energy degay.
        E(i) = 1/ 2*norm(Phi(x) - y)^2 + lambda*norm(x, 1)
    hh = plot(log10((E(1: end/ 4)-E(end))/ (E(1)-E(end))))
    set(hh, 'LineWidth', 2)
    axis('tight')
    set_label('i', 'log_{10}(E-E*)')


def exo2():
    """
    Impement the relaxed FB algorithm and display its convergence rate for several values of $\mu$.
    isplay
    """
    mu_list = [-.5 0 .5 .8 .9 .95 .99]
    Erelax = []
    for imu  in  1: length(mu_list):
        mu = mu_list(imu)
        lgd{imu} = ['\mu = ' num2str(mu)]
        x = y; z = y
    for i in 1: niter:
            xold = x
            told = t
            x = proxg(z - gamma*gradf(z), gamma)
            z = x + mu*(x-xold)
            % record the energy degay.
            Erelax(i, imu) = 1/ 2*norm(Phi(x) - y)^2 + lambda*norm(x, 1)
    Emin = min(Erelax(: ))
    clf; hold on
    hh = plot(log10(Erelax(1: end/ 4, : )-Emin))
    set(hh, 'LineWidth', 2)
    axis('tight')
    set_label('i', 'log_{10}(E-E*)')
    legend(lgd)


def exo3():
    """
    Compute the solution of L1 deconvolution using FISTA.
    Keep track of the degay of the energy $E = f+g$.
    isplay energy decay
    """
    gamma = 1/ L
    x = y
    z = y
    t = 1
    Efista = []
    for i in 1: niter:
        xold = x
        told = t
        x = proxg(z - gamma*gradf(z), gamma)
        t = (1 + sqrt(1 + 4*t^2))/ 2
        z = x + (told - 1)/ t*(x-xold)
        % record the energy degay.
        Efista(i) = 1/ 2*norm(Phi(x) - y)^2 + lambda*norm(x, 1)
    Emin = min([min(E) min(Efista) min(Erelax(: ))])
    clf; hold on
    hh = plot(log10(E(1: end/ 4)-Emin), 'b')
    set(hh, 'LineWidth', 2)
    hh = plot(log10(Erelax(1: end/ 4, end)-Emin), 'g')
    set(hh, 'LineWidth', 2)
    hh = plot(log10(Efista(1: end/ 4)-Emin), 'r')
    set(hh, 'LineWidth', 2)
    axis('tight')
    set_label('i', 'log_{10}(E-E*)')
    legend('FB', ['FB-relax, \mu = ' num2str(mu_list(end))], 'FISTA')


