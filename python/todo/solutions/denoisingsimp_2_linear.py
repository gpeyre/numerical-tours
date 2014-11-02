def exo1():
    """
    Display a denoised signal for several values of $\mu$.
    """
    mu_list = linspace(.5, 6, 6)
    for i in 1: length(mu_list):
        mu = mu_list(i)
        subplot(3, 2, i)
        plot(denoise(y, mu))
        title(['\mu = ' num2str(mu)])
        axis([1 N -.05 1.05])


def exo2():
    """
    Display the evolution of the oracle denoising error
    $ \norm{y-x_0} $ as a function of $\mu$.
    Set $\mu$ to the value of the optimal parameter.
    etrieve the best denoising result
    """
    mulist = linspace(.1, 3.5, 31)
    err = arrayfun(lambda mu: norm(x0-denoise(y, mu), 'fro'), mulist)
    h = plot(mulist, err); axis('tight')
    set_label('\mu', '|y-x_0|')
    [~, i] = min(err)
    mu = mulist(i)


