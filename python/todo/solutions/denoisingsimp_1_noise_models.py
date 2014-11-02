def exo1():
    """
    What is the optimal threshold |T| to remove as much as possible of noise
    ? Try several values of |T|.
    """
    k = 100
    Tlist = linspace(0, 1, k)
    err = []
    for i in 1: k:
        x1 = x .* (abs(x) >Tlist(i))
        err(i) = norm(x0-x1, 'fro')
    plot(Tlist, err); axis('tight')
    set_label('T', 'Error')


def exo2():
    """
    The theory predicts that the maximum of |n| Gaussian variable of variance |sigma^2|
    is smaller than |sqrt(2*log(n))| with large probability (that tends to 1
    when |n| increases). This is also a sharp result. Check this numerically
    by computing with Monte Carlo sampling the maximum with |n| increasing
    (in power of 2). Check also the deviation of the maximum when you
    perform several trial with |n| fixed.
    """
    ntrial = 5000
    nlist = 2.^(5: 13)
    gauss_max = []; gauss_dev = []
    for i in 1: length(nlist):
        w = randn(nlist(i), ntrial)
        m = compute_max(abs(w), 1)
        gauss_max(i) = mean(m)
        gauss_dev(i) = std(m)
    gauss_max_th = sqrt(2*log(nlist))
    subplot(2, 1, 1)
    hold on
    plot(gauss_max_th, gauss_max, 'r.-')
    plot(gauss_max_th, gauss_max_th, 'b: ')
    hold off
    axis('tight')
    set_label('sqrt(2*log(n))', 'empirical max')
    title('Maximum of n dimensional Gaussian noise')
    subplot(2, 1, 2)
    plot(gauss_max_th, gauss_dev, '.-')
    axis('tight')
    set_label('sqrt(2*log(n))', 'Deviation from max')


