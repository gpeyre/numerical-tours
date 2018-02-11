def exo1():
    """
    Find the optimal solution |fL2| by testing several value of
    |lambda|.
    """
    eps_list = linspace(2*1e-3, .03, 40)
    err = []
    for i in 1: length(eps_list):
        lambda = eps_list(i)
        fL2 = real(ifft2(yF .* hF ./ (abs(hF).^2 + lambda)))
        err(i) = snr(f0, fL2)
    hh = plot(eps_list, err); axis tight
    set_label('lambda', 'SNR')
    if using_matlab()
        set(hh, 'LineWidth', 2)
    [tmp, i] = max(err)
    lambda = eps_list(i)
    fL2 = real(ifft2(yF .* hF ./ (abs(hF).^2 + lambda)))


def exo2():
    """
    Find the optimal solution |fSob| by testing several value of
    |lambda|.
    """
    lambda_list = linspace(.03, .2, 40)
    err = []
    for i in 1: length(lambda_list):
        lambda = lambda_list(i)
        fSob = real(ifft2(yF .* hF ./ (abs(hF).^2 + lambda*S)))
        err(i) = snr(f0, fSob)
    hh = plot(lambda_list, err); axis tight
    set_label('lambda', 'SNR')
    if using_matlab()
        set(hh, 'LineWidth', 2)
    [tmp, i] = max(err)
    lambda = lambda_list(i)
    fSob = real(ifft2(yF .* hF ./ (abs(hF).^2 + lambda*S)))


def exo3():
    """
    Perform the deblurring by a  gradient descent.
    Keep track of the function being minimized.
    isplay energy
    """
    tau = 1.9 / (1 + lambda * 8 / epsilon)
    fTV = y
    E = []
    for i in 1: niter:
        % Compute the gradient of the smoothed TV functional.
        Gr = grad(fTV)
        d = sqrt(epsilon^2 + sum3(Gr.^2, 3))
        G = -div(Gr./ repmat(d, [1 1 2]))
        % step
        e = Phi(fTV, h)-y
        fTV = fTV - tau*(Phi(e, h) + lambda*G)
        % energy
        E(i) = 1/ 2*norm(e, 'fro')^2 + lambda*sum(d(: ))
    plot(E); axis('tight')
    set_label('Iteration #', 'Energy')


def exo4():
    """
    Explore the different values of |lambda| to find the optimal solution.
    Display the SNR as a function of |lambda|.
    """
    niter = 400
    lambda_list = linspace(1e-6, .01, 20)
    tau = 1.9 / (1 + max(lambda_list) * 8 / epsilon)
    fBest = y; fTV = y
    err = []
    for it in 1: length(lambda_list):
        lambda = lambda_list(it)
    for i in 1: niter:
            % Compute the gradient of the smoothed TV functional.
            Gr = grad(fTV)
            d = sqrt(epsilon^2 + sum3(Gr.^2, 3))
            G = -div(Gr./ repmat(d, [1 1 2]))
            % step
            e = Phi(fTV, h)-y
            fTV = fTV - tau*(Phi(e, h) + lambda*G)
        err(it) = snr(f0, fTV)
        if err(it) >snr(f0, fBest)
            fBest = fTV
    plot(lambda_list, err)
    axis('tight')
    xlabel('\lambda'); ylabel('SNR')
    fTV = fBest


def exo5():
    """
    Compare sparsity, Sobolev and TV deblurring.
    """


