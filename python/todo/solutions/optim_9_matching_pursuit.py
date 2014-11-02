def exo1():
    """
    Perform the Matching pursuit. Record the evolution of the solution $x^{(\ell)}$
    as well as the maximum correlation $ \la^{(\ell)} = \abs{c_i} $.
    Display the decay of the error $ E(\ell) = \norm{y-\Phi x^{(\ell)}}. $
    """
    M = 2*N
    x = zeros(N, 1)
    X = []; lambda = []; E = []
    for k in 1: M:
        E(k) = norm(y-Phi*x)
        c = Phi'*(y-Phi*x)
        [lambda(k), i] = max(abs(c))
        x(i) = x(i) + c(i)
        % record
        X(: , k) = x
    h = plot(log10(E)); axis tight
    set(h, 'LineWidth', 2)
    title('log(E)')


def exo2():
    """
    Perform the Orthogonal Matching Pursuit algorithm,
    and compare it with the Matching Pursuit.
    -- MP --- %
    -- OMP --- %
    isplay
    """
    M = 2*N
    x = zeros(N, 1)
    X_mp = []; lambda_mp = []; E_mp = []
    for k in 1: M:
        E_mp(k) = norm(y-Phi*x)
        c = Phi'*(y-Phi*x)
        [lambda_mp(k), i] = max(abs(c))
        x(i) = x(i) + c(i)
        % record
        X_mp(: , k) = x
    M = P
    x = zeros(N, 1)
    X_omp = []; lambda_omp = []; E_omp = []
    for k in 1: M:
        E_omp(k) = norm(y-Phi*x)
        c = Phi'*(y-Phi*x)
        [lambda_omp(k), i] = max(abs(c))
        x(i) = x(i) + c(i)
        % projection    
        I = find(x~ = 0)
        x(I) = pinv(Phi(: , I))*y
        % record
        X_omp(: , k) = x
    I = find(x0~ = 0); J = setdiff(1: N, I)
    subplot(2, 1, 1); hold on
    h = plot(E_mp, X_mp(I, : )', '-'); set(h, 'LineWidth', 2)
    h = plot(E_mp, X_mp(J, : )', 'k-'); set(h, 'LineWidth', 2)
    title('MP'); axis tight;  box on
    subplot(2, 1, 2); hold on
    h = plot(E_omp, X_omp(I, : )', '-'); set(h, 'LineWidth', 2)
    h = plot(E_omp, X_omp(J, : )', 'k-'); set(h, 'LineWidth', 2)
    title('OMP'); axis tight; box on


