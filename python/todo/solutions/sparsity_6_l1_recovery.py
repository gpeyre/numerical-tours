def exo1():
    """
    Display how the average coherence of a random matrix
    decays with the redundancy $\eta = P/N$ of
    the matrix $\Phi$. Can you derive an empirical law between
    $P$ and the maximal sparsity?
    = plot(eta_list, k_mean ); set(h, 'LineWidth', 2);
    abel('\eta'); ylabel('1/2(1+1/E(\mu))');
    """
    eta_list = linspace(.1, 1, 10)
    ntrials = 20
    mu_mean = []
    mu_std = []
    N = 500
    for i in 1: length(eta_list):
        eta = eta_list(i)
        P = round(eta*N)
        c = []
    for k in 1: ntrials:
            c(end + 1) = mu(PhiRand(P, N))
        mu_mean(end + 1) = mean(c)
        mu_std(end + 1) = std(c)
    k_mean = 1/ 2*(1 + 1./ mu_mean)
    subplot(2, 1, 1)
    h = plot(eta_list, mu_mean); set(h, 'LineWidth', 2)
    axis tight
    xlabel('\eta'); ylabel('\mu')
    subplot(2, 1, 2)
    h = plot(log10(eta_list), log10(k_mean)); set(h, 'LineWidth', 2)
    xlabel('log_{10}(\eta)'); ylabel('log_{10}(1/ 2(1 + 1/ \mu))')
    axis tight


def exo2():
    """
    Show that this inequality holds on a given matrix.
    What can you conclude about the sharpness of these criteria ?
    """
    N = 2000
    P = N-10
    Phi = PhiRand(N, P)
    s = zeros(N, 1); s(1: 6) = 1
    I = supp(s)
    k = length(I)
    fprintf('N = %d, P = %d, |I| = %d\n', N, P, k)
    fprintf('F(s)    = %.2f\n', F(Phi, s))
    fprintf('ERC(I)  = %.2f\n', erc(Phi, I))
    fprintf('w-ERC(s) = %.2f\n', werc(Phi, I))
    fprintf('Coh(|s|) = %.2f\n', Coh(Phi, k))


def exo3():
    """
    For a given matrix $\Phi$ generated using |PhiRand|, draw as a function of the sparsity $k$
    the probability that a random sign vector $s$ of sparsity
    $\norm{s}_0=k$ satisfies the conditions $\text{F}(x_0)<1$,
    $\text{ERC}(x_0)<1$ and $\text{w-ERC}(x_0)<1$
    """
    N = 600
    P = N/ 2
    Phi = PhiRand(P, N)
    klist = round(linspace(1, P/ 7, 20))
    ntrials = 60
    proba = []
    for i in 1: length(klist):
        proba(i, 1: 3) = 0
        k = klist(i)
    for j in 1: ntrials:
            s = zeros(N, 1); I = randperm(N); I = I(1: k)
            s(I) = sign(randn(k, 1))
            proba(i, 1) = proba(i, 1) + (F(Phi, s) <1)
            proba(i, 2) = proba(i, 2) + (erc(Phi, I) <1)
            proba(i, 3) = proba(i, 3) + (werc(Phi, I) >0 & werc(Phi, I) <1)
    h = plot(klist, proba/ ntrials); set(h, 'LineWidth', 2)
    xlabel('k');  legend('F <1', 'ERC <1', 'w-ERC <1')
    title(sprintf('N = %d, P = %d', N, P))
    axis tight


def exo4():
    """
    Display, for an increasing value of $k$ the histogram of repartition
    of the eigenvalues $A^* A$ where $A$ is a Gaussian matrix of size $(P,k)$ and
    variance $1/P$. For this, accumulate the eigenvalues for many
    realization of $A$.
    """
    klist = [10 30 50]
    P = 200
    ntrials = 200
    tmin = 0; tmax = 2.5
    q = 50
    t = linspace(tmin, tmax, q)
    t1 = linspace(tmin, tmax, 1000)
    dt = (tmax-tmin)/ q
    for j in 1: length(klist):
        k = klist(j)
        % simulation    
        v = []
    for i in 1: ntrials:
            v = [v; svd(randn(P, k)/ sqrt(P)).^2]
        subplot(length(klist), 1, j)
        hold on
        h = hist(v, t)
        h = h/ sum(h)/ dt
        bar(t, h)
        % theoritical law
        beta = k/ P
        a = (1-sqrt(beta))^2
        b = (1 + sqrt(beta))^2
        z = sqrt(max(t1-a, 0).*max(b-t1, 0))./ (2*pi*beta*t1)
        hh = plot(t1, z, 'r'); set(hh, 'LineWidth', 2)
        axis([tmin tmax 0 max(h)*1.05]); % max([z(: ); h(: )])])
        title(sprintf('P = %d, k = %d', P, k))


def exo5():
    """
    Estimate numerically lower bound on $\de_k^1,\de_k^2$ by Monte-Carlo
    sampling of sub-matrices.
    """
    N = 600
    P = N/ 2
    Phi = PhiRand(P, N)
    klist = round(linspace(1, P/ 15, 10))
    ntrials = 2000
    rip_val = []
    for i in 1: length(klist):
        rip_val(i, 1: 2) = 0
        k = klist(i)
    for j in 1: ntrials:
            I = randperm(N); I = I(1: k)
            [a, b] = ric(Phi(: , I))
            rip_val(i, 1: 2) = max(rip_val(i, 1: 2), [a b])
    hold on
    h = plot(klist, rip_val);  set(h, 'LineWidth', 2)
    h = plot(klist, klist*0 + sqrt(2)-1, 'r: ');  set(h, 'LineWidth', 2)
    legend('\delta^2_k', '\delta^2_k', '0.41'); xlabel('k')
    title(sprintf('N = %d, P = %d', N, P))
    axis tight


def exo6():
    """
    Plot the evolution of the criteria F, ERC and Coh as a function of $d$.
    Do the same plot for other signs patterns for $x_0$.
    Do the same plot for a Dirac comb with a varying spacing $d$.
    """
    criter = []
    dlist = 1: N/ 20
    for i in 1: length(dlist):
        s = twosparse(dlist(i))
        I = supp(s)
        criter(i, : ) = [F(Phi, s) erc(Phi, I) werc(Phi, I)]
    criter(criter <0) = Inf
    clf; hold on
    h = plot(dlist, criter); set(h, 'LineWidth', 2)
    h = plot(dlist, dlist*0 + 1, 'k--'); set(h, 'LineWidth', 2)
    xlabel('d');  legend('F', 'ERC', 'w-ERC')
    axis([1 max(dlist) min(criter(: )) 1.5])
    axis tight


