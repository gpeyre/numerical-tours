def exo1():
    """
    Check that each subsampled gradient $L_i$ is indeed a tight frame, and
    determine the value of $b$.
    You can for instance apply the operators to random vector fields.
    """
    randu = randn(N/ 2, N/ 2, 2); % a random vector field
    b = 2
    for i in 1: 4:
        LLs_u = Li{i}(LiS{i}(randu))
        % relative error should be very small
        norm(abs(LLs_u(: ) - b*randu(: )))/ norm(randu(: ))


def exo2():
    """
    The parameter $ \la $ does not appear explicitely in the iterations of the generalized forward-backward algorithm. Where does it step in ?
    t scales the functionals $G_i$ and $\text{prox}_{n \gamma G_i}$
    hould be replaced by $\text{prox}_{n \lambda\gamma G_i}$.
    """


def exo3():
    """
    Compute 100 iterations of the generalized forward-backward,
    while monitoring the value $E(x_t,\la)$ of the objective
    at iteration $t$. Display the evolution of the objective along
    iterations: it must decrease.
    = plot(log10(ObjList(1:round(end*.7))-min(ObjList)));
    itle( 'log_{10}(E(x_t,\lambda) - E(x^{*},\lambda))' );
    """
    nIter = 100
    ObjList = zeros(1, nIter)
    for it in 1: nIter:
    	ObjList(it) = E(x, lambda)
    	forward = 2*x - gamma*nablaF(x)
    for i in 1: 4:
    		z(: , : , i) = z(: , : , i) + proxGi(forward-z(: , : , i), n*gamma*lambda, i) - x
    	x = mean(z, 3)
    h = plot(ObjList)
    set(h, 'LineWidth', 2)
    xlabel('t')
    title('E(x_t, \lambda)')
    axis tight


def exo4():
    """
    Display the resulting SNR as a function of $\la$.
    Take the best regularization parameter and display the
    corresponding recovery.
    """
    nIter = 100
    lambdaListLen = length(lambdaList)
    SNRlist = zeros(1, lambdaListLen)
    SNRmax = 0
    bestLambda = 0
    for l in 1: lambdaListLen:
    	lambda = lambdaList(l)
    	x = zeros(N, N); % initialization minimizer
    	z = zeros(N, N, n); % initialization auxiliary variables
    for it in 1: nIter:
            	forward = 2*x - gamma*nablaF(x)
    for i in 1: 4:
    		    z(: , : , i) = z(: , : , i) + proxGi(forward-z(: , : , i), n*gamma*lambda, i) - x
    		x = mean(z, 3)
    	SNRlist(l) = snr(f0, x)
    	if SNRlist(l) > SNRmax
    		recov = x
    		SNRmax = SNRlist(l)
    		bestLambda = lambda
    h = semilogx(lambdaList, SNRlist)
    set(h, 'LineWidth', 2)
    xlabel('\lambda')
    ylabel('SNR')
    axis tight


def exo5():
    """
    Solve the composite regularization model. Keep the previous value of
    $\la$, set $\mu = 10^{-3}$,
    and perform 500 iterations. Display the results and compare
    visually to the previous one. Is the SNR significantly improved ?
    Conclude on the SNR as a quality criterium, and on the usefulness of mixing different regularizations priors.
    """
    nIter = 500
    lambda = bestLambda
    mu = 1e-3
    x = zeros(N, N); % initialization minimizer
    z = zeros(N, N, n); % initialization auxiliary variables
    for it in 1: nIter:
    	forward = 2*x - gamma*nablaF(x)
    for i in 1: 4:
    		z(: , : , i) = z(: , : , i) + proxGi(forward-z(: , : , i), n*gamma*lambda, i) - x
    	z(: , : , 5) = z(: , : , 5) + proxG5(forward-z(: , : , 5), n*gamma*mu) - x
    	x = mean(z, 3)
    imageplot(x)
    title(sprintf('\\lambda = %.1e; \\mu = %.1e; SNR = %.2fdB', lambda, mu, snr(f0, x)))


