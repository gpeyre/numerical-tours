def exo1():
    """
    Implement the full homotopy algorithm.
    Record in |X(:,i)|, |Lambda(i)|, |Sparsity(i)| respectively
    the solution $x_\lambda$, the value of $\lambda$
    and the size of the support $I$ at step |i| of the algorithm.
    or bookeeping
    nitialization
    """
    X = []; Lambda = []; Sparsity = []
    C = Phi'*y
    [lambda, I] = max(abs(C))
    x = zeros(N, 1)
    X(: , end + 1) = x; Lambda(end + 1) = lambda; Sparsity(end + 1) = length(I)
    niter = 1000
    for it in 1: niter:
        % complementary support
        J = setdiff(1: N, I)
        % correlations
        c = Phi'*(y-Phi*x)
        if norm(c) <1e-6
            break
        if it <= 6
            % Display to check that x_lambda is indeed a solution
            subplot(2, 3, it)
            hold on
            stem(I, c(I)/ lambda, 'b.')
            stem(J, c(J)/ lambda, 'r.')
            plot([1 N], [1 1], 'k--')
            plot([1 N], -[1 1], 'k--')
            axis([1 N -1.05 1.05])
            title(['Step ' num2str(it)])
        % update direction
        d = zeros(N, 1)
        d(I) = (Phi(: , I)'*Phi(: , I)) \ sign(c(I))
        % useful vector
        v = Phi(: , I)*d(I)
        % Compute minimum |gamma| so that situation 1) is in force.
        w = (lambda-c(J)) ./ (1 - Phi(: , J)'*v)
        gamma1 = min(w(w >0))
        if not(isempty(gamma1))
            i1 = J(w = =gamma1)
        % Compute minimum |gamma| so that situation 2) is in force.
        w = (lambda + c(J)) ./ (1 + Phi(: , J)'*v)
        gamma2 = min(w(w >0))
        if not(isempty(gamma2))
            i2 = J(w = =gamma2)
        % Compute minimum |gamma| so that situation 3) is in force.
        w = -x(I)./ d(I)
        gamma3 = min(w(w >0))
        if not(isempty(gamma3))
            i3 = I(w = =gamma3)
        % any condition is in force
        gamma = min([gamma1 gamma2 gamma3])
        if isempty(gamma)
            break;  % final solution reached
        % new solution
        x = x + gamma*d
        lambda = lambda - gamma
        if gamma = =gamma1
            I = [I i1]
        elseif gamma = =gamma2
            I = [I i2]
        elseif gamma = =gamma3
            I(I = =i3) = []
            x(i3) = 0
        % record sparsity and lambda
        X(: , end + 1) = x; Lambda(end + 1) = lambda
        Sparsity(end + 1) = length(I)


def exo2():
    """
    Test the algorithm for $x_0$ of varying sparsity $s$.
    What can you conclude in term of noise robustness ?
    """
    slist = round([P/ 16 P/ 8 P/ 4 P/ 2])
    slist = round([2 3 6 10])
    for it_s in 1: length(slist):
        s = slist(it_s)
        %
        I = randperm(N); I = I(1: s)
        x0 = zeros(N, 1)
        % x0(1: s) = 1
        x0(I) = sign(randn(s, 1))
        sigma = 0.1 * norm(Phi*x0)/ sqrt(P)
        y = Phi*x0 + sigma*randn(P, 1)
        %
        [X, Lambda] = perform_homotopy(Phi, y)
        % Fuchs criteria
        J = setdiff(1: N, I)
        U = Phi(: , J)' * pinv(Phi(: , I))'
        c = U * sign(x0(I))
        Fuchs = max(abs(c))
        % ERC 
        ERC = norm(U, 'inf')
        %
        subplot(2, 2, it_s)
        title(['s = ' num2str(s), ', F = ' num2str(Fuchs, 2) ', ERC = ' num2str(ERC, 2)])
        I = find(x0~ = 0); J = setdiff(1: N, I)
     	hold on
        h = plot(Lambda, X(I, : )', '-'); set(h, 'LineWidth', 2)
        h = plot(Lambda, X(J, : )', 'k-'); set(h, 'LineWidth', 2)
        axis tight
        box on


