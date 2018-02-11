def exo1():
    """
    Compute $w_0$ and $w_1$. Compute and display the segmentation when
    $\la=0$.
    """
    compute_w = lambda I, c: sum((I - repmat(reshape(c, [1 1 3]), [n n 1])).^2, 3)
    w0 = compute_w(I, c0)
    w1 = compute_w(I, c1)
    Omega = w0 <w1
    display_segmentation = lambda u: repmat(reshape(c0, [1 1 3]), [n n 1]) .* repmat(u >.5, [1 1 3]) + ...
            repmat(reshape(c1, [1 1 3]), [n n 1]) .* repmat(u <.5, [1 1 3])
    imageplot(display_segmentation(Omega))


def exo2():
    """
    Implement the DR iterative algorithm on |niter| iterations.
    Keep track of the evolution of the minimized energy
    $$ E(f) = \dotp{w}{f} + \la \norm{\nabla f}_1 $$
    during the iterations.
    _Remark:_ to speedup the convergence, you can use a "clever"
    initialization.
    """
    tf = double(w >0)
    tu = Grad(tf)
    E = []
    lun = lambda u: sum(sum(sqrt(sum(u.^2, 3))))
    for i in 1: niter:
        tf1 = tf; tu1 = tu
        % rProxG
        [tf2, tu2] = ProxG(tf1, tu1, gamma)
        tf1 = 2*tf2 - tf1; tu1 = 2*tu2 - tu1
        % rProxF
        [tf2, tu2] = ProxF(tf1, tu1, gamma)
        tf1 = 2*tf2 - tf1; tu1 = 2*tu2 - tu1
        % Average
        tf = (1-mu/ 2)*tf + mu/ 2*tf1
        tu = (1-mu/ 2)*tu + mu/ 2*tu1
        % Converging sequence
        [f, u] = ProxG(tf, tu, gamma)
        % 
        E(end + 1) = sum(w(: ).*f(: )) + lambda*lun(u)
    plot(log10(E(1: end/ 2)-min(E)))
    axis tight
    title('log_{10}(E(f) - E^*)')


def exo3():
    """
    Test with different value of the $\lambda$ parameter.
    """
    lambda_list = [0.01 .1 .5 1]
    tf = double(w >0)
    tu = Grad(tf)
    lun = lambda u: sum(sum(sqrt(sum(u.^2, 3))))
    for k in 1: length(lambda_list):
        E = []
        %
        lambda = lambda_list(k)
        ProxF = lambda f, u, gamma: deal(ProxF0(f, gamma), ProxL1(u, gamma*lambda))
    for i in 1: niter:
            tf1 = tf; tu1 = tu
            % rProxG
            [tf2, tu2] = ProxG(tf1, tu1, gamma)
            tf1 = 2*tf2 - tf1; tu1 = 2*tu2 - tu1
            % rProxF
            [tf2, tu2] = ProxF(tf1, tu1, gamma)
            tf1 = 2*tf2 - tf1; tu1 = 2*tu2 - tu1
            % Average
            tf = (1-mu/ 2)*tf + mu/ 2*tf1
            tu = (1-mu/ 2)*tu + mu/ 2*tu1
            % Converging sequence
            [f, u] = ProxG(tf, tu, gamma)
            %
            E(end + 1) = sum(w(: ).*f(: )) + lambda*lun(u)
        subplot(2, 2, k)
        imageplot(display_segmentation(f))
        title(['\lambda = ' num2str(lambda)])


