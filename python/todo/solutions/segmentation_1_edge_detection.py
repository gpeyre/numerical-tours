def exo1():
    """
    Test blurring with several blurring size $\si$.
    """
    slist = [1 2 5 10]
    for i in 1: length(slist):
        sigma = slist(i)
        imageplot(blur(f0, sigma), ['\sigma = ' num2str(sigma)], 2, 2, i)


def exo2():
    """
    For $\si=1$, study the influence of the threshold value $t$.
    """
    t_list = max(d(: )) * [1/ 4 1/ 5 1/ 10 1/ 20]
    for i in 1: length(t_list):
        t = t_list(i)
        imageplot(double(d >t), ['t = ' num2str(t, 2)] , 2, 2, i)


def exo3():
    """
    Study the influence of $\si$.
    """
    slist = [1 2 4 6]
    for i in 1: length(slist):
        sigma = slist(i)
        d = sqrt(sum(nabla(blur(f0, sigma)).^2, 3))
        t = max(d(: )) * 1/ 5
        imageplot(double(d >t), ['\sigma = ' num2str(sigma)], 2, 2, i)


def exo4():
    """
    Study the influence of $\si$.
    """
    slist = [4 6 10 15]
    for i in 1: length(slist):
        sigma = slist(i)
        subplot(2, 2, i)
        plot_levelset(delta(blur(f0, sigma)) , 0, f0)
        title(['\sigma = ' num2str(sigma)])


def exo5():
    """
    Study the influence of $\si$.
    """
    slist = [4 6 10 15]
    for i in 1: length(slist):
        sigma = slist(i)
        %
        g = grad(blur(f0, sigma))
        h = hessian(blur(f0, sigma))
        a = h(: , : , 1: 2).*repmat(g(: , : , 1), [1 1 2]) + ...
            h(: , : , 2: 3).*repmat(g(: , : , 2), [1 1 2])
        %
        subplot(2, 2, i)
        plot_levelset(sum(a.*g, 3) , 0, f0)
        title(['\sigma = ' num2str(sigma)])


