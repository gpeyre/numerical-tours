def exo1():
    """
    Compute the rank filter for several values of $\beta$.
    """
    beta_list = linspace(0, 1, 6)
    for i in 1: length(beta_list):
        beta = beta_list(i)
        imageplot(phi(f, beta), ['\beta = ' num2str(beta)], 2, 3, i)


def exo2():
    """
    Compute a closing followed by an opening.
    """
    closingopening = lambda f: opening(closing(f))
    imageplot(closingopening(f))


def exo3():
    """
    Compute an opening followed by a closing.
    """
    openingclosing = lambda f: closing(opening(f))
    imageplot(openingclosing(f))


def exo4():
    """
    Perform iterated opening and closing.
    """
    f1 = f
    for i in 1: 4:
        f1 = opening(f1)
        imageplot(f1, ['iteration ' num2str(i)], 2, 4, i)
    f1 = f
    for i in 1: 4:
        f1 = closing(f1)
        imageplot(f1, ['iteration ' num2str(i)], 2, 4, 4 + i)


def exo5():
    """
    Perform iterated median filtering, and store the output in |f1|.
    """
    f1 = f
    for i in 1: 6:
        f1 = medfilt(f1)
        imageplot(f1, ['iteration ' num2str(i)], 2, 3, i)


