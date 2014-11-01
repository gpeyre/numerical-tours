def exo1():
    """
    Display the result $g_{\tau(m)}$ for different values of $m$.
    """
    mlist = N*[10 5 2 1]
    for i in 1: length(mlist):
        subplot(4, 1, i)
        m = mlist(i)
        S1 = sparse(Gamma(S, tau(m)))
        plot(PsiS(S1*Psi(f)))
        axis tight
        title(['m/ N = ' num2str(m/ N)])


def exo2():
    """
    Study the speed gain as a function of $m$ of using the sparse multiplication with respect
    to the direct computation of $T f$.
    o correction for this exercise.
    """


