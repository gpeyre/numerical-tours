def exo1():
    """
    Compute and display the histogram of $f$ for an increasing number of bins.
    """
    plist = [10 30 100 1000]
    for i in 1: length(plist):
        subplot(length(plist), 1, i)
        [h, t] = hist(f(: ), plist(i))
        bar(t, h*plist(i)/ n^2)
        axis([0 1 0 2.5])


def exo2():
    """
    Compare the two histograms.
    """
    Q = 100
    [h0, t] = hist(f(: ), Q)
    [h1, t] = hist(g(: ), Q)
    subplot(2, 1, 1)
    bar(t, h0*Q/ n^2); axis([0 1 0 6])
    subplot(2, 1, 2)
    bar(t, h1*Q/ n^2); axis([0 1 0 6])


def exo3():
    """
    Display the progression of the interpolation of the histograms.
    """
    p = 100
    tlist = linspace(0, 1, 5)
    for i in 1: length(tlist):
        a = ft(tlist(i))
        subplot(length(tlist), 1, i)
        [h, t] = hist(a(: ), p)
        bar(t, h*p/ n^2)
        axis([0 1 0 6])


