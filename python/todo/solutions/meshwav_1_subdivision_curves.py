def exo1():
    """
    Perform several step of subdivision. Display the different curves.
    """
    Jmax = 3
    f = f0
    for j in 0: Jmax:
        f = subdivide(f, h)
        subplot(2, 2, j + 1)
        hold on
        myplot(f, 'k.-')
        myplot(f0, 'r.--')
        myaxis(0)


def exo2():
    """
    Compute the scaling function $\phi$
    associated to the subdivision.
    """
    n = 6
    f = zeros(n, 1); f(n/ 2 + 1) = 1
    for i in 1: 5:
        f = subdivide(f, h)
    plot(linspace(-1/ 2, 1/ 2, length(f)), f); axis([-1/ 2 1/ 2 -.01 max(f)*1.03])


def exo3():
    """
    Test with different configurations of control points.
    """
    f0 = [0 + 0i; 1 + 0i; 1 + 1i; 0 + 1i]
    f = f0
    for j in 0: Jmax:
        f = subdivide(f, h)
        subplot(2, 2, j + 1)
        hold on
        myplot(f, 'k.-')
        myplot(f0, 'r.--')
        myaxis(.03)


def exo4():
    """
    Test the corner-cutting for $w=3$.
    """
    h = hcc(3)
    f = f0
    for j in 0: Jmax:
        f = subdivide(f, h)
        subplot(2, 2, j + 1)
        hold on
        myplot(f, 'k.-')
        myplot(f0, 'r.--')
        myaxis(.03)


def exo5():
    """
    Test the corner-cutting for vaious values of $w$.
    """
    wlist = [1 2 3 6]
    subd = lambda f, w: f
    for i in 1: 5:
        subd = lambda f, w: subdivide(subd(f, w), hcc(w))
    lgd = {}
    F = []
    for i in 1: length(wlist):
        w = wlist(i)
        F(: , end + 1) = subd(f0, w)
        lgd{i} = ['w = ' num2str(w*16)]
    hold on
    plot(F([1: end 1], : ), 'LineWidth', 2)
    myplot(f0, 'r.--')
    legend(lgd)
    myaxis(.03)


def exo6():
    """
    Perform the interpolating subdivision
    for $w=1/16$.
    """
    w = 1/ 16
    h = h4pt(w)
    f = f0
    for j in 0: Jmax:
        f = subdivide(f, h)
        subplot(2, 2, j + 1)
        hold on
        myplot(f, 'k.-')
        myplot(f0, 'r.--')
        myaxis(.13)
        hold off


def exo7():
    """
    Test the influence of $w$.
    """
    wlist = [.5 1 1.5 2]/ 16
    subd = lambda f, w: f
    for i in 1: 5:
        subd = lambda f, w: subdivide(subd(f, w), h4pt(w))
    lgd = {}
    F = []
    for i in 1: length(wlist):
        w = wlist(i)
        F(: , end + 1) = subd(f0, w)
        lgd{i} = ['w = ' num2str(w*16) '/ 16']
    hold on
    plot(F([1: end 1], : ), 'LineWidth', 2)
    myplot(f0, 'r.--')
    legend(lgd)
    axis tight; axis off; axis equal


def exo8():
    """
    Compare the result of the quadratic B-spline, cubic B-spline,
    and 4-points interpolating.
    """
    hh = []
    hh{end + 1} = [1 3 3 1]/ 4
    hh{end + 1} = [1 4 6 4 1]/ 8
    hh{end + 1} = [-1, 0, 9, 1, 9, 0, -1]/ 16
    hh{end}((end + 1)/ 2) = 1
    lgd = {'Quadratic', 'Cubic', 'Interpolating'}
    col = {'r' 'g' 'b'}
    clf; hold on
    for k in 1: length(hh):
        h = hh{k}
        f = f0
    for j in 0: 7:
            f = subdivide(f, h)
        myplot(f, col{k})
    myplot(f0, 'k.--')
    axis tight; axis off
    legend(lgd)


def exo9():
    """
    Display the scaling function associated to these Deslauriers-Dubuc filters.
    """
    n = 8
    for k in 1: 4:
        h = hdd(k)
        f = zeros(n, 1); f(n/ 2 + 1) = 1
    for i in 1: 5:
            f = subdivide(f, h)
        subplot(4, 1, k)
        plot(linspace(-n/ 2, n/ 2, length(f)), f)
        axis([-n/ 2 n/ 2 -.15 1.03])


def exo10():
    """
    Perform an approximation $f$ of the curve using a uniform sampling with $N_0=20$
    points.
    """
    p = 28
    t0 = (0: 1/ n: 1-1/ n)'
    t = (0: 1/ p: 1-1/ p)'
    f0 = interp1(t0, F, t)
    f = f0
    Jmax = ceil(log2(n/ p))
    for j in 0: Jmax:
        f = subdivide(f, h)
    clf; hold on
    myplot(F, 'k')
    myplot(f0, 'k.')
    myplot(f, 'r')
    myaxis(0)


def exo11():
    """
    Display the decay of the Hausdorff approximation error as the number $N_0$ of
    sampling points increases.
    """
    q = 15
    plist = round(linspace(8, 80, q))
    d = []
    for k in 1: length(plist):
        p = plist(k)
        % sample
        t0 = (0: 1/ n: 1-1/ n)'
        t = (0: 1/ p: 1-1/ p)'
        f0 = interp1(t0, F, t)
        % interpolate
        f = f0
        Jmax = ceil(log2(n/ p))
    for j in 0: Jmax:
            f = subdivide(f, h)
        % record distance    
        d(end + 1) = hausdorff(F, f)
    plot(plist, d, 'LineWidth', 2)
    xlabel('N_0'); ylabel('d')
    axis tight


def exo12():
    """
    Perform curve subdivision in 3D space.
    ontrol mesh.
    """
    f0 = rand(12, 3)
    Jmax = 4; ms = 20; lw = 1.5
    f = f0
    for j in 0: Jmax:
        f = cat(2, upsampling(f(: , 1)), upsampling(f(: , 2)), upsampling(f(: , 3)))
        f = cat(2, cconv(f(: , 1), h), cconv(f(: , 2), h), cconv(f(: , 3), h))
    subplot(1, 2, 1)
    hold on
    hh = plot3([f(: , 1); f(1, 1)], [f(: , 2); f(1, 2)], [f(: , 3); f(1, 3)], 'k-')
    set(hh, 'MarkerSize', ms)
    set(hh, 'LineWidth', lw)
    hh = plot3([f0(: , 1); f0(1, 1)], [f0(: , 2); f0(1, 2)], [f0(: , 3); f0(1, 3)], 'r.--')
    set(hh, 'LineWidth', lw)
    axis('tight'); box('on'); view(3);  % axis('off')
    subplot(1, 2, 2)
    hold on
    hh = plot3([f(: , 1); f(1, 1)], [f(: , 2); f(1, 2)], [f(: , 3); f(1, 3)], 'k-')
    set(hh, 'MarkerSize', ms)
    set(hh, 'LineWidth', lw)
    hh = plot3([f0(: , 1); f0(1, 1)], [f0(: , 2); f0(1, 2)], [f0(: , 3); f0(1, 3)], 'r.--')
    set(hh, 'LineWidth', lw)
    axis('tight'); box('on')
    view(70, 25)


