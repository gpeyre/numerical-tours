from numpy import *
from matplotlib.pyplot import *

from nt_toolbox.general import *
from nt_toolbox.signal import *
from nt_toolbox.graph import *

def exo1(x0,W):

    """
        Implement the Dijkstra algorithm by iterating these step while the
        stack |I| is non empty.
        Display from time to time the front that propagates.
    """
    n = W.shape[0]
    pstart = transpose(array([x0]))
    [D,Dsvg,Ssvg] = perform_dijstra_fm(W, pstart, inf,'dijstr', 'sym',n*6)
    plt.clf;
    for i in arange(0,4):
        plt.subplot(2, 2, i+1)
        d = Dsvg[:,:,i]
        d[d==inf] = 0
        imageplot(d)
        set_cmap('jet')
    return D

def exo2(x0,W):

    """
        Implement the FM algorithm by iterating these step while the
        stack |I| is non empty.
        Display from time to time the front that propagates.
    """
    n = W.shape[0]
    pstart = transpose(array([x0]))
    [D,Dsvg,Ssvg] = perform_dijstra_fm(W, pstart, inf,'fm', 'sym',n*6)
    clf;
    for i in arange(0,4):
        subplot(2, 2, i+1)
        d = Dsvg[:,:,i]
        d[d==Inf] = 0
        imageplot(d)
        set_cmap('jet')
    return D

def exo3(x0,W):
    """
    Compute the distance map to these starting point using the FM algorithm.
    """
    n = W.shape[0]
    pstart = transpose(array([x0]))
    [D,Dsvg,Ssvg] = perform_dijstra_fm(W, pstart, inf,'fm', 'sym',n*6)
    # display
    k = 8
    displ = lambda D: cos(2*pi*k*D/ max(D.flatten()))
    clf
    imageplot(displ(D))
    set_cmap('jet')
    return D

def exo4(tau,x0,x1,G):
    """
    Perform the full geodesic path extraction by iterating the gradient
    descent. You must be very careful when the path become close to
    $x_0$, because the distance function is not differentiable at this
    point. You must stop the iteration when the path is close to $x_0$.
    """
    n = G.shape[0]
    Geval = lambda G,x: bilinear_interpolate(G[:,:,0], imag(x), real(x) ) + 1j * bilinear_interpolate(G[:,:,1],imag(x), real(x))
    niter = 1.5*n/tau;
    # init gamma
    gamma = [x1]
    xtgt = x0[0] + 1j*x0[1]
    for i in arange(0,niter):
        g = Geval(G, gamma[-1] )
        gamma.append( gamma[-1] - tau*g )
        if abs(gamma[-1]-xtgt)<1:
            break
    gamma.append( xtgt )
    return gamma
