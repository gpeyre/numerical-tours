    """
    Perform the full geodesic path extraction by iterating the gradient
    descent. You must be very careful when the path become close to
    $x_0$, because the distance function is not differentiable at this
    point. You must stop the iteration when the path is close to $x_0$.
    """
    niter = 1.5*n/tau;
    gamma = zeros(niter,1);
    gamma = x1
    for i in arange(0,niter):
        gamma[:,i+1] = gamma[:,i] - tau*Geval( G, gamma[:,i] )
        if norm(gamma[:,i+1]-x0)<1:
            break
    gamma[:,i+1] = x0
    gamma = gamma[:,0:i+1];