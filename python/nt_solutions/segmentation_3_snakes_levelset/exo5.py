gD = grad(phi, order=2)
d = np.maximum(eps*np.ones([n,n]), np.sqrt(np.sum(gD**2, 2)))
g = gD/np.repeat(d[:,:,np.newaxis], 2, 2)
G = - W*d*div(g[:,:,0], g[:,:,1], order=2) - np.sum(gW*gD,2)