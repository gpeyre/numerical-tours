gD = grad(phi, order=2)
d = np.maximum(eps*np.ones([n,n]), np.sqrt(np.sum(gD**2, 2)))
g = gD/np.repeat(d[:,:,np.newaxis], 2, 2)
G = d*div(g[:,:,0], g[:,:,1], order=2) - lambd*(f0-c1)**2 + lambd*(f0-c2)**2