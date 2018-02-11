import numpy as np

def compute_normal(vertex,face):
    """
        compute_normal - compute the normal of a triangulation
        
          [normal,normalf] = compute_normal(vertex,face);
        
          normal(i,:) is the normal at vertex i.
          normalf(j,:) is the normal at face j.
        
          Copyright (c) 2004 Gabriel Peyre
    """

    eps = np.finfo(float).eps
    nface = np.shape(face)[1]
    nvert = np.shape(vertex)[1]
    normal = np.zeros([3,nvert])
    
    # unit normals to the faces
    normalf = crossp(vertex[:,face[1,:]] - vertex[:,face[0,:]],
                     vertex[:,face[2,:]] - vertex[:,face[0,:]])
    d = np.sqrt(np.sum(normalf**2,0))
    d[d < eps] = 1
    normalf = normalf/np.tile(d,(3,1))
    
    # unit normal to the vertex
    for i in range(nface):
        f = face[:,i]
        for j in range(3):
            normal[:,f[j]] = normal[:,f[j]] + normalf[:,i]

    # normalize
    d = np.sqrt(np.sum(normal**2,0))
    d[d < eps] = 1
    normal = normal/np.tile(d,(3,1))

    # enforce that the normal are outward
    v = vertex - np.tile(np.mean(vertex,0),(3,1))
    s = np.sum(v*normal,1)
    if np.sum(s > 0) < np.sum(s < 0):
        # flip
        normal = -normal
        normalf = -normalf
        
    return normal

##############################################################
def crossp(x,y):
        
    # x and y are (m,3) dimensional
    z = np.copy(x)
    z[0,:] = x[1,:]*y[2,:] - x[2,:]*y[1,:]
    z[1,:] = x[2,:]*y[0,:] - x[0,:]*y[2,:]
    z[2,:] = x[0,:]*y[1,:] - x[1,:]*y[0,:]
     
    return z
