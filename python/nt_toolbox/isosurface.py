import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D
    
def isosurface(M,v,step,title=""):
    """
    returns the isosurface of value v of M, subsetting M with the steps argument
    """
    sel = np.arange(0,np.shape(M)[0],step)
    
    verts, faces = measure.marching_cubes(M[np.ix_(sel,sel,sel)], v, spacing=(1.0, 1.0, 1.0))
    
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2], lw=.1, cmap="jet")
    ax.axis("off")
    ax.view_init(elev=35, azim=70)
    plt.title(title)
    plt.show()
