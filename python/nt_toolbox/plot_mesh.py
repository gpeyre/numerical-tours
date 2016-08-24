import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D
    
def plot_mesh(X, F, subplot = [1,1,1], title="", el=90, az=-90, lwdt=.1, dist=6, c="grey"):
    """
    plots the mesh defined by X : vertices and F : faces
    """
    ax = plt.subplot(subplot[0],subplot[1],subplot[2], projection='3d')
    ax.plot_trisurf(X[0,:], X[1,:], X[2,:], triangles=np.transpose(F), lw=lwdt, color=c, alpha=1)
    ax.axis("off")
    ax.view_init(elev=el, azim=az)
    ax.dist = dist
    plt.title(title)
