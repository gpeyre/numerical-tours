import numpy as np
import matplotlib.pyplot as plt
from numpy import random

def plot_vf(velocities):
    """
        velocities is supposed to be of shape nxnx2
    """
    n = np.shape(velocities)[0]
    u = velocities[:,:,0]
    v = velocities[:,:,1]
    x,y = np.meshgrid(np.arange(n), np.arange(n))
    plt.quiver(x,y,u,v,color="b")
    plt.xlim(0,n)
    plt.ylim(0,n)
    plt.axis("off")
    plt.show()