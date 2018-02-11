import numpy as np
import matplotlib.pyplot as plt
from nt_toolbox.signal import imageplot

def plot_levelset(Z, level=0, f=[]):
    """
        f is supposed to be of the same shape as Z
    """
    if len(f) == 0:
        f = np.copy(Z)
        
    n,p = np.shape(Z)
    X,Y = np.meshgrid(np.arange(0,n),np.arange(0,p))
    plt.contour(X, Y, Z,[level],linewidths=2, colors="red")
    imageplot(f)
