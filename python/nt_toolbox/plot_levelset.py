import numpy as np
import matplotlib.pyplot as plt
from nt_toolbox.signal import imageplot

def plot_levelset(Z,level,f):
    """
	f is supposed to be of the shape as Z
    """
    n,p = np.shape(Z)
    X,Y = np.meshgrid(np.arange(0,n),np.arange(0,p))
    plt.contour(X, Y, Z,[level],linewidths=2, colors="red")
    imageplot(f)