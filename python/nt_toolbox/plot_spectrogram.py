import matplotlib.pyplot as plt
import numpy as np

def plot_spectrogram(S,title):

    S = abs(S[:np.shape(S)[0]/2,:])
    S = np.log(S + 1e-4)
    plt.imshow(S,cmap = plt.get_cmap("jet"), interpolation = "nearest")
    plt.title(title)
    plt.show()