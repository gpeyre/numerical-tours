# import numpy as np
# import matplotlib.pyplot as plt
# from skimage import measure
# from mpl_toolkits.mplot3d import Axes3D

using PyCall
@pyimport skimage.measure as ms
@pyimport numpy as np


function isosurface(M, v, step, title)
    """
    returns the isosurface of value v of M, subsetting M with the steps argument
    """
    sel = 1 : step : size(M)[1]

    ix = np.ix_(sel, sel, sel)
    ix1 = ix[1][:]
    ix2 = ix[2][:]
    ix3 = ix[3][:]

    verts, faces = ms.marching_cubes(M[ix1, ix2, ix3], v, spacing = (1.0, 1.0, 1.0))

    fig = figure(figsize = (10,7))
    ax = fig[:add_subplot](111, projection = "3d")
    ax[:plot_trisurf](verts[:, 1], verts[:, 2], faces, verts[:, 3], lw = .1, cmap = "jet")
    ax[:axis]("off")
    ax[:view_init](elev = 35, azim = 70)
    # title(title)
    show()
end
