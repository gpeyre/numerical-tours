# import numpy as np
# import matplotlib.pyplot as plt
# from skimage import measure
# from mpl_toolkits.mplot3d import Axes3D

using Meshing

function isosurface(M, v, step, title = "")
    """
    returns the isosurface of value v of M, subsetting M with the steps argument
    """
    sel = 1 : step : size(M)[1]

    (verts, faces) = marching_cubes(M[sel, sel, sel], v, (1.0, 1.0, 1.0))

    fig = figure(figsize = (10,7))
    ax = add_subplot(111, projection = "3d")
    plot_trisurf(verts[:, 1], verts[:, 2], faces, verts[:, 3], lw = .1, cmap = "jet")
    axis("off")
    view_init(elev = 35, azim = 70)
    title(title)
    show()
end
