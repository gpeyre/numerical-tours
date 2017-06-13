function plot_mesh(X, F; sub = [1,1,1], titre="", el=90, az=-90, lwdt=.1, dist=6, c="grey")
    """
    plots the mesh defined by X : vertices and F : faces
    """
    subplot(sub[1],sub[2],sub[3], projection="3d")
    plot_trisurf(X[1,:], X[2,:], X[3,:], triangles=F'-1, lw=lwdt, color=c, alpha=1)
    axis("off")
    gca()[:view_init](el, az)
    gca()[:dist]=dist
    title(titre)
end
