function plot_vf(velocities)
    """
        velocities is supposed to be of shape nxnx2
    """
    n = size(velocities,1)
    u = velocities[:,:,1]
    v = velocities[:,:,2]
    x,y = meshgrid(0:n-1, 0:n-1)
    quiver(x,y,u,v,color="b")
    xlim(0,n)
    ylim(0,n)
    axis("off")
    show()
end
