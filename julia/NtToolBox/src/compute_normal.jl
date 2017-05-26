function compute_normal(vertex,face)
    """
        compute_normal - compute the normal of a triangulation

          [normal,normalf] = compute_normal(vertex,face);

          normal(i,:) is the normal at vertex i.
          normalf(j,:) is the normal at face j.

          Copyright (c) 2004 Gabriel Peyre
    """
    epsilon = eps(Float64)
    nface = size(face,2)
    nvert = size(vertex,2)
    normal = zeros(3,nvert)

    # unit normals to the faces
    normalf = crossp(vertex[:,face[2,:]] - vertex[:,face[1,:]],
                     vertex[:,face[3,:]] - vertex[:,face[1,:]])
    d = sqrt(sum(normalf.^2,1))
    d[d .< epsilon] = 1
    normalf = normalf./repeat(d,outer=(3,1))

# unit normal to the vertex
    for i in 1:nface
        f = face[:,i]
        for j in 1:3
            normal[:,f[j]] = normal[:,f[j]] + normalf[:,i]
        end
    end

    # normalize
    d = sqrt(sum(normal.^2,1))
    d[d .< epsilon] = 1
    normal = normal./repeat(d,outer=(3,1))

    # enforce that the normal are outward
    v = vertex - repeat(mean(vertex,1),outer=(3,1))
    s = sum(v.*normal,2)
    if sum(s .> 0) < sum(s .< 0)
        # flip
        normal = -normal
        normalf = -normalf
    end
    return normal
end

###################################################
###################################################
###################################################

function crossp(x,y)

    # x and y are (m,3) dimensional
    z = copy(x)
    z[1,:] = x[2,:].*y[3,:] - x[3,:].*y[2,:]
    z[2,:] = x[3,:].*y[1,:] - x[1,:].*y[3,:]
    z[3,:] = x[1,:].*y[2,:] - x[2,:].*y[1,:]

    return z
end
