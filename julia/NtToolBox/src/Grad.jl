function Grad(M, bound = "sym", order = 1)
    """
        grad - gradient, forward differences

          [gx,gy] = grad(M, options);
        or
          g = grad(M, options);

          options.bound = 'per' or 'sym'
          options.order = 1 (backward differences)
                        = 2 (centered differences)

          Works also for 3D array.
          Assme that the function is evenly sampled with sampling step 1.

          See also: div.

          Copyright (c) Gabriel Peyre
    """


    # retrieve number of dimensions
    nbdims = ndims(M)


    if bound == "sym"
        nx = size(M)[1]
        if order == 1
            fx = M[vcat(collect(2 : nx), nx), :] - M
        else
            fx = (M[vcat(collect(2 : nx), nx), :] - M[vcat(1, collect(1 : nx - 1)), :])./2.
            # boundary
            fx[1, :] = M[2, :] - M[1, :]
            fx[nx, :] = M[nx, :] - M[nx - 1, :]
        end

        if nbdims >= 2
            ny = size(M)[2]
            if order == 1
                fy = M[:, vcat(collect(2 : ny), ny)] - M
            else
                fy = (M[:, vcat(collect(2 : ny), ny)] - M[:, vcat(1, collect(1 : ny - 1))])./2.
                # boundary
                fy[:, 1] = M[:, 2] - M[:, 1]
                fy[:, ny] = M[:, ny] - M[:, ny - 1]
            end
        end

        if nbdims >= 3
            nz = size(M)[3]
            if order == 1
                fz = M[:, :, vcat(collect(2 : nz), nz)] - M
            else
                fz = (M[:, :, vcat(collect(2 : nz), nz)] - M[:, :, vcat(1, collect(1 : nz - 1))])./2.
                # boundary
                fz[:, :, 1] = M[:, :, 2] - M[:, :, 1]
                fz[:, :, ny] = M[:, :, nz] - M[:, :, nz - 1]
            end
        end
    else
        nx = size(M)[1]
        if order == 1
            fx = M[vcat(collect(2 : nx), 1), :] - M
        else
            fx = (M[vcat(collect(2 : nx), 1), :] - M[vcat(nx, collect(1 : nx - 1)), :])./2.
        end

        if nbdims >= 2
            ny = size(M)[2]
            if order == 1
                fy = M[:, vcat(collect(2 : ny), 1)] - M
            else
                fy = (M[:, vcat(collect(2 : ny), 1)] - M[:, vcat(ny, collect(1 : ny - 1))])./2.
            end
        end

        if nbdims >= 3
            nz = size(M)[3]
            if order == 1
                fz = M[:, :, vcat(collect(2 : nz), 1)] - M
            else
                fz = (M[:, :, vcat(collect(2 : nz), 1)] - M[:, :, vcat(nz, collect(1 : nz - 1))])./2.
            end
        end
    end

    if nbdims==2
        fx = cat(3, fx[:, :], fy[:, :])
    elseif nbdims==3
        fx = cat(4, (fx[:, :, :], fy[:, :, :], fz[:, :, :]))
    end

    fx[end, :, 1] = fx[end - 1, :, 1]
    fx[:, end, 2] = fx[:, end - 1, 2]


    return fx
end