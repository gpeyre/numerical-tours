function Div(Px, Py, bound = "sym", order = 1)
    """
        div - divergence operator

        fd = div(Px,Py, options);
        fd = div(P, options);

          options.bound = 'per' or 'sym'
          options.order = 1 (backward differences)
                        = 2 (centered differences)

          Note that the -div and grad operator are adjoint
          of each other such that
              <grad(f),g>=<f,-div(g)>

          See also: grad.

        Copyright (c) 2007 Gabriel Peyre
    """

    # retrieve number of dimensions
    nbdims = ndims(Px)

    if nbdims >= 3
        if nbdims == 3
            Py = P[:, :, 2]
            Px = P[:, :, 1]
            nbdims = 2
        else
            Pz = P[:, :, :, 3]
            Py = P[:, :, :, 2]
            Px = P[:, :, :, 1]
            nbdims = 3
        end
    end

    if bound == "sym"
        nx = size(Px)[1]
        if order == 1
            fx = Px .- Px[vcat(1, collect(1, nx - 1)), :]
            fx[1, :] = Px[1, :]                        # boundary
            fx[nx, :] = - Px[nx - 1, :]

            if nbdims >= 2
                ny = size(Py)[2]
                fy = Py - Py[:, vcat(1, collect(1, ny - 1))]
                fy[:, 1] = Py[:, 1]                    # boundary
                fy[:, ny] = - Py[:, ny - 1]
            end

            if nbdims >= 3
                nz = size(Pz)[3]
                fz = Pz - Pz[:, :, vcat(1, collect(1 : nz - 1))]
                fz[:, :, 1] = Pz[:, :, 1]                # boundary
                fz[:, :, nz] = - Pz[:, :, nz - 1]
            end
        else
            fx = (Px[vcat(collect(2 : nx), nx), :] - Px[vcat([1], collect(1, nx-1)), :])./2.
            fx[1, :] = + Px[1, :]./2. + Px[1, :]           # boundary
            fx[1, :] = + Px[2, :]./2. - Px[1, :]
            fx[nx, :] = - Px[nx, :] - Px[nx - 1, :]./2.
            fx[nx - 1, :] = + Px[nx, :] - Px[nx - 2, :]./2.

            if nbdims >= 2
                ny = size(Py)[2]
                fy = (Py[:, vcat(collect(2 : ny), ny)] - Py[:, vcat(1, collect(1 : ny-1))])./2.
                fy[:, 1] = + Py[:, 2]./2. + Py[:, 1]       # boundary
                fy[:, 2] = + Py[:, 3]./2. - Py[:,1]
                fy[:, ny] = - Py[:, ny] - Py[:, ny-1]./2.
                fy[:, ny - 1] = + Py[:, ny] - Py[:, ny - 2]./2.
            end

            if nbdims >= 3
                nz = size(Pz)[3]
                fz = (Pz[:, :, vcat(collect(2 : nz), nz)] - Pz[:, :, vcat(1, collect(1 : nz - 1))])./2.
                fz[:, :, 1] = + Pz[:, :, 2]./2. + Pz[:, :, 1] # boundary
                fz[:, :, 2] = + Pz[:, :, 3]./2. - Pz[:, :, 1]
                fz[:, :, ny] = - Pz[:, :, nz] - Pz[:, :, nz - 1]./2.
                fz[:, :, ny - 1] = + Pz[:, :, nz] - Pz[:, :, nz - 2]./2.
            end
        end

    else
        if order == 1
            nx = size(Px)[1]
            fx = Px - Px[vcat(nx, collect(1 : nx-1)), :]

            if nbdims >= 2
                ny = size(Py)[2]
                fy = Py - Py[:, vcat(ny, collect(1 : ny-1))]
            end

            if nbdims>=3
                nz = size(Pz)[3]
                fz = Pz - Pz[:, :, vcat(nz, collect(1 : nz-1))]
            end

        else
            nx = size(Px)[1]
            fx = (Px[vcat(collect(2 : nx), 1), :]) - (Px[vcat(nx, collect(1, nx - 1)), :])

            if nbdims >= 2
                ny = size(Py)[2]
                fy = (Py[:, vcat(collect(2 : ny), 1)]) - (Py[:, vcat(ny, collect(1, ny - 1))])
            end

            if nbdims >= 3
                nz = size(Pz)[3]
                fz = (Pz[:, :, vcat(collect(2, nz), 1)]) - (Pz[:, :, vcat(nz, collect(1 : nz - 1))])
            end
        end
    end

    # gather result
    if nbdims == 3
        fd = fx + fy + fz

    elseif nbdims == 2
        fd = fx + fy

    else
        fd = fx
    end

    return fd
end
