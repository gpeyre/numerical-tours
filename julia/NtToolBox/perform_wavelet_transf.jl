function perform_wavelet_transf(f, Jmin, dir, filter = "9-7",separable = 0, ti = 0)
  """""
  perform_wavelet_transf - peform fast lifting transform

  y = perform_wavelet_transf(x, Jmin, dir, filter = "9-7",separable = 0, ti = 0);

  Implement 1D and 2D symmetric wavelets with symmetric boundary treatements, using
  a lifting implementation.

  filter gives the coefficients of the lifting filter.
  You can use h='linear' or h='7-9' to select automatically biorthogonal
  transform with 2 and 4 vanishing moments.

  You can set ti=1 to compute a translation invariant wavelet transform.

  You can set separable=1 to compute a separable 2D wavelet
  transform.

  Copyright (c) 2008 Gabriel Peyre
  """

  # Copy f.
  x = copy(f)

  # Convert Jmin to Int.
  Jmin = Int(Jmin)

  # Number of dimensions of x.
  d = ndims(x)

  if filter in ["linear","5-3"]
    h = [1/2, 1/4, sqrt(2)]

  elseif filter in ["9-7","7-9"]
    h = [1.586134342, -.05298011854, -.8829110762, .4435068522, 1.149604398]

  else
    error("Unknown filter")   # Raising an error if the input filter is unknown.

  end

  if (d == 2) & (separable == 1)
    ti = 0
    if ti == 1
      warn("Separable does not works for translation invariant transform")

      # perform a separable wavelet transform
    end

    n = size(x)[1]

    if dir == 1
      for i in 1:n
        x[:,i] = perform_wavelet_transf(x[:,i], Jmin, dir, filter, separable, ti)
      end
      for i in 1:n
        x[i,:] = transpose(perform_wavelet_transf(np.transpose(x[i,:]), Jmin, dir, filter, separable, ti))
      end

    else
      for i in 1:n
        x[i,:] = transpose(perform_wavelet_transf(np.transpose(x[i,:]), Jmin, dir, filter, separable, ti))
      end
      for i in 1:n
        x[:,i] = perform_wavelet_transf(x[:,i], Jmin, dir, filter, separable, ti)
      end
    end

    # number of lifting steps
    if ndims(x) == 1
      n = length(x)
    else
      n = size(x)[2]
    end
    m = Int((length(h) - 1)/2)
    Jmax = Int(log2(n) - 1)
    jlist = Jmax:-1:Jmin

    if dir == -1
      jlist = Jmin:Jmax
    end

    if ti == 0
      # subsampled
      for j in jlist
        if d == 1
          x[1:2^(j+1), :] = lifting_step(x[1:2^(j+1)], h, dir) # A vérifier les indices.
        else
          x[1:2^(j+1), 1:2^(j+1)] = lifting_step(x[1:2^(j+1), 1:2^(j+1)], h, dir) # A vérifier les indices.
          x[1:2^(j+1), 1:2^(j+1)] = transpose(lifting_step(transpose(x[1:2^(j+1), 1:2^(j+1)]), h, dir)) # A vérifier les indices.
        end
      end

    else
      # TI
      nJ = Jmax - Jmin + 1
      if (dir == 1) & (d == 1)
        x = repeat(x, outer = [nJ + 1, 1, 1])   # A vérifier.
      elseif (dir == 1) & (d == 2)
        x = repeat(x, outer = [3*nJ + 1,1,1])
      end
      for j in jlist
        dist = 2^(Jmax - j)

        if d == 1
          if dir == 1
            x[1:(j - Jmin + 2), :, :] = lifting_step_ti(x[0, :, :], h, dir, dist)
          else
            x[1, :, :] = lifting_step_ti(x[1:(j - Jmin + 2), :, :], h, dir, dist)
          end
        else
          dj = 3*(j - Jmin)

          if dir == 1
#              x[1 : dj + 1, :, :] = lifting_step_ti(x[1, :, :], h, dir, dist) # A vérifier avec précautions ligne 104
              x[[1, dj + 1], :, :] = lifting_step_ti(x[1, :, :], h, dir, dist)

#              x[1:dj + 2, :, :] = lifting_step_ti(transpose(x[1, :, :]), h, dir, dist) # A vérifier avec précautions ligne 106
              x[[1, dj + 2],:,:] = lifting_step_ti(transpose(x[1, :, :]), h, dir, dist)
              x[1, :, :] = transpose(x[1, :, :])
              x[dj + 3, :, :] = transpose(x[dj + 3, :, :])

#              x[2 + dj : 3 + dj, :, :] = lifting_step_ti(transpose(x[dj, :, :]), h, dir, dist)
              x[[2 + dj, 3 + dj], :, :] = lifting_step_ti(transpose(x[dj + 2, :, :]), h, dir, dist)
              x[dj + 2, :, :] = transpose(x[dj + 2, :, :])
              x[dj + 4, :, :] = transpose(x[dj + 4, :, :])
          else

              x[dj + 2, :, :] = transpose(x[dj + 2, :, :])
              x[dj + 4, :, :] = transpose(x[dj + 4, :, :])

#              x[dj + 2, :, :] = transpose(lifting_step_ti(x[[2 + dj : 3 + dj], :, :], h, dir, dist))
              x[dj + 2, :, :] = transpose(lifting_step_ti(x[[2 + dj, 3 + dj], :, :], h, dir, dist))

              x[1, :, :] = transpose(x[1, :, :])
              x[dj + 3, :, :] = transpose(x[dj + 3, :, :])
#              x[1, :, :] = transpose(lifting_step_ti(x[[1 : dj + 2], :, :], h, dir, dist))
              x[1, :, :] = transpose(lifting_step_ti(x[[1, dj+2],:,:], h, dir, dist))

#              x[1, :, :] = lifting_step_ti(x[[1, dj+1], :, :], h, dir, dist)
              x[1, :, :] = lifting_step_ti(x[[1 : dj+1], :, :], h, dir, dist)
            end
          end
        end
      end
      if dir == -1
        x = x[1, :, :]
      end
    end
    return x
end

function lifting_step(x0, h, dir)

    #copy x
    x = copy(x0)

    # number of lifting steps
    m = Int((length(h) - 1)/2)

    if dir == 1
        # split
#        d = x[1::2,]
        d = x[collect(i for i in 1:length(x) if i%2 == 0)]
#        x = x[0::2,]
        x = x[collect(i for i in 1:length(x) if i%2 == 1)]

        for i in 1 : m
            d = d - h[2*i] * (x + [x[1 : end, :]'; x[end, :]']) # I think we have to transpose.
            x = x + h[2*i + 1] * (d + [d[1, :]'; d[1 : end - 1,]']) #Be very cautiaus about transposing.
        end
        x = [x*h[end]; d/h[end]]

    else
        # retrieve detail coefs
        fin = length(x)
        d = x[fin/2 + 1 : end, :]*h[end]
        x = x[1:fin/2, :]/h[end]
        for i in m:-1:1
            x = x - h[2*i - 1] * (d + [d[1, :]; d[1:end - 1, :]])
            d = d + h[2*i - 2] * (x + [x[2:end, :]; x[end, :]])
        end
        # merge
        x1 = [x; x]
#        x1[::2,] = x
        x1[collect(i for i in 1:length(x[:, 1]) if i%2 == 1)] = x
#        x1[1::2,] = d
        x1[collect(i for i in 1:length(x[:, 1]) if i%2 == 0)] = d
        x = x1
    end

    return x
end

function lifting_step_ti(x0, h, dir, dist)

    #copy x
    x = copy(x0)

    # number of lifting steps
    m = Int((len(h) - 1)/2)
    n = size(x[0])[0]

    s1 = collect(i for i in 1:n) + dist
    s2 = collect(i for i in 1:n) - dist

    # boundary conditions
    s1[collect(i for i in 1:length(s1) if s1[i] > n)] = 2*n - s1[collect(i for i in 1:length(s1) if s1[i] > n)]
    s1[collect(i for i in 1:length(s1) if s1[i] < 1)] = 2 - s1[collect(i for i in 1:length(s1) if s1[i] < 1)]

    s2[collect(i for i in 1:length(s2) if s2[i] > n)] = 2*n - s2[collect(i for i in 1:length(s2) if s2[i] > n)]
    s2[s2 < 1] = 2   - s2[s2 < 1]
    s2[collect(i for i in 1:length(s2) if s2[i] < 1)] = 2 - s2[collect(i for i in 1:length(s2) if s2[i] < 1)]


    if dir == 1
        # split
        d = x
        for i in 1:m
            if ndims(x) == 2
                x = repeat(x, outer = [1, 1, 1])
            end
            d = d - h[2*i] * (x[:, s1, :] + x[:, s2, :])
            x = x + h[2*i + 1] * (d[:, s1, :] + d[:, s2, :])
        end

        #merge
#        x = np.concatenate((x*h[end],d/h[end]))
        x = [x*h[end]; d/h[end]]

    else
        # retrieve detail coefs

        d = x[1, :, :]*h[end]
        x = x[0, :, :]/h[end]

        for i in m:-1:1
            x = x - h[2*i - 1] * (d[s1, :] + d[s2, :])
            d = d + h[2*i - 2] * (x[s1, :] + x[s2, :])
        end

        # merge
        x = (x + d)/2
    end

    return x
end
