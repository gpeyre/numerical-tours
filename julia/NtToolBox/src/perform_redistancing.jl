
function perform_redistancing(D)
    """
        perform_redistancing - redistance a function

          D1 = perform_redistancing(D, options);

          Compute a signed distance function D1 that has the same 0 level set as D (2D-matrix).

          Copyright (c) 2007 Gabriel Peyre
    """

    n = size(D,1)

    # % horizontal
    P1 = D[1 : end - 1, :]
    P2 = D[2 : end, :]
    P = (P1.*P2) .<= 0
    d = abs(P1-P2)

    l = collect((i, j) for i in 1:size(d)[1] for j in 1:size(d)[2] if d[i, j] < eps())
    for i in 1:size(d)[1]
        for j in 1:size(d)[2]
            if (i, j) in l
                d[i, j] =1
            end
        end
    end

    v1 = abs(P1)./d
    v2 = abs(P2)./d
    Ah = ([P; zeros(1,n)] + [zeros(1,n); P]) .> 0
    Vh = max([v1; zeros(1, n)], [zeros(1, n); v2])
    # % vertical
    P1 = D[:, 1 : end - 1]
    P2 = D[:, 2 : end]
    P = (P1.*P2) .<= 0;
    d = abs(P1-P2)
    l = collect((i, j) for i in 1:size(d)[1] for j in 1:size(d)[2] if d[i, j] < eps())
    for i in 1:size(d)[1]
        for j in 1:size(d)[2]
            if (i, j) in l
                d[i, j] =1
            end
        end
    end
    v1 = abs(P1)./d
    v2 = abs(P2)./d
    Av = ([P zeros(n, 1)] + [zeros(n, 1) P]) .> 0
    Vv = max([v1 zeros(n, 1)], [zeros(n,1) v2])

    V = zeros(n, n)
    I = find(Ah .> 0)
    V[I] = Vh[I]
    I = find(Av .> 0)
    V[I] = max(V[I], Vv[I])

    I = find(V .!= 0)
    x, y = ind2sub(size(D), I)
    start_points = [x[:]'; y[:]']
    start_points[1, :], start_points[2, :] = start_points[2, :], start_points[1, :]

    D1 = NtToolBox.perform_fast_marching(ones(n, n), start_points)
    D1 = D1.*n
    D1[D .< 0] = -D1[D .< 0]

    return D1
end
