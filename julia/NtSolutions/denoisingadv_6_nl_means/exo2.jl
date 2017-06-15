figure(figsize = (10, 10))

tau = .03
q_list = [10, 20]
w_list = [3, 6]
ind_plot = 0
for i_q in 1 : length(q_list)
    for i_w in 1 : length(w_list)

        w = w_list[i_w]
        q = q_list[i_q]
        ind_plot += 1

        #patch
        w1 = 2*w + 1

        (Y, X, dX, dY) = ndgrid(1 : n, 1 : n, -w : w, -w : w)
        X = X + dX
        Y = Y + dY

        X[X .< 1] = 2 .- X[X .< 1]
        Y[Y .< 1] = 2 .- Y[Y .< 1]
        X[X .> n] = 2*n .- X[X .> n]
        Y[Y .> n] = 2*n .- Y[Y .> n]

        I = X .+ (Y .- 1)*n
        for k in 1 : Base.div(n, w)
            for l in 1 : Base.div(n, w)
                I[k, l, :, :] = transpose(I[k, l, :, :])
            end
        end

        function patch(f) # There is no problem with patch function
            return f'[I]
        end

        P = patch(test)

        #PCA
        resh = P -> transpose((reshape(P, (n*n,w1*w1))))
        remove_mean = Q -> Q - repeat(mean(Q, 1), inner = (w1*w1, 1))

        P1 = remove_mean(resh(P))
        C = P1*transpose(P1)
        (D, V) = eig(C)
        D = D[end : -1 : 1]
        V = V[:, end : -1 : 1]

        iresh = Q -> reshape(Q', (n, n, d))
        descriptor = f -> iresh(V[: , 1 : d]'*remove_mean(resh(P)))

        H = descriptor(f)

        #NL_means
        i = [84, 73]

        function distance_0(i,sel)
            H1 = (H[sel[1, :], :, :])
            H2 = (H1[:, sel[2, :], :])
            return sum((H2 - repeat(reshape(H[i[1], i[2], :], (1, 1, length(H[i[1], i[2], :]))), inner = [length(sel[1, :]), length(sel[1, :]), 1])).^2, 3)/w1*w1
        end

        distance = i -> distance_0(i, selection(i))
        kernel = (i, tau) -> normalize(exp(-distance(i)./(2*tau^2)))
        selection = i -> [clamP(i[1] - q : i[1] + q, 1, n)'; clamP(i[2] - q : i[2] + q, 1, n)']

        function NLval_0(K,sel)
        #     f_temp = f[sel[1, :], :]
            f_temp = test[sel[1, :], :]
            return sum(K.*f_temp[:, sel[2, :]])
        end

        NLval = (i, ta) -> NLval_0(kernel(i, tau), selection(i))
        NLmeans = tau -> arrayfun((i1, i2) -> NLval([i1,i2], tau), X, Y)
        f1 = NLmeans(tau)

        imageplot(clamP(f1), @sprintf("q = %i, w = %i, SNR = %.1f dB" , q, w, snr(f0,f1)), [2, 2, ind_plot])
    end
end
