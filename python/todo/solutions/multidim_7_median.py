def exo1():
    """
    A first way to denoise the image is to apply the local median filter
    implemented with the function |perform_median_filtering| on
    each channel |M(:,:,i)| of the image, to get a denoised image |Mindep| with SNR |pindep|.
    """
    Mindep = zeros(n, n, 3)
    for i in 1: 3:
        Mindep(: , : , i) = perform_median_filtering(M(: , : , i), k)
    pindep = snr(M0, Mindep)
    imageplot(clamp(M), strcat(['Noisy, SNR = ' num2str(pnoisy)]), 1, 2, 1)
    imageplot(clamp(Mindep), strcat(['Denoised, SNR = ' num2str(pindep)]), 1, 2, 2)


def exo2():
    """
    Compute the median |med| of the points in |X|
    using the iterative reweighted least squares algorithm.
    This computed median |med| should be stored in the result as
    |Mmed(x,y,:)| (you need to reshape |med| so that its size is |[1 1 3]|).
    """
    med = mean(X, 2)
    niter = 8; energy = []
    for i in 1: niter:
        % comute the distance from med to the points
        dist = sqrt(sum((X-repmat(med, [1 w*w])).^2))
        % compute the weight, take care of not dividing by 0
        weight = 1./ max(dist, 1e-10); weight = weight/ sum(weight)
        % compute the weighted mean
        med = sum(repmat(weight, [3 1]).*X, 2)
        energy(end + 1) = sum(dist)


def exo3():
    """
    Implement the 3D median filter by looping through all the pixel |(x,y)|.
    isplay the results
    """
    Mmed = zeros(n, n, 3)
    niter = 15; t = 0
    for x in 1: n:
        selx = x-k: x + k; selx = mod(selx-1, n) + 1
    for y in 1: n:
            t = t + 1; % progressbar(t, n*n)
            % extract patch
            sely = y-k: y + k; sely = mod(sely-1, n) + 1
            X = reshape(M(selx, sely, : ), [w*w 3])'
            % compute median
            med = mean(X, 2)
    for i in 1: niter:
                dist = sqrt(sum((X-repmat(med, [1 w*w])).^2))
                weight = 1./ max(dist, 1e-10); weight = weight/ sum(weight)
                med = sum(repmat(weight, [3 1]).*X, 2)
            % store result
            Mmed(x, y, : ) = reshape(med, [1 1 3])
    pmed = snr(M0, Mmed)
    imageplot(clamp(Mindep), strcat(['1D median, SNR = ' num2str(pindep)]), 1, 2, 1)
    imageplot(clamp(Mmed), strcat(['3D median, SNR = ' num2str(pmed)]), 1, 2, 2)


