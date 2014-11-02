def exo1():
    """
    Compare the log of the modulus of the Fourier transforms of the input image $f$
    and its periodic component $p$. What do you observe ?
    """
    mydisp = lambda x: log10(abs(fftshift(x)) + 1e-5)
    imageplot(mydisp(fft2(f)), 'Original', 1, 2, 1)
    imageplot(mydisp(fft2(p)), 'Periodic layer', 1, 2, 2)


def exo2():
    """
    Compare the histograms of the input and synthesized textures pixel empirical distributions.
    """
    k = 100
    t = linspace(0, .5, k)
    subplot(2, 1, 1)
    bar(t, histc(f0(: ), t)*k/ n^2); axis tight
    title('Input')
    subplot(2, 1, 2)
    bar(t, histc(f(: ), t)*k/ n^2);  axis tight
    title('Synthesized')


def exo3():
    """
    Perform the color texture synthesis.
    """
    f = real(ifft2(fft2(f0).*repmat(fft2(w), [1 1 size(f0, 3)])))
    u = f0; u(1, 1, : ) = 0; u(2, 1, : ) = 1
    imageplot(clamp(u), 'Input', 1, 2, 1)
    u = f; u(1, 1, : ) = 0; u(2, 1, : ) = 1
    imageplot(clamp(u), 'Synthesized', 1, 2, 2)


def exo4():
    """
    Compute several realizations of the color texture synthesis.
    """
    for i in 1: 4:
        w = randn(n)/ n;  w = w-mean(w(: )) + 1/ n^2
        f = real(ifft2(fft2(f0).*repmat(fft2(w), [1 1 size(f0, 3)])))
        u = f; u(1, 1, : ) = 0; u(2, 1, : ) = 1
        subplot(2, 2, i)
        imageplot(clamp(u))


def exo5():
    """
    Perform the color texture synthesis using this extended spot noise.
    """
    w = randn(n0)/ n0;  w = w-mean(w(: )) + 1/ n0^2
    f = real(ifft2(fft2(f1).*repmat(fft2(w), [1 1 size(f0, 3)])))
    u = ones(n0, n0, 3)
    u(end/ 2-n/ 2 + 1: end/ 2 + n/ 2, end/ 2-n/ 2 + 1: end/ 2 + n/ 2, : ) = f0
    u(1, 1, : ) = 0; u(2, 1, : ) = 1
    imageplot(clamp(u), 'Input', 1, 2, 1)
    u = f; u(1, 1, : ) = 0; u(2, 1, : ) = 1
    imageplot(clamp(u), 'Synthesized', 1, 2, 2)


