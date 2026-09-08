import numpy as np
import pylab as pyl
from scipy import signal


def perform_blurring(M, sigma, bound="sym"):
    """
    perform_blurring - gaussian blurs an image

    M = perform_blurring(M, sigma, options);

    M is the original data
    sigma is the width of blurs (in pixels)

    Copyright (c) 2007 Gabriel Peyre
    """

    if np.all(sigma == 0):
        return M

    if np.ndim(M) > 2:
        for i in range(np.shape(M)[2]):
            M[:, :, i] = perform_blurring(M[:, :, i], sigma, bound)

    n = max(np.shape(M))

    eta = 4
    p = np.round((sigma * eta) / 2.0) * 2 + 1
    p = np.minimum(p, (round(n // 2.0) * 2 - 1) * np.ones(len(p)))

    A = np.array([1.0, 1.0])
    if np.ndim(M) == 1:
        A = 1  # 1D

    h = compute_gaussian_filter(p * A, sigma / (4.0 * n), n * A)
    M = perform_convolution(M, h, bound)

    return M


def compute_gaussian_filter(n, s, N):
    """
    compute_gaussian_filter - compute a 1D or 2D Gaussian filter.

      f = compute_gaussian_filter(n,s,N);

      'n' is the size of the filter, odd for no phase in the filter.
          (if too small it will alterate the filter).
          use n=[n1,n2] for a 2D filter or n = [n1] for a 1D filter
      's' is the standard deviation of the filter.
      'N' is the size of the big signal/image (supposed to lie in [0,1] or [0,1]x[0,1]).
          use N=[N1,N2] for a 2D filter or N = [N1] for a 1D filter

      The equation (in 1D) is
          f[k] = exp( -(x(k)^2/(2*s^2)) );
      where x spans [-1/2,1/2].

      The filter is normalised so that it sums to 1.

      Copyright (c) 2004 Gabriel Peyre
    """
    nd = 1
    if len(n) > 1 and n[1] > 1:
        nd = 2

    if nd == 2 and len(s) == 1:
        s = np.hstack((s, s))

    if nd == 2 and len(N) == 1:
        N = np.hstack((N, N))

    if nd == 1:
        f = build_gaussian_filter_1d(n, s, N)
    else:
        f = build_gaussian_filter_2d(n, s, N)
    return f


def build_gaussian_filter_2d(n, s, N=[]):
    """
    build_gaussian_filter_2d - compute a 2D Gaussian filter.

    f = build_gaussian_filter_2d(n,s,N);

    'n' is the size of the filter, odd for no phase in the filter.
        (if too small it will alterate the filter).
    's' is the standard deviation of the filter.
    'N' is the size of the big image (supposed to lie in [0,1]x[0,1]).

    The filter is normalised so that it sums to 1.

    Copyright (c) 2004 Gabriel Peyre
    """

    n = np.asarray(n)
    s = np.asarray(s)
    N = np.asarray(N)

    if len(N) == 0:
        N = n

    if len(N) == 1 or N[0] == 1:
        N = np.hstack((N, N))

    if len(s) == 1 or s[0] == 1:
        s = np.hstack((s, s))

    if len(s[s <= 0]) > 0:
        f = np.zeros(n)
        f[np.round((n - 1) / 2).astype(int)] = 1
        return f

    x = (np.arange(0, n[0]) - (n[0] - 1) / 2.0) / (N[0] - 1)
    y = (np.arange(0, n[1]) - (n[1] - 1) / 2.0) / (N[1] - 1)
    [Y, X] = np.meshgrid(y, x)
    f = np.exp(-(X**2 / (2 * s[0] ** 2)) - (Y**2 / (2 * s[1] ** 2)))
    f = f / np.sum(f)
    return f


def build_gaussian_filter_1d(n, s, N=[]):
    """
    build_gaussian_filter_1d - compute a Gaussian filter.

    f = build_gaussian_filter_1d(n,s,N);

    Copyright (c) 2004 Gabriel Peyre
    """
    if len(N) == 0:
        N = n

    n = n[0]
    s = s[0]
    N = N[0]

    if s <= 0:
        f = np.zeros(n)
        f[np.round((n - 1) / 2)] = 1
        return f

    x = (np.arange(0, n) - (n - 1) / 2.0) / (N - 1)
    f = np.exp(-(x**2) / (2 * s**2))
    f = f / np.sum(f)
    return f


def perform_convolution(x, h, bound="sym"):
    """Convolve with a centered kernel and symmetric or periodic boundaries.

    Color channels are processed independently. The input is never modified.
    """
    from scipy.ndimage import convolve

    x = np.asarray(x)
    h = np.asarray(h)
    if bound not in {"sym", "per"}:
        raise ValueError("bound must be 'sym' or 'per'")
    mode = "reflect" if bound == "sym" else "wrap"
    if x.ndim == 3 and h.ndim == 2:
        return np.stack(
            [convolve(x[..., k], h, mode=mode) for k in range(x.shape[-1])], axis=-1
        )
    if x.ndim != h.ndim:
        raise ValueError("The input and kernel must have compatible dimensions")
    return convolve(x, h, mode=mode)
