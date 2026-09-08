"""Numerical invariants for shared routines used throughout the tours."""

from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from nt_toolbox.general import crop, circshift
from nt_toolbox.signal import grad, div, bilinear_interpolate, perform_wavortho_transf
from nt_toolbox.compute_wavelet_filter import compute_wavelet_filter
from nt_toolbox.perform_wavelet_transf import perform_wavelet_transf
from nt_toolbox.perform_stft import perform_stft
from nt_toolbox.perform_linprog import perform_linprog
from nt_toolbox.perform_blurring import perform_convolution


def test_gradient_divergence_adjoint():
    rng = np.random.default_rng(4)
    x = rng.normal(size=(17, 13))
    v = rng.normal(size=(*x.shape, 2))
    np.testing.assert_allclose(np.sum(grad(x) * v), -np.sum(x * div(v)), atol=1e-12)


def test_bilinear_affine_and_boundary():
    y, x = np.mgrid[:5, :7]
    image = 2 * x + 3 * y + 1
    np.testing.assert_allclose(
        bilinear_interpolate(image, [1.5, 6, 20], [2.5, 4, 20]), [11.5, 25, 25]
    )


def test_crop_and_independent_shifts():
    x = np.arange(36).reshape(6, 6)
    np.testing.assert_array_equal(crop(x, 6), x)
    assert crop(x).shape == (3, 3)
    np.testing.assert_array_equal(
        circshift(x, [1, 2]), np.roll(x, (-1, -2), axis=(0, 1))
    )


def test_wavelet_perfect_reconstruction():
    x = np.random.default_rng(1).normal(size=(32, 32))
    h = compute_wavelet_filter("Daubechies", 4)
    w = perform_wavortho_transf(x, 0, 1, h)
    np.testing.assert_allclose(np.linalg.norm(w), np.linalg.norm(x), rtol=1e-12)
    np.testing.assert_allclose(perform_wavortho_transf(w, 0, -1, h), x, atol=1e-12)
    w = perform_wavelet_transf(x, 0, 1)
    np.testing.assert_allclose(perform_wavelet_transf(w, 0, -1), x, atol=1e-10)


def test_stft_energy_and_reconstruction():
    x = np.random.default_rng(2).normal(size=1024)
    spectrum = perform_stft(x, 64, 16, len(x))
    assert np.iscomplexobj(spectrum)
    np.testing.assert_allclose(np.linalg.norm(spectrum), np.linalg.norm(x), rtol=1e-12)
    np.testing.assert_allclose(perform_stft(spectrum, 64, 16, len(x)), x, atol=1e-12)


def test_linear_programming_feasibility_and_failure():
    A = np.array([[1.0, 1.0]])
    x = perform_linprog(A, [1.0], [1.0, 2.0])
    np.testing.assert_allclose(A @ x, [1.0])
    np.testing.assert_allclose(x, [1.0, 0.0])
    with pytest.raises(RuntimeError, match="failed"):
        perform_linprog(A, [-1.0], [1.0, 2.0])


@pytest.mark.parametrize("boundary", ["per", "sym"])
def test_convolution_preserves_constants_and_input(boundary):
    x = np.ones((9, 11, 3))
    h = np.ones((3, 3)) / 9
    y = perform_convolution(x, h, boundary)
    np.testing.assert_allclose(y, 1.0)
    np.testing.assert_array_equal(x, np.ones_like(x))
