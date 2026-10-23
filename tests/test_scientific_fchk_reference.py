"""Absolute FCHK basis amplitudes from independent radial/angular formulas."""

import math
from pathlib import Path

import numpy as np
import pytest
from conftest import load_plugin, mocks_with_real_numpy

PATH = (
    Path(__file__).resolve().parents[1]
    / "plugins/Gaussian_MO_Analyzer/gaussian_fchk_mo_analyzer/analyzer.py"
)
with mocks_with_real_numpy():
    module = load_plugin(PATH)


def engine(shell_type, exponents, coefficients, p_coefficients=None):
    data = {
        "Shell types": np.array([shell_type]),
        "Number of primitives per shell": np.array([len(exponents)]),
        "Shell to atom map": np.array([1]),
        "Primitive exponents": np.array(exponents),
        "Contraction coefficients": np.array(coefficients),
        "Current cartesian coordinates": np.array([0.2, -0.3, 0.4]),
        "Pure/Cartesian d shells": np.array([0]),
        "Pure/Cartesian f shells": np.array([0]),
    }
    if p_coefficients is not None:
        data["P(S=P) Contraction coefficients"] = np.array(p_coefficients)
    return module.BasisSetEngine(data)


@pytest.mark.parametrize("degree", [2, 3, 4], ids=["5D", "7F", "9G"])
@pytest.mark.parametrize("alpha", [0.19, 1.0, 3.7])
def test_absolute_signed_spherical_amplitudes(degree, alpha):
    special = pytest.importorskip("scipy.special")
    basis = engine(-degree, [alpha], [1.0])
    # Different radii distinguish angular shape from radial normalization.
    points = np.array(
        [[0.3, 0.8, -0.5], [-0.9, 0.2, 0.7], [1.2, -0.4, 0.6], [-0.2, -0.3, -0.1]]
    )
    radius = np.linalg.norm(points, axis=1)
    theta = np.arccos(points[:, 2] / radius)
    azimuth = np.arctan2(points[:, 1], points[:, 0])
    radial_norm = math.sqrt(
        2 * (2 * alpha) ** (degree + 1.5) / math.gamma(degree + 1.5)
    )
    order = [0] + [m for n in range(1, degree + 1) for m in (n, -n)]
    for component, m in enumerate(order):
        y = (
            special.sph_harm_y(degree, abs(m), theta, azimuth)
            if hasattr(special, "sph_harm_y")
            else special.sph_harm(abs(m), degree, azimuth, theta)
        )
        angular = (
            y.real
            if m == 0
            else math.sqrt(2) * (-1) ** abs(m) * (y.real if m > 0 else y.imag)
        )
        expected = radial_norm * radius**degree * np.exp(-alpha * radius**2) * angular
        coefficients = np.eye(basis.n_basis)[component]
        got = basis.evaluate_mo_on_grid(0, points + [0.2, -0.3, 0.4], coefficients)
        np.testing.assert_allclose(got, expected, rtol=2e-12, atol=1e-13)


def test_sp_contraction_uses_separate_signed_s_and_p_channels():
    exponents = [0.23, 1.7, 4.1]
    cs, cp = [0.2, -0.5, 0.7], [-0.8, 0.4, 0.1]
    basis = engine(-1, exponents, cs, cp)
    delta = np.array([[0.7, -0.3, 0.2], [-0.2, 0.8, -0.4]])
    r2 = np.sum(delta**2, axis=1)
    expected = np.zeros((len(delta), 4))
    for alpha, s_coeff, p_coeff in zip(exponents, cs, cp):
        s = (2 * alpha / math.pi) ** 0.75 * np.exp(-alpha * r2)
        expected[:, 0] += s_coeff * s
        expected[:, 1:] += p_coeff * (2 * math.sqrt(alpha)) * s[:, None] * delta
    for component in range(4):
        actual = basis.evaluate_mo_on_grid(
            0, delta + [0.2, -0.3, 0.4], np.eye(4)[component]
        )
        np.testing.assert_allclose(actual, expected[:, component], atol=1e-13)


@pytest.mark.parametrize("degree", [2, 3, 4])
def test_mixed_cartesian_and_pure_shells_cannot_share_one_definition(degree):
    # Without Number of basis functions, the old global flag could make both
    # signed types consume the same number of coefficients and silently pass.
    data = {
        "Shell types": np.array([degree, -degree]),
        "Number of primitives per shell": np.array([1, 1]),
        "Shell to atom map": np.array([1, 1]),
        "Primitive exponents": np.array([1.0, 0.5]),
        "Contraction coefficients": np.array([1.0, 1.0]),
        "Current cartesian coordinates": np.array([0.0, 0.0, 0.0]),
        "Pure/Cartesian d shells": np.array([0]),
        "Pure/Cartesian f shells": np.array([0]),
        "Alpha MO coefficients": np.zeros(2 * (2 * degree + 1)),
    }
    with pytest.raises(module.UnsupportedBasisError, match="Mixed Cartesian"):
        module.BasisSetEngine(data)
