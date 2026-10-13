from __future__ import annotations

import unittest

import numpy as np

from dsipm.pulse_shape import (
    histogram_observables,
    partial_recharge_weights,
    scope_kernel,
    scope_waveform,
)


class PulseShapeTests(unittest.TestCase):
    def test_partial_recharge_tracks_each_cell_independently(self):
        tau = 10.0
        cells = np.asarray([2, 1, 2, 1])
        times = np.asarray([0.0, 0.0, tau * np.log(2.0), 100.0])
        weights = partial_recharge_weights(cells, times, tau)
        np.testing.assert_allclose(weights, [1.0, 1.0, 0.5, 1.0 - np.exp(-10.0)])

    def test_scope_kernel_has_unit_area(self):
        kernel = scope_kernel(0.01)
        self.assertAlmostEqual(float(np.sum(kernel)), 1.0)
        self.assertTrue(np.all(kernel >= 0.0))

    def test_scope_waveform_preserves_integral_away_from_right_edge(self):
        impulses = np.zeros(2000)
        impulses[100] = 2.0
        impulses[300] = 3.0
        waveform = scope_waveform(impulses, bin_width_ns=0.01)
        self.assertAlmostEqual(float(np.sum(waveform)), 5.0, places=10)

    def test_histogram_observables(self):
        hist = np.asarray([1.0, 2.0, 1.0, 0.0])
        edges = np.arange(5, dtype=np.float64)
        values = histogram_observables(hist, edges, prompt_end_ns=2.0)
        self.assertEqual(values["integral"], 4.0)
        self.assertAlmostEqual(values["prompt_fraction"], 0.75)
        self.assertAlmostEqual(values["mean_time_ns"], 1.5)
        self.assertAlmostEqual(values["central80_width_ns"], 2.0)


if __name__ == "__main__":
    unittest.main()
