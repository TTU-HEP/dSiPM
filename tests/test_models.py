from __future__ import annotations
import unittest
import numpy as np
from dsipm.sensors import (AnalogSiPMScenario, FullDetectorPhotonEvent, FullDetectorScenario,
    digitize_analog_sipm_event, digitize_full_detector_event, global_cell_ids,
    histogram_observables, partial_recharge_weights, scope_kernel, scope_waveform,
    ChannelReadoutConfig, channel_readout)
from dsipm.geometry import (ORIGINAL_DREAM_FIBER_DIAMETER_MM, ORIGINAL_DREAM_TOWER_RADIUS_MM,
    ORIGINAL_DREAM_TOWER_SPACING_MM, hgdream_cherenkov_channel_ids, hgdream_cherenkov_channel_layout,
    original_dream_pmt_centers, original_dream_pmt_channel_ids)

def make_event(*, sensor_id, x_mm, time_ns, rod_ix=None, rod_iy=None) -> FullDetectorPhotonEvent:
    x = np.asarray(x_mm, dtype=np.float64)
    time = np.asarray(time_ns, dtype=np.float64)
    sensor = np.asarray(sensor_id, dtype=np.int64)
    zeros = np.zeros(time.size, dtype=np.float64)
    rod_x = (
        np.zeros(time.size, dtype=np.int64)
        if rod_ix is None
        else np.asarray(rod_ix, dtype=np.int64)
    )
    rod_y = (
        np.zeros(time.size, dtype=np.int64)
        if rod_iy is None
        else np.asarray(rod_iy, dtype=np.int64)
    )
    return FullDetectorPhotonEvent(
        event_index=4,
        beam_energy_gev=100.0,
        c_signal=10.0,
        s_signal=1.0,
        ecalo_truth_gev=99.0,
        eleak_truth_gev=1.0,
        raw_photon_count=time.size,
        local_x_mm=x,
        local_y_mm=zeros,
        sensor_x_mm=zeros,
        sensor_y_mm=zeros,
        rod_ix=rod_x,
        rod_iy=rod_y,
        sensor_id=sensor,
        time_ns=time,
        wavelength_nm=np.full(time.size, 450.0),
        production_z_mm=zeros,
    )


class SensorComparisonTests(unittest.TestCase):
    def test_hgdream_cherenkov_layout_matches_installed_channel_counts(self):
        layout = hgdream_cherenkov_channel_layout()
        self.assertEqual(layout.n_channels, 448)
        self.assertEqual(int(np.count_nonzero(layout.is_core)), 64)
        self.assertEqual(int(np.count_nonzero(~layout.is_core)), 384)
        np.testing.assert_array_equal(layout.rods_per_channel[layout.is_core], 3)
        np.testing.assert_array_equal(layout.rods_per_channel[~layout.is_core], 12)
        np.testing.assert_array_equal(layout.sensor_width_mm[layout.is_core], 3.0)
        np.testing.assert_array_equal(layout.sensor_width_mm[~layout.is_core], 6.0)

    def test_readout_projections_assign_central_photons(self):
        event = make_event(sensor_id=[10, 11], x_mm=[0.0, 0.0], time_ns=[1.0, 1.1])
        layout = hgdream_cherenkov_channel_layout()
        hgdream_ids = hgdream_cherenkov_channel_ids(event, layout)
        self.assertTrue(np.all(hgdream_ids >= 0))
        centers = original_dream_pmt_centers()
        self.assertEqual(centers.shape, (19, 2))
        center_id = int(np.argmin(np.sum(centers**2, axis=1)))
        np.testing.assert_array_equal(
            original_dream_pmt_channel_ids(event, centers),
            center_id,
        )
        nearest_distance = np.min(
            np.linalg.norm(
                centers[:, np.newaxis, :] - centers[np.newaxis, :, :],
                axis=2,
            )
            + np.eye(19) * 1.0e9,
        )
        self.assertAlmostEqual(nearest_distance, ORIGINAL_DREAM_TOWER_SPACING_MM)
        self.assertAlmostEqual(2.0 * ORIGINAL_DREAM_TOWER_RADIUS_MM, 80.0)
        self.assertAlmostEqual(ORIGINAL_DREAM_FIBER_DIAMETER_MM, 0.8)

    def test_original_dream_projection_rejects_photons_outside_detector(self):
        event = make_event(
            sensor_id=[10, 11],
            x_mm=[0.0, 0.0],
            time_ns=[1.0, 1.1],
            rod_ix=[0, 60],
            rod_iy=[0, 0],
        )
        channel_ids = original_dream_pmt_channel_ids(event)
        self.assertGreaterEqual(int(channel_ids[0]), 0)
        self.assertEqual(int(channel_ids[1]), -1)

    def test_same_local_spad_in_different_fibers_has_different_global_id(self):
        event = make_event(sensor_id=[10, 11], x_mm=[0.001, 0.001], time_ns=[1.0, 1.1])
        ids = global_cell_ids(event, 20.0)
        self.assertNotEqual(int(ids[0]), int(ids[1]))

    def test_dead_time_is_independent_between_physical_fibers(self):
        event = make_event(sensor_id=[10, 11], x_mm=[0.001, 0.001], time_ns=[1.0, 1.1])
        scenario = FullDetectorScenario(pde=1.0, pitch_um=20.0, recovery_ns=40.0)
        result = digitize_full_detector_event(event, scenario, seed=7)
        self.assertEqual(result.n_detected, 2)
        self.assertEqual(result.n_accepted, 2)

    def test_analog_sipm_partial_recharge_has_expected_charge(self):
        tau = 82.0
        event = make_event(
            sensor_id=[10, 10],
            x_mm=[0.001, 0.001],
            time_ns=[0.0, tau * np.log(2.0)],
        )
        scenario = AnalogSiPMScenario(
            name="test",
            pde=1.0,
            microcell_size_um=35.0,
            n_microcells=4774,
            recovery_tau_ns=tau,
        )
        result = digitize_analog_sipm_event(event, scenario, seed=7)
        self.assertEqual(result.n_detected, 2)
        np.testing.assert_allclose(result.charge_per_detected_photon, [1.0, 0.5])
        self.assertAlmostEqual(result.charge_pe, 1.5)

    def test_analog_sipm_recharge_is_independent_between_fibers(self):
        event = make_event(sensor_id=[10, 11], x_mm=[0.001, 0.001], time_ns=[1.0, 1.1])
        scenario = AnalogSiPMScenario(
            name="test",
            pde=1.0,
            microcell_size_um=35.0,
            n_microcells=4774,
            recovery_tau_ns=82.0,
        )
        result = digitize_analog_sipm_event(event, scenario, seed=7)
        self.assertEqual(result.n_detected, 2)
        self.assertAlmostEqual(result.charge_pe, 2.0)

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


class SharedReadoutTests(unittest.TestCase):
    def test_channel_buffer_is_independent_for_each_readout_channel(self):
        result = channel_readout(
            np.asarray([0, 1, 2, 3]),
            np.asarray([0, 0, 0, 0]),
            np.asarray([1.0, 2.0, 3.0, 4.0]),
            n_spads_axis=4,
            config=ChannelReadoutConfig(
                spads_per_channel_side=2,
                timestamp_buffer_depth=1,
                channel_resolving_time_ns=0.0,
                count_bits=8,
                timestamp_bits=12,
                address_bits=4,
                event_overhead_bits=16,
            ),
        )
        self.assertEqual(result.active_channels, 2)
        self.assertEqual(result.counted_hits, 4)
        self.assertEqual(result.timestamp_records, 2)
        self.assertEqual(result.overflowed_channels, 2)
        self.assertEqual(result.payload_bits, 16 + 2 * (4 + 8) + 2 * 12)

    def test_channel_resolving_time_and_count_saturation_are_separate(self):
        result = channel_readout(
            np.asarray([0, 1, 2, 3]),
            np.asarray([0, 0, 0, 0]),
            np.asarray([1.00, 1.02, 1.08, 1.20]),
            n_spads_axis=4,
            config=ChannelReadoutConfig(
                spads_per_channel_side=4,
                timestamp_buffer_depth=8,
                channel_resolving_time_ns=0.05,
                count_bits=2,
            ),
        )
        self.assertEqual(result.accepted_hits, 4)
        self.assertEqual(result.counted_hits, 3)
        self.assertEqual(result.count_saturated_channels, 1)
        self.assertEqual(result.timestamp_candidates, 3)
        self.assertEqual(result.resolving_lost_hits, 1)
        self.assertEqual(result.timestamp_records, 3)

if __name__ == '__main__':
    unittest.main()
