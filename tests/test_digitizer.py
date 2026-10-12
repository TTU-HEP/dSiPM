from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import unittest

import numpy as np

from dsipm.architecture import ChannelReadoutConfig, channel_readout
from dsipm.active_area import (
    active_area_mask,
    active_width_um,
    junction_pde_probability,
)
from dsipm.characterization import (
    expected_unique_cells,
    scan_pde_pitch_recovery,
)
from dsipm.cpad_plots import pitch_only_projection
from dsipm.config import load_device_catalog
from dsipm.digitizer import (
    PhotonEvent,
    _cell_indices,
    digitize_event,
    legacy_pde,
)
from dsipm.requirements import (
    accepted_counts_by_recovery,
    interpolated_pitch_limit,
    scan_pitch_and_recovery,
)
from dsipm.root_io import transform_coordinates


REPO_DIR = Path(__file__).resolve().parents[1]


def make_event(x_mm, time_ns, wavelength_nm=None) -> PhotonEvent:
    x = np.asarray(x_mm, dtype=np.float64)
    time = np.asarray(time_ns, dtype=np.float64)
    if wavelength_nm is None:
        wavelength_nm = np.full(time.size, 450.0)
    return PhotonEvent(
        event_index=7,
        raw_photon_count=int(time.size),
        x_mm=x,
        y_mm=np.zeros(time.size, dtype=np.float64),
        time_ns=time,
        wavelength_nm=np.asarray(wavelength_nm, dtype=np.float64),
        production_z_mm=np.zeros(time.size, dtype=np.float64),
    )


class DigitizerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = load_device_catalog(REPO_DIR / "device_configs.json")

    def test_catalog_profiles_resolve(self):
        for profile in self.catalog.profiles:
            selected = self.catalog.select(profile=profile)
            self.assertTrue(selected, profile)

    def test_cpad_pitch_projection_interpolates_historical_points(self):
        pitch = np.asarray([10.0, 25.0, 50.0])
        efficiency = np.asarray([0.99, 0.90, 0.70])
        self.assertAlmostEqual(
            pitch_only_projection(37.5, pitch, efficiency),
            0.80,
        )
        with self.assertRaises(ValueError):
            pitch_only_projection(5.0, pitch, efficiency)

    def test_cpad_benchmark_pitches_match_device_catalog(self):
        with (REPO_DIR / "cpad_benchmarks.json").open() as handle:
            benchmarks = json.load(handle)
        for benchmark in benchmarks["competitors"]:
            configured = self.catalog.devices[benchmark["name"]]
            self.assertAlmostEqual(
                float(benchmark["pitch_um"]),
                configured.spad_pitch_um,
                msg=benchmark["name"],
            )

    def test_legacy_pde_is_bounded_and_zero_outside_fit_range(self):
        values = legacy_pde(np.asarray([250.0, 400.0, 450.0, 620.0, 700.0]))
        self.assertEqual(values[0], 0.0)
        self.assertEqual(values[-1], 0.0)
        self.assertTrue(np.all((values >= 0.0) & (values <= 1.0)))
        self.assertAlmostEqual(values[2], 0.2146637262, places=8)

    def test_square_active_area_separates_fill_factor_from_pitch(self):
        mask = active_area_mask(
            np.asarray([0.005, 0.001, 0.009, 0.015]),
            np.asarray([0.005, 0.001, 0.009, 0.015]),
            pitch_um=10.0,
            fill_factor=0.25,
            detector_min_mm=0.0,
            detector_max_mm=0.02,
        )
        np.testing.assert_array_equal(mask, [True, False, False, True])
        self.assertAlmostEqual(active_width_um(20.0, 0.25), 10.0)

    def test_circular_active_area_rejects_impossible_fill_factor(self):
        with self.assertRaises(ValueError):
            active_width_um(20.0, 0.90, shape="circle")

    def test_spectral_proxy_has_requested_peak_and_documented_support(self):
        wavelength = np.linspace(250.0, 700.0, 4501)
        probability = junction_pde_probability(
            wavelength,
            peak_probability=0.30,
        )
        self.assertAlmostEqual(float(np.max(probability)), 0.30, places=8)
        self.assertEqual(float(probability[0]), 0.0)
        self.assertEqual(float(probability[-1]), 0.0)

    def test_event_long_deadtime_keeps_first_photon_per_spad(self):
        config = self.catalog.devices["legacy_20um"]
        event = make_event([0.001, 0.002, 0.021], [1.0, 2.0, 3.0])
        result = digitize_event(
            event,
            config,
            seed=1,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
        )
        self.assertEqual(result.n_detected, 3)
        self.assertEqual(result.n_accepted, 2)
        self.assertEqual(result.n_lost_deadtime, 1)
        np.testing.assert_allclose(result.accepted_time_ns, [1.0, 3.0])

    def test_finite_recovery_allows_a_later_hit(self):
        base = self.catalog.devices["legacy_10um"]
        config = replace(base, dead_time_ns=4.3)
        event = make_event([0.001, 0.001, 0.001], [1.0, 3.0, 6.0])
        result = digitize_event(
            event,
            config,
            seed=1,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
        )
        self.assertEqual(result.n_accepted, 2)
        np.testing.assert_allclose(result.accepted_time_ns, [1.0, 6.0])

    def test_historical_three_mm_sensor_scale_matches_slide_segmentation(self):
        historical = self.catalog.devices["legacy_10um"]
        unscaled = replace(historical, coordinate_scale=1.0)
        event = make_event([0.001, 0.005], [1.0, 2.0])
        historical_result = digitize_event(
            event,
            historical,
            seed=1,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
        )
        unscaled_result = digitize_event(
            event,
            unscaled,
            seed=1,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
        )
        self.assertEqual(historical_result.n_accepted, 2)
        self.assertEqual(unscaled_result.n_accepted, 1)

    def test_common_random_numbers_make_pde_scan_monotonic(self):
        base = self.catalog.devices["legacy_10um"]
        low = replace(base, pde_model="constant", pde_value=0.2, dead_time_ns=None)
        high = replace(base, pde_model="constant", pde_value=0.5, dead_time_ns=None)
        event = make_event(
            np.linspace(0.001, 0.999, 1000),
            np.linspace(1.0, 30.0, 1000),
        )
        low_result = digitize_event(
            event,
            low,
            seed=123,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
        )
        high_result = digitize_event(
            event,
            high,
            seed=123,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
        )
        self.assertLessEqual(low_result.n_detected, high_result.n_detected)

    def test_first_last_policy_exports_no_more_than_two_records(self):
        base = self.catalog.devices["aspides_2025"]
        config = replace(
            base,
            pde_model="unity",
            pde_value=None,
            timing_fwhm_ps=None,
            tdc_bin_ps=None,
        )
        event = make_event([0.01, 0.04, 0.08], [1.0, 2.0, 3.0])
        result = digitize_event(
            event,
            config,
            seed=1,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
        )
        self.assertEqual(result.n_timing_records, 2)
        np.testing.assert_allclose(result.timestamp_ns, [1.0, 3.0])

    def test_timestamp_policy_uses_observed_time_after_jitter(self):
        base = self.catalog.devices["aspides_2025"]
        config = replace(
            base,
            pde_model="unity",
            pde_value=None,
            dead_time_ns=None,
            timing_fwhm_ps=5000.0,
            tdc_bin_ps=None,
            timestamp_policy="first_global",
        )
        event = make_event([0.01, 0.04, 0.08], [1.0, 1.1, 1.2])
        result = digitize_event(
            event,
            config,
            seed=3,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
        )
        self.assertEqual(result.n_timing_records, 1)
        self.assertAlmostEqual(result.timestamp_ns[0], np.min(result.accepted_time_ns))

    def test_active_timing_groups_is_not_timestamp_record_count(self):
        base = self.catalog.devices["aspides_2025"]
        config = replace(
            base,
            pde_model="unity",
            pde_value=None,
            dead_time_ns=None,
            timing_fwhm_ps=None,
            tdc_bin_ps=None,
            timestamp_policy="first_global",
            timing_group_spads_x=1,
            timing_group_spads_y=1,
        )
        event = make_event([0.01, 0.04, 0.08], [1.0, 2.0, 3.0])
        result = digitize_event(
            event,
            config,
            seed=1,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
        )
        self.assertEqual(result.n_timing_records, 1)
        self.assertEqual(result.n_active_timing_groups, 3)

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

    def test_cell_assignment_rejects_out_of_sensor_coordinates(self):
        with self.assertRaises(ValueError):
            _cell_indices(
                np.asarray([-0.01, 0.5]),
                np.asarray([0.5, 0.5]),
                0.01,
                0.0,
                1.0,
            )

    def test_old_geometry_is_only_shift_and_unit_conversion(self):
        x, y = transform_coordinates(
            np.asarray([-3.2908]),
            np.asarray([3.6878]),
            np.asarray([0]),
            "old_2025",
        )
        np.testing.assert_allclose(x, [0.0], atol=1e-12)
        np.testing.assert_allclose(y, [0.0], atol=1e-12)

    def test_recovery_scan_counts_repeated_spad_hits(self):
        cell_ids = np.asarray([4, 4, 4, 9], dtype=np.int64)
        time_ns = np.asarray([1.0, 3.0, 6.0, 2.0])
        counts = accepted_counts_by_recovery(cell_ids, time_ns, [0.0, 2.0, 4.3, 40.0])
        np.testing.assert_array_equal(counts, [4, 4, 3, 2])

    def test_requirement_scan_and_interpolated_limit(self):
        event = make_event([0.001, 0.002, 0.021], [1.0, 2.0, 3.0])
        points = scan_pitch_and_recovery(
            [event],
            pitches_um=[10.0, 30.0],
            recovery_times_ns=[40.0],
            coordinate_scale=1.0,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
        )
        self.assertEqual(points[0].accepted_photons, 2)
        self.assertEqual(points[1].accepted_photons, 1)
        self.assertAlmostEqual(
            interpolated_pitch_limit(
                points,
                recovery_ns=40.0,
                target_retention=0.5,
            ),
            20.0,
        )

    def test_uniform_flash_saturation_matches_exact_small_cases(self):
        expectation = expected_unique_cells([0, 1, 2], n_cells=10)
        np.testing.assert_allclose(expectation, [0.0, 1.0, 1.9])
        np.testing.assert_allclose(
            expected_unique_cells([0, 1, 100], n_cells=1),
            [0.0, 1.0, 1.0],
        )
        with self.assertRaises(ValueError):
            expected_unique_cells([1], n_cells=0)

    def test_pde_pitch_recovery_scan_separates_yield_from_survival(self):
        event = make_event(
            np.linspace(0.001, 0.999, 1000),
            np.linspace(1.0, 30.0, 1000),
        )
        points = scan_pde_pitch_recovery(
            [event],
            pdes=[0.2, 0.5, 1.0],
            pitches_um=[1.0],
            recovery_times_ns=[40.0],
            seed=123,
            coordinate_scale=1.0,
            detector_min_mm=0.0,
            detector_max_mm=1.0,
            bootstrap_samples=0,
        )
        self.assertEqual(len(points), 3)
        self.assertTrue(
            all(
                point.accepted_photons
                <= point.detected_photons
                <= point.incident_photons
                for point in points
            )
        )
        self.assertLessEqual(points[0].detected_photons, points[1].detected_photons)
        self.assertLessEqual(points[1].detected_photons, points[2].detected_photons)
        self.assertLessEqual(points[0].usable_fraction, points[1].usable_fraction)
        self.assertLessEqual(points[1].usable_fraction, points[2].usable_fraction)
        self.assertGreater(points[0].recovery_retention, 0.95)


if __name__ == "__main__":
    unittest.main()
