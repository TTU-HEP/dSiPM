from __future__ import annotations

import unittest

import numpy as np

from dsipm.full_detector import FullDetectorPhotonEvent
from dsipm.reduction import HistogramAccumulator, reduce_event, scenario_column


class ReductionTests(unittest.TestCase):
    def _event(self) -> FullDetectorPhotonEvent:
        return FullDetectorPhotonEvent(
            event_index=17,
            beam_energy_gev=100.0,
            c_signal=12.0,
            s_signal=18.0,
            ecalo_truth_gev=99.0,
            eleak_truth_gev=1.0,
            raw_photon_count=4,
            local_x_mm=np.asarray([0.01, 0.011, 0.50, -0.40]),
            local_y_mm=np.asarray([0.01, 0.011, 0.50, -0.40]),
            sensor_x_mm=np.asarray([0.01, 0.011, 4.50, -4.40]),
            sensor_y_mm=np.asarray([0.01, 0.011, 4.50, -4.40]),
            rod_ix=np.asarray([0, 0, 1, -1]),
            rod_iy=np.asarray([0, 0, 1, -1]),
            sensor_id=np.asarray([10, 10, 11, 12]),
            time_ns=np.asarray([10.0, 10.1, 12.0, 20.0]),
            wavelength_nm=np.asarray([420.0, 420.0, 430.0, 440.0]),
            production_z_mm=np.asarray([100.0, 200.0, 300.0, 400.0]),
        )

    def test_event_row_closes_reference_scan(self) -> None:
        event = self._event()
        row, detected, accepted, lost = reduce_event(
            event,
            sample="eplus_100gev",
            source_root="input.root",
            source_entry=0,
            pde=1.0,
            pitches_um=[20.0],
            recovery_times_ns=[4.3],
            reference_pitch_um=20.0,
            reference_recovery_ns=4.3,
            seed=4307,
        )
        self.assertEqual(row["event"], 17)
        self.assertEqual(row["detected_reference"], 4)
        self.assertEqual(row["accepted_reference"], int(np.count_nonzero(accepted)))
        self.assertEqual(
            row[scenario_column(20.0, 4.3)],
            row["accepted_reference"],
        )
        self.assertEqual(
            int(np.count_nonzero(detected)),
            int(np.count_nonzero(accepted)) + int(np.count_nonzero(lost)),
        )

    def test_histogram_payload_is_additive_and_closes(self) -> None:
        event = self._event()
        _, detected, accepted, lost = reduce_event(
            event,
            sample="eplus_100gev",
            source_root="input.root",
            source_entry=0,
            pde=1.0,
            pitches_um=[20.0],
            recovery_times_ns=[4.3],
            reference_pitch_um=20.0,
            reference_recovery_ns=4.3,
            seed=4307,
        )
        histograms = HistogramAccumulator()
        histograms.add(event, detected, accepted, lost)
        closure = histograms.closure()
        self.assertEqual(closure["incident_photons_in_time_histogram"], 4)
        self.assertEqual(closure["detected_photons_in_time_histogram"], 4)
        self.assertEqual(
            closure["detected_photons_in_time_histogram"],
            closure["accepted_photons_in_time_histogram"]
            + closure["lost_recovery_photons_in_time_histogram"],
        )
        payload = histograms.payload()
        self.assertEqual(payload["local_xy_incident"].shape, (120, 120))
        self.assertEqual(payload["global_xy_incident"].shape, (200, 200))


if __name__ == "__main__":
    unittest.main()
