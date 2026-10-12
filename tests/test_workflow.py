from __future__ import annotations

import unittest

from dsipm.workflow import DEFAULT_WORKFLOW_PATH, load_sensor_workflow


class SensorWorkflowTests(unittest.TestCase):
    def test_canonical_workflow_loads(self):
        workflow = load_sensor_workflow()
        self.assertEqual(workflow.schema_version, 1)
        self.assertEqual(workflow.source_path, DEFAULT_WORKFLOW_PATH.resolve())
        self.assertEqual(len(workflow.sha256), 64)

    def test_reference_point_is_in_every_relevant_scan(self):
        workflow = load_sensor_workflow()
        self.assertIn(
            workflow.reference.spad_pitch_um,
            workflow.scans.full_detector_spad_pitch_um,
        )
        self.assertIn(
            workflow.reference.spad_pitch_um,
            workflow.scans.one_fiber_spad_pitch_um,
        )
        self.assertIn(
            workflow.reference.spad_recovery_ns,
            workflow.scans.full_detector_recovery_ns,
        )
        self.assertIn(
            workflow.factorized_proxy.fill_factor,
            workflow.scans.fill_factor,
        )


if __name__ == "__main__":
    unittest.main()
