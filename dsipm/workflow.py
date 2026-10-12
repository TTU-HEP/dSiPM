"""Canonical sensor assumptions shared by every CPAD/JINST plot workflow."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any


DEFAULT_WORKFLOW_PATH = Path(__file__).resolve().parents[1] / "sensor_workflow.json"


@dataclass(frozen=True)
class ReferenceScenario:
    name: str
    detection_model: str
    effective_pde: float
    spad_pitch_um: float
    spad_recovery_ns: float
    target_readout_channel_side_um: float
    timestamp_buffer_depth: int
    channel_resolving_time_ps: float
    count_bits: int
    seed: int
    note: str

    def validate(self) -> None:
        if self.detection_model != "constant_effective_pde":
            raise ValueError(f"Unsupported reference detection model {self.detection_model!r}")
        if not 0.0 <= self.effective_pde <= 1.0:
            raise ValueError("effective_pde must be between zero and one")
        if self.spad_pitch_um <= 0 or self.spad_recovery_ns < 0:
            raise ValueError("SPAD pitch must be positive and recovery non-negative")
        if self.target_readout_channel_side_um <= 0:
            raise ValueError("target readout-channel side must be positive")
        if self.timestamp_buffer_depth < 0 or self.channel_resolving_time_ps < 0:
            raise ValueError("readout timing resources cannot be negative")
        if self.count_bits < 0:
            raise ValueError("count_bits cannot be negative")


@dataclass(frozen=True)
class FactorizedDetectionProxy:
    name: str
    active_area_shape: str
    fill_factor: float
    junction_peak_probability: float
    spectral_model: str
    note: str

    def validate(self) -> None:
        if self.active_area_shape not in {"square", "circle"}:
            raise ValueError(f"Unsupported active-area shape {self.active_area_shape!r}")
        if not 0.0 <= self.fill_factor <= 1.0:
            raise ValueError("fill_factor must be between zero and one")
        if not 0.0 <= self.junction_peak_probability <= 1.0:
            raise ValueError("junction_peak_probability must be between zero and one")
        if self.spectral_model not in {"flat", "legacy_shape"}:
            raise ValueError(f"Unsupported spectral model {self.spectral_model!r}")


@dataclass(frozen=True)
class ScanGrids:
    full_detector_spad_pitch_um: tuple[float, ...]
    full_detector_recovery_ns: tuple[float, ...]
    one_fiber_spad_pitch_um: tuple[float, ...]
    fill_factor: tuple[float, ...]
    spads_per_readout_channel_side: tuple[int, ...]

    @classmethod
    def from_dict(cls, values: dict[str, Any]) -> "ScanGrids":
        return cls(
            full_detector_spad_pitch_um=tuple(values["full_detector_spad_pitch_um"]),
            full_detector_recovery_ns=tuple(values["full_detector_recovery_ns"]),
            one_fiber_spad_pitch_um=tuple(values["one_fiber_spad_pitch_um"]),
            fill_factor=tuple(values["fill_factor"]),
            spads_per_readout_channel_side=tuple(
                values["spads_per_readout_channel_side"]
            ),
        )

    def validate(self) -> None:
        for name in (
            "full_detector_spad_pitch_um",
            "one_fiber_spad_pitch_um",
        ):
            values = getattr(self, name)
            if not values or any(value <= 0 for value in values):
                raise ValueError(f"{name} must contain positive values")
        if not self.full_detector_recovery_ns or any(
            value < 0 for value in self.full_detector_recovery_ns
        ):
            raise ValueError("full_detector_recovery_ns must be non-negative")
        if not self.fill_factor or any(
            not 0.0 <= value <= 1.0 for value in self.fill_factor
        ):
            raise ValueError("fill_factor scan must stay between zero and one")
        if not self.spads_per_readout_channel_side or any(
            value < 1 for value in self.spads_per_readout_channel_side
        ):
            raise ValueError("SPAD grouping factors must be positive integers")


@dataclass(frozen=True)
class SensorWorkflow:
    schema_version: int
    reference: ReferenceScenario
    factorized_proxy: FactorizedDetectionProxy
    scans: ScanGrids
    source_path: Path
    sha256: str


def load_sensor_workflow(path: str | Path = DEFAULT_WORKFLOW_PATH) -> SensorWorkflow:
    source = Path(path).expanduser().resolve()
    payload = source.read_bytes()
    raw = json.loads(payload)
    if int(raw.get("schema_version", -1)) != 1:
        raise ValueError("Unsupported sensor-workflow schema version")
    reference = ReferenceScenario(**raw["reference_scenario"])
    factorized = FactorizedDetectionProxy(**raw["factorized_detection_proxy"])
    scans = ScanGrids.from_dict(raw["scan_grids"])
    reference.validate()
    factorized.validate()
    scans.validate()
    return SensorWorkflow(
        schema_version=1,
        reference=reference,
        factorized_proxy=factorized,
        scans=scans,
        source_path=source,
        sha256=hashlib.sha256(payload).hexdigest(),
    )
