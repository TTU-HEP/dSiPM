"""Device configuration loading and validation.

The catalog deliberately records whether a number is measured, designed,
targeted, or merely a simulation proxy.  Device comparisons should never turn
an unknown public parameter into an apparently measured value.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any


PDE_MODELS = {"unity", "constant", "legacy_polynomial"}
TIMESTAMP_POLICIES = {
    "all",
    "first_global",
    "first_last_global",
    "first_per_group",
    "histogram_per_group",
}


@dataclass(frozen=True)
class DeviceConfig:
    """Parameters needed by the first device-aware response model."""

    name: str
    label: str
    category: str
    maturity: str
    source_url: str
    spad_pitch_um: float
    timing_group_spads_x: int
    timing_group_spads_y: int
    dead_time_ns: float | None
    pde_model: str
    pde_value: float | None
    pde_scale: float
    pde_basis: str
    timing_fwhm_ps: float | None
    tdc_bin_ps: float | None
    timestamp_policy: str
    histogram_bins: int
    count_bits: int
    timestamp_bits: int
    address_bits: int
    event_overhead_bits: int
    notes: str
    coordinate_scale: float = 1.0

    @classmethod
    def from_dict(cls, name: str, values: dict[str, Any]) -> "DeviceConfig":
        config = cls(name=name, **values)
        config.validate()
        return config

    def validate(self) -> None:
        if self.spad_pitch_um <= 0:
            raise ValueError(f"{self.name}: spad_pitch_um must be positive")
        if self.coordinate_scale <= 0:
            raise ValueError(f"{self.name}: coordinate_scale must be positive")
        if self.timing_group_spads_x < 1 or self.timing_group_spads_y < 1:
            raise ValueError(f"{self.name}: timing group dimensions must be positive")
        if self.dead_time_ns is not None and self.dead_time_ns < 0:
            raise ValueError(f"{self.name}: dead_time_ns cannot be negative")
        if self.pde_model not in PDE_MODELS:
            raise ValueError(f"{self.name}: unsupported pde_model {self.pde_model!r}")
        if self.pde_model == "constant" and self.pde_value is None:
            raise ValueError(f"{self.name}: constant PDE requires pde_value")
        if self.pde_value is not None and not 0 <= self.pde_value <= 1:
            raise ValueError(f"{self.name}: pde_value must be between zero and one")
        if self.pde_scale < 0:
            raise ValueError(f"{self.name}: pde_scale cannot be negative")
        if self.timestamp_policy not in TIMESTAMP_POLICIES:
            raise ValueError(
                f"{self.name}: unsupported timestamp_policy {self.timestamp_policy!r}"
            )
        if self.timestamp_policy == "histogram_per_group" and self.histogram_bins < 1:
            raise ValueError(f"{self.name}: histogram mode needs at least one bin")
        for field_name in (
            "histogram_bins",
            "count_bits",
            "timestamp_bits",
            "address_bits",
            "event_overhead_bits",
        ):
            if getattr(self, field_name) < 0:
                raise ValueError(f"{self.name}: {field_name} cannot be negative")

    def to_dict(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.__dict__.items()
        }


@dataclass(frozen=True)
class DeviceCatalog:
    devices: dict[str, DeviceConfig]
    profiles: dict[str, list[str]]
    source_path: Path

    def select(
        self,
        *,
        profile: str | None = None,
        device_names: list[str] | None = None,
    ) -> list[DeviceConfig]:
        if device_names:
            names = device_names
        else:
            selected_profile = profile or "default"
            if selected_profile not in self.profiles:
                available = ", ".join(sorted(self.profiles))
                raise KeyError(
                    f"Unknown profile {selected_profile!r}. Available profiles: {available}"
                )
            names = self.profiles[selected_profile]

        missing = [name for name in names if name not in self.devices]
        if missing:
            raise KeyError(f"Unknown device configuration(s): {', '.join(missing)}")
        return [self.devices[name] for name in names]


def load_device_catalog(path: str | Path) -> DeviceCatalog:
    catalog_path = Path(path).expanduser().resolve()
    with catalog_path.open() as handle:
        raw = json.load(handle)

    devices = {
        name: DeviceConfig.from_dict(name, values)
        for name, values in raw["devices"].items()
    }
    profiles = {
        name: list(device_names)
        for name, device_names in raw["profiles"].items()
    }
    catalog = DeviceCatalog(devices=devices, profiles=profiles, source_path=catalog_path)

    for profile_name, names in profiles.items():
        missing = [name for name in names if name not in devices]
        if missing:
            raise ValueError(
                f"Profile {profile_name!r} refers to missing devices: {', '.join(missing)}"
            )
    return catalog
