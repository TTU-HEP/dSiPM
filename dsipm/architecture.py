"""Channel-level readout model for FLASH architecture trade studies.

SPAD recovery and readout-channel capacity are distinct effects.  A photon can
produce an accepted avalanche while its timestamp is subsequently merged or
dropped by electronics shared across several SPADs.  This module keeps the
count and timestamp paths separate so a large timing channel cannot appear
lossless merely because its ideal counter still integrates every avalanche.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .digitizer import _apply_nonparalyzable_dead_time


@dataclass(frozen=True)
class ChannelReadoutConfig:
    """Resources shared by a square group of SPADs."""

    spads_per_channel_side: int
    timestamp_buffer_depth: int = 8
    channel_resolving_time_ns: float = 0.05
    count_bits: int = 8
    timestamp_bits: int = 12
    address_bits: int = 12
    event_overhead_bits: int = 64

    def validate(self) -> None:
        if self.spads_per_channel_side < 1:
            raise ValueError("spads_per_channel_side must be positive")
        if self.timestamp_buffer_depth < 0:
            raise ValueError("timestamp_buffer_depth cannot be negative")
        if self.channel_resolving_time_ns < 0:
            raise ValueError("channel_resolving_time_ns cannot be negative")
        for field_name in (
            "count_bits",
            "timestamp_bits",
            "address_bits",
            "event_overhead_bits",
        ):
            if getattr(self, field_name) < 0:
                raise ValueError(f"{field_name} cannot be negative")


@dataclass(frozen=True)
class ChannelReadoutResult:
    """Count and timing output from one event after SPAD recovery."""

    accepted_hits: int
    active_channels: int
    counted_hits: int
    count_saturated_channels: int
    timestamp_candidates: int
    timestamp_records: int
    resolving_lost_hits: int
    buffer_lost_records: int
    overflowed_channels: int
    max_hits_in_channel: int
    payload_bits: int

    @property
    def count_survival(self) -> float:
        return self.counted_hits / self.accepted_hits if self.accepted_hits else 1.0

    @property
    def timestamp_survival(self) -> float:
        return (
            self.timestamp_records / self.accepted_hits
            if self.accepted_hits
            else 1.0
        )

    @property
    def any_buffer_overflow(self) -> bool:
        return self.overflowed_channels > 0


def _keep_first_per_group(
    group_ids: np.ndarray,
    maximum_per_group: int,
) -> np.ndarray:
    """Return a mask retaining the first ``maximum_per_group`` group entries."""

    if group_ids.size == 0 or maximum_per_group <= 0:
        return np.zeros(group_ids.size, dtype=bool)
    order = np.argsort(group_ids, kind="stable")
    sorted_groups = group_ids[order]
    starts = np.r_[0, np.flatnonzero(sorted_groups[1:] != sorted_groups[:-1]) + 1]
    ends = np.r_[starts[1:], sorted_groups.size]
    keep = np.zeros(group_ids.size, dtype=bool)
    for start, end in zip(starts, ends, strict=True):
        keep[order[start : min(end, start + maximum_per_group)]] = True
    return keep


def channel_readout(
    spad_ix: np.ndarray,
    spad_iy: np.ndarray,
    time_ns: np.ndarray,
    *,
    n_spads_axis: int,
    config: ChannelReadoutConfig,
) -> ChannelReadoutResult:
    """Apply shared-channel resolving time, counters, buffers, and payload.

    Inputs are avalanches already accepted by the individual-SPAD recovery
    model.  Each active channel exports one address and count word plus up to
    ``timestamp_buffer_depth`` timestamp words.  The counter integrates every
    accepted avalanche until its finite bit depth saturates.
    """

    config.validate()
    spad_ix = np.asarray(spad_ix, dtype=np.int64)
    spad_iy = np.asarray(spad_iy, dtype=np.int64)
    time_ns = np.asarray(time_ns, dtype=np.float64)
    if spad_ix.shape != spad_iy.shape or spad_ix.shape != time_ns.shape:
        raise ValueError("SPAD indices and time_ns must have the same shape")
    if n_spads_axis < 1:
        raise ValueError("n_spads_axis must be positive")
    if spad_ix.size and (
        np.any(spad_ix < 0)
        or np.any(spad_iy < 0)
        or np.any(spad_ix >= n_spads_axis)
        or np.any(spad_iy >= n_spads_axis)
    ):
        raise ValueError("SPAD indices fall outside the sensor grid")

    accepted_hits = int(time_ns.size)
    if not accepted_hits:
        return ChannelReadoutResult(
            accepted_hits=0,
            active_channels=0,
            counted_hits=0,
            count_saturated_channels=0,
            timestamp_candidates=0,
            timestamp_records=0,
            resolving_lost_hits=0,
            buffer_lost_records=0,
            overflowed_channels=0,
            max_hits_in_channel=0,
            payload_bits=config.event_overhead_bits,
        )

    order = np.argsort(time_ns, kind="stable")
    time_ns = time_ns[order]
    spad_ix = spad_ix[order]
    spad_iy = spad_iy[order]

    group_side = config.spads_per_channel_side
    n_channels_axis = int(np.ceil(n_spads_axis / group_side))
    channel_ids = (
        (spad_iy // group_side) * n_channels_axis + (spad_ix // group_side)
    )
    _, channel_hit_counts = np.unique(channel_ids, return_counts=True)
    active_channels = int(channel_hit_counts.size)
    max_hits = int(np.max(channel_hit_counts))

    maximum_count = (1 << config.count_bits) - 1 if config.count_bits else 0
    counted_hits = int(np.minimum(channel_hit_counts, maximum_count).sum())
    count_saturated = int(np.count_nonzero(channel_hit_counts > maximum_count))

    resolving_mask = _apply_nonparalyzable_dead_time(
        channel_ids,
        time_ns,
        config.channel_resolving_time_ns,
    )
    candidate_channels = channel_ids[resolving_mask]
    timestamp_candidates = int(candidate_channels.size)
    _, candidate_counts = np.unique(candidate_channels, return_counts=True)
    overflowed_channels = int(
        np.count_nonzero(candidate_counts > config.timestamp_buffer_depth)
    )
    buffer_mask = _keep_first_per_group(
        candidate_channels,
        config.timestamp_buffer_depth,
    )
    timestamp_records = int(np.count_nonzero(buffer_mask))

    payload_bits = (
        config.event_overhead_bits
        + active_channels * (config.address_bits + config.count_bits)
        + timestamp_records * config.timestamp_bits
    )
    return ChannelReadoutResult(
        accepted_hits=accepted_hits,
        active_channels=active_channels,
        counted_hits=counted_hits,
        count_saturated_channels=count_saturated,
        timestamp_candidates=timestamp_candidates,
        timestamp_records=timestamp_records,
        resolving_lost_hits=accepted_hits - timestamp_candidates,
        buffer_lost_records=timestamp_candidates - timestamp_records,
        overflowed_channels=overflowed_channels,
        max_hits_in_channel=max_hits,
        payload_bits=int(payload_bits),
    )
