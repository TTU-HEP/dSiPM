"""Corrected detector coordinates, readout projections and ideal-fiber reference.

HG-DREAM/DREAM layouts are transverse projection proxies, not optical-coupling
or prototype simulations. Projection inputs use the legacy rod-centred frame;
the raw reader uses the explicitly pinned physical world-to-calo transform.
"""
from __future__ import annotations
from dataclasses import dataclass
from enum import IntEnum
import numpy as np
from .depth import C_MM_NS, DepthContract

FIBER_OFFSETS_MM = np.array([[0., 0.], [.81*np.cos(np.pi/6), .405],
                             [-.81*np.cos(np.pi/6), .405], [0., -.81]])


def rotation(x_deg=2., y_deg=2.):
    if not np.isfinite([x_deg, y_deg]).all():
        raise ValueError("Nonfinite rotation")
    x, y = np.deg2rad([x_deg, y_deg])
    cx, sx, cy, sy = np.cos(x), np.sin(x), np.cos(y), np.sin(y)
    return np.array([[cy, 0., sy], [0., 1., 0.], [-sy, 0., cy]]) @ np.array([[1., 0., 0.], [0., cx, -sx], [0., sx, cx]])


def world_to_calorimeter(points_mm, *, x_deg=2., y_deg=2.):
    p = np.asarray(points_mm, float)
    if p.ndim != 2 or p.shape[1] != 3 or not np.isfinite(p).all():
        raise ValueError("Need finite N x 3 world coordinates in mm")
    return p @ rotation(x_deg, y_deg).T


def localize_endpoints(points_world_mm, fiber, *, x_deg=2., y_deg=2.):
    calo = world_to_calorimeter(points_world_mm, x_deg=x_deg, y_deg=y_deg)
    f = np.asarray(fiber)
    if f.shape != (len(calo),) or f.dtype.kind not in "iu" or np.any((f < 0) | (f > 3)):
        raise ValueError("Invalid fiber identity")
    offset = FIBER_OFFSETS_MM[f]
    rod = 4*np.rint((calo[:, :2]-offset-2)/4)+2
    local = calo[:, :2]-offset-rod
    index = np.rint((rod-2)/4).astype(np.int64)
    if (np.any(abs(calo[:, 2]-1000) > 1e-4) or np.any(np.linalg.norm(local, axis=1) > .40001) or
            np.any(index[:, 0] < -45) or np.any(index[:, 0] > 44) or
            np.any(index[:, 1] < -40) or np.any(index[:, 1] > 39)):
        raise ValueError("Downstream plane, fiber radius or rod limits mismatch")
    ids = ((index[:, 0]+2048) << 16) | ((index[:, 1]+2048) << 2) | f
    return local, rod+offset, index, ids, calo


def dreamsim_front_contract(*, gun_world_mm=(25., -25., -1050.), x_deg=2., y_deg=2., **sensor_optics):
    """Prompt world-z beam front and exact entrance time for pinned geometry.

    front time uses SENSOR CENTRES in the physical calo frame, an explicitly
    labelled sub-fiber transverse approximation (<~0.1 ps here). A general
    off-axis/slow shower secondary need not lie on this beta=1 plane.
    Returns (contract, entrance_time_ns, entrance_world_mm).
    """
    r = rotation(x_deg, y_deg)
    gun = np.asarray(gun_world_mm, float)
    if gun.shape != (3,) or not np.isfinite(gun).all() or r[2, 2] <= 0:
        raise ValueError("Invalid world-z beam front")
    distance = (-1000.-(r @ gun)[2])/r[2, 2]
    if distance < 0:
        raise ValueError("Gun must precede entrance plane")
    face = gun+np.array([0., 0., distance])
    contract = DepthContract(front_z=float(r[2, 2]), front_x=float(r[0, 2]), front_y=float(r[1, 2]),
        front_offset_mm=float(-r[2, 2]*1000-face[2]),
        geometry_id=f"DREAMSim-54b2447-RyRx-{x_deg:g}-{y_deg:g}-worldZ-front-v1", **sensor_optics)
    if (contract.length_mm != 2000. or contract.core_radius_mm != .39 or
            contract.fiber_outer_radius_mm != .4):
        raise ValueError("Pinned DREAMSim geometry is L=2000 mm, core radius=0.39 mm, outer radius=0.40 mm")
    return contract, float(distance/C_MM_NS), face


ROD_PITCH_CM = 0.4


ORIGINAL_DREAM_FIBER_DIAMETER_MM = 0.8


HGDREAM_OUTER_BOARD_ORIGINS = (
    (-13.5, 3.5),
    (-9.5, 7.5),
    (-5.5, 7.5),
    (-1.5, 9.5),
    (2.5, 7.5),
    (6.5, 7.5),
    (10.5, 3.5),
    (-9.5, -0.5),
    (-5.5, -0.5),
    (-1.5, -2.5),
    (2.5, -0.5),
    (6.5, -0.5),
)


HGDREAM_CORE_BOARD_ORIGINS = ((-1.5, 1.875), (0.5, 1.875))


ORIGINAL_DREAM_TOWER_APEX_TO_APEX_MM = 80.0


ORIGINAL_DREAM_TOWER_RADIUS_MM = 0.5 * ORIGINAL_DREAM_TOWER_APEX_TO_APEX_MM


ORIGINAL_DREAM_TOWER_SPACING_MM = (
    np.sqrt(3.0) * ORIGINAL_DREAM_TOWER_RADIUS_MM
)


@dataclass(frozen=True)
class TransverseChannelLayout:
    """Rectangular calorimeter channels and their coupled sensor sizes."""

    center_x_mm: np.ndarray
    center_y_mm: np.ndarray
    width_x_mm: np.ndarray
    width_y_mm: np.ndarray
    sensor_width_mm: np.ndarray
    rods_per_channel: np.ndarray
    is_core: np.ndarray

    @property
    def n_channels(self) -> int:
        return int(self.center_x_mm.size)


def hgdream_cherenkov_channel_layout() -> TransverseChannelLayout:
    """Return the 448-channel Cherenkov-side HG-DREAM segmentation.

    The layout contains 384 outer channels that each collect 12 rods onto a
    6 x 6 mm2 SiPM and 64 core channels that each collect three rods onto a
    3 x 3 mm2 SiPM.  Board origins are copied from the current test-beam
    channel map; the positions are converted to millimetres at the absorber
    face.
    """

    rows: list[tuple[float, float, float, float, float, int, bool]] = []
    for origin_x, origin_y in HGDREAM_OUTER_BOARD_ORIGINS:
        for local_x in range(4):
            for local_y in range(8):
                rows.append(
                    (
                        12.0 * (origin_x + local_x),
                        16.0 * (origin_y - local_y),
                        12.0,
                        16.0,
                        6.0,
                        12,
                        False,
                    )
                )
    for origin_x, origin_y in HGDREAM_CORE_BOARD_ORIGINS:
        for local_x in range(2):
            for local_y in range(16):
                rows.append(
                    (
                        12.0 * (origin_x + local_x),
                        16.0 * (origin_y - local_y / 4.0),
                        12.0,
                        4.0,
                        3.0,
                        3,
                        True,
                    )
                )
    values = np.asarray(rows, dtype=np.float64)
    return TransverseChannelLayout(
        center_x_mm=values[:, 0],
        center_y_mm=values[:, 1],
        width_x_mm=values[:, 2],
        width_y_mm=values[:, 3],
        sensor_width_mm=values[:, 4],
        rods_per_channel=values[:, 5].astype(np.int64),
        is_core=values[:, 6].astype(bool),
    )


def hgdream_cherenkov_channel_ids(
    event: object,
    layout: TransverseChannelLayout | None = None,
) -> np.ndarray:
    """Map each selected photon to the current HG-DREAM C channel."""

    layout = hgdream_cherenkov_channel_layout() if layout is None else layout
    photon_x_mm = 10.0 * ROD_PITCH_CM * event.rod_ix
    photon_y_mm = 10.0 * ROD_PITCH_CM * event.rod_iy
    channel_ids = np.full(event.time_ns.size, -1, dtype=np.int64)
    for channel_id in range(layout.n_channels):
        half_x = 0.5 * layout.width_x_mm[channel_id]
        half_y = 0.5 * layout.width_y_mm[channel_id]
        inside = (
            (photon_x_mm >= layout.center_x_mm[channel_id] - half_x)
            & (photon_x_mm < layout.center_x_mm[channel_id] + half_x)
            & (photon_y_mm >= layout.center_y_mm[channel_id] - half_y)
            & (photon_y_mm < layout.center_y_mm[channel_id] + half_y)
        )
        if np.any(inside & (channel_ids >= 0)):
            raise ValueError("HG-DREAM channel rectangles overlap")
        channel_ids[inside] = channel_id
    return channel_ids


def original_dream_pmt_centers(
    *, spacing_mm: float = ORIGINAL_DREAM_TOWER_SPACING_MM,
) -> np.ndarray:
    """Return the 19 tower centers of the original DREAM PMT readout."""

    if spacing_mm <= 0:
        raise ValueError("PMT tower spacing must be positive")
    axial_coordinates = [
        (q, r)
        for q in range(-2, 3)
        for r in range(-2, 3)
        if max(abs(q), abs(r), abs(q + r)) <= 2
    ]
    return np.asarray(
        [
            (
                spacing_mm * (q + 0.5 * r),
                spacing_mm * np.sqrt(3.0) * 0.5 * r,
            )
            for q, r in axial_coordinates
        ],
        dtype=np.float64,
    )


def original_dream_pmt_channel_ids(
    event: object,
    centers_mm: np.ndarray | None = None,
    *,
    tower_radius_mm: float = ORIGINAL_DREAM_TOWER_RADIUS_MM,
) -> np.ndarray:
    """Project the full-detector photons onto the original 19 C PMTs.

    The available event was generated with the HG-DREAM geometry, so this is
    a transverse nearest-tower projection rather than a separate simulation
    of the original DREAM module.
    """

    centers = original_dream_pmt_centers() if centers_mm is None else np.asarray(
        centers_mm, dtype=np.float64
    )
    if centers.shape != (19, 2):
        raise ValueError("original DREAM layout must contain 19 x-y centers")
    if tower_radius_mm <= 0:
        raise ValueError("original DREAM tower radius must be positive")
    photon_xy = np.column_stack(
        (
            10.0 * ROD_PITCH_CM * event.rod_ix,
            10.0 * ROD_PITCH_CM * event.rod_iy,
        )
    )
    distance_squared = np.sum(
        (photon_xy[:, np.newaxis, :] - centers[np.newaxis, :, :]) ** 2,
        axis=2,
    )
    channel_ids = np.argmin(distance_squared, axis=1).astype(np.int64)

    # Keep only photons inside the nearest pointy-top regular hexagon.  This
    # prevents an edge tower from absorbing photons that fall outside the
    # physical 19-tower DREAM footprint.
    local_xy = photon_xy - centers[channel_ids]
    abs_x = np.abs(local_xy[:, 0])
    abs_y = np.abs(local_xy[:, 1])
    inside = (
        (abs_x <= 0.5 * np.sqrt(3.0) * tower_radius_mm)
        & (abs_y <= tower_radius_mm - abs_x / np.sqrt(3.0))
    )
    channel_ids[~inside] = -1
    return channel_ids


# Ideal-cylinder diagnostic only; does not replace stored Geant4 photons.

class RayStatus(IntEnum):
    GUIDED_FACE = 0
    DIRECT_FACE = 1
    ESCAPED_WALL = 2
    BACKWARD = 3


@dataclass(frozen=True)
class RayResult:
    exit_xy_mm: np.ndarray
    arrival_ns: np.ndarray
    status: np.ndarray
    reflections: np.ndarray


def trace_cylinder(position_mm, direction, *, length_mm=2000., radius_mm=.4,
                   n_phase=1.49, n_clad=1.42, n_group=1.49, emission_ns=0.):
    p, d = np.asarray(position_mm, float), np.asarray(direction, float)
    if (p.ndim != 2 or p.shape[1] != 3 or d.shape != p.shape or
            not np.isfinite(p).all() or not np.isfinite(d).all()):
        raise ValueError("Need finite aligned N x 3 positions/directions")
    if (not np.isfinite([length_mm, radius_mm, n_phase, n_clad, n_group]).all() or
            not (length_mm > 0 and radius_mm > 0 and n_phase > n_clad > 0 and n_group > 0)):
        raise ValueError("Invalid optical model")
    if np.any(abs(np.linalg.norm(d, axis=1)-1) > 1e-10):
        raise ValueError("Directions must be unit vectors")
    if np.any(np.linalg.norm(p[:, :2], axis=1) >= radius_mm) or np.any((p[:, 2] < 0) | (p[:, 2] > length_mm)):
        raise ValueError("Source must be inside cylinder and within its axial length")
    emit = np.broadcast_to(np.asarray(emission_ns, float), (len(p),))
    if not np.isfinite(emit).all():
        raise ValueError("Invalid emission time")
    n = len(p)
    xy, time = np.full((n, 2), np.nan), np.full(n, np.nan)
    count = np.zeros(n, dtype=np.int64)
    status = np.full(n, int(RayStatus.BACKWARD), dtype=np.int8)
    forward = np.flatnonzero(d[:, 2] > 0)
    for j in forward:
        dz, pt = d[j, 2], np.linalg.norm(d[j, :2])
        remaining = length_mm-p[j, 2]
        if pt < 1e-15:
            xy[j], status[j] = p[j, :2], RayStatus.DIRECT_FACE
        else:
            u = d[j, :2]/pt
            pos = p[j, :2]
            dot = float(pos @ u)
            wall = -dot+np.sqrt(max(0., dot*dot+radius_mm**2-float(pos @ pos)))
            distance = remaining*pt/dz
            if distance <= wall:
                xy[j], status[j] = pos+distance*u, RayStatus.DIRECT_FACE
            else:
                b = pos[0]*u[1]-pos[1]*u[0]
                guided = np.sqrt(dz*dz+(pt*b/radius_mm)**2) >= n_clad/n_phase
                if not guided:
                    status[j] = RayStatus.ESCAPED_WALL
                    continue
                hit = pos+wall*u
                reflected = u-2*float(u @ hit)*hit/radius_mm**2
                chord = 2*np.sqrt(max(0., radius_mm**2-b*b))
                # At zero impact the next wall is OPPOSITE, not unchanged.
                angle = (1. if b >= 0 else -1.)*2*np.arccos(np.clip(abs(b)/radius_mm, 0, 1))
                rem = distance-wall
                turns = int(np.floor(rem/chord))
                residual = rem-turns*chord
                theta = np.remainder(turns*angle, 2*np.pi)
                rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
                xy[j] = rot @ (hit+residual*reflected)
                count[j], status[j] = turns+1, RayStatus.GUIDED_FACE
        time[j] = emit[j]+remaining*n_group/(C_MM_NS*dz)
    return RayResult(xy, time, status, count)
