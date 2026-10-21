"""Ideal straight-cylinder ray reference, NOT a replacement Geant4 transport.

No absorption, scattering, dispersion, coupling, end refraction or detection.
Photons are generated INSIDE the core. Direct face arrivals do not encounter
a wall and therefore do not require TIR. Group index sets time, phase indices
set wall guidance. This library does not generate new events.
"""
from dataclasses import dataclass
from enum import IntEnum

import numpy as np

from .depth_reconstruction import C_MM_NS


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
