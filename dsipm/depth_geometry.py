"""Explicit DREAMSim 54b2447 geometry adapter, independent of legacy maps.

Port of the independently checked 2026-10-15 correction. Geant4 passive
placement convention: world-to-calo is Ry @ Rx. All public lengths are mm.
Do not silently call the legacy full_detector coordinate transformation.
"""
import numpy as np

from .depth_reconstruction import C_MM_NS, DepthContract

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
    if contract.length_mm != 2000. or contract.core_radius_mm != .4:
        raise ValueError("Pinned DREAMSim geometry is L=2000 mm, core radius=0.4 mm")
    return contract, float(distance/C_MM_NS), face
