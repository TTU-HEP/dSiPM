"""Corrected detector coordinates, readout projections and ideal-fiber reference.

HG-DREAM/DREAM layouts are transverse projection proxies, not optical-coupling
or prototype simulations. Projection inputs use the legacy rod-centred frame;
the raw reader uses the explicitly pinned physical world-to-calo transform.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass
from enum import IntEnum
import hashlib
from pathlib import Path
import re
import subprocess
import numpy as np
from .depth import C_MM_NS, DepthContract, digest

FIBER_OFFSETS_MM = np.array([[0., 0.], [.81*np.cos(np.pi/6), .405],
                             [-.81*np.cos(np.pi/6), .405], [0., -.81]])


@dataclass(frozen=True)
class DreamsimDetector:
    """Source-derived upstream detector, not a second independently tuned model.

    Only the current constant-index, seven-fiber replica geometry and fixed
    world-z gun are supported. Unsupported source forms fail closed. A new
    source revision, material, placement, beam or sampling scope changes the
    contract identity. This is software provenance, not transport certification.
    """
    revision: str
    source_hashes: tuple
    length_mm: float
    rod_pitch_mm: float
    n_rods: int
    n_layers: int
    hole_radius_mm: float
    core_radius_mm: float
    outer_radius_mm: float
    fibers: tuple  # (plastic/quartz/scintillation, copy number, x_mm, y_mm)
    indices: tuple  # (material, constant core index, constant cladding index)
    rotation_deg: tuple
    gun_world_mm: tuple
    sample_rod_layer: tuple

    def __post_init__(self):
        if not re.fullmatch('[0-9a-f]{40}', self.revision):
            raise ValueError('Need the full DREAMSim source revision')
        if (not self.source_hashes or len(dict(self.source_hashes)) != len(self.source_hashes)
                or any(not re.fullmatch('[0-9a-f]{64}', h) for _, h in self.source_hashes)):
            raise ValueError('Need unique source-file hashes')
        numeric = [self.length_mm, self.rod_pitch_mm, self.hole_radius_mm,
                   self.core_radius_mm, self.outer_radius_mm, *self.rotation_deg, *self.gun_world_mm]
        scope = self.sample_rod_layer
        valid_scope = (len(scope) == 2 and all(isinstance(n, int) for n in scope) and
                       all(isinstance(n, int) for n in (self.n_rods, self.n_layers)) and
                       (tuple(scope) == (-1, -1) or
                        (0 <= scope[0] < self.n_rods and 0 <= scope[1] < self.n_layers)))
        if (not np.isfinite(numeric).all() or min(numeric[:5]) <= 0 or
                not 0 < self.core_radius_mm <= self.outer_radius_mm < self.hole_radius_mm or
                any(not isinstance(n, int) or n <= 0 or n % 2 for n in (self.n_rods, self.n_layers)) or
                len(self.rotation_deg) != 2 or len(self.gun_world_mm) != 3 or
                not valid_scope):
            raise ValueError('Invalid upstream detector dimensions or sampling scope')
        if len(self.fibers) != 7 or len({(f[0], f[1]) for f in self.fibers}) != 7:
            raise ValueError('Need seven distinct material/copy identities')
        expected = {'plastic': {0, 1}, 'quartz': {2, 3}, 'scintillation': {1, 2, 3}}
        if any({f[1] for f in self.fibers if f[0] == k} != v for k, v in expected.items()):
            raise ValueError('Unsupported DREAMSim fiber material/copy layout')
        xy = np.asarray([f[2:] for f in self.fibers], float)
        if not np.isfinite(xy).all() or np.any(np.linalg.norm(xy, axis=1)+self.outer_radius_mm > self.hole_radius_mm+1e-9):
            raise ValueError('Fiber outside source-defined rod hole')
        distance = np.linalg.norm(xy[:, None]-xy[None, :], axis=2)
        if np.any(distance[np.triu_indices(7, 1)] < 2*self.outer_radius_mm-1e-9):
            raise ValueError('Overlapping DREAMSim fiber placements; fix the producer, not the reader')
        if {i[0] for i in self.indices} != set(expected) or len(self.indices) != 3:
            raise ValueError('Missing material optics')
        if any(not np.isfinite([n, c]).all() or not n > c > 0 for _, n, c in self.indices):
            raise ValueError('Invalid source material indices')

    @property
    def sha256(self):
        return digest(asdict(self))

    @property
    def tracks_all_rods(self):
        """Explicit producer opt-in, never inferred from available birth records."""
        return tuple(self.sample_rod_layer) == (-1, -1)

    @classmethod
    def from_checkout(cls, checkout, macro=None):
        """Read a clean committed DREAMSim tree and the actual run-parameter macro.

        This does not compile or execute C++. A macro outside the checkout is
        allowed as a run input; the caller records its separate full hash.
        """
        root = Path(checkout)
        def git(*args):
            return subprocess.check_output(['git', '-C', str(root), *args], text=True).strip()
        if git('status', '--porcelain', '--untracked-files=no'):
            raise ValueError('Commit DREAMSim source changes before deriving a detector')
        paths = ['sim/src/CaloXDetectorConstruction.cc', 'sim/src/CaloXSteppingAction.cc',
                 'sim/src/CaloXTree.cc', 'sim/src/CaloXPrimaryGeneratorAction.cc']
        sources = {p: (root/p).read_text() for p in paths}
        run_macro = root/'sim/paramBatch03_single.mac' if macro is None else Path(macro)
        return cls.from_sources(git('rev-parse', 'HEAD'), sources, run_macro.read_text())

    @classmethod
    def from_sources(cls, revision, sources, macro):
        """Parse the supported source form; fixtures can exercise the same parser.

        No eval, C++ execution, guessed defaults, or dispersion-to-constant fit.
        A changed source form requires an explicit reader review.
        """
        clean = lambda s: re.sub(r'/\*.*?\*/|//[^\n]*', '', s, flags=re.S)
        g = clean(sources['sim/src/CaloXDetectorConstruction.cc'])
        step = clean(sources['sim/src/CaloXSteppingAction.cc'])
        def one(pattern, text=g):
            matches = re.findall(pattern, text, re.S)
            if len(matches) != 1:
                raise ValueError('Unsupported/ambiguous DREAMSim source: '+pattern)
            return matches[0]
        number = r'([0-9]+(?:\.[0-9]*)?|\.[0-9]+)'
        def scalar(name, unit=None):
            return float(one(r'double\s+'+name+r'\s*=\s*'+number+
                             (r'\s*\*\s*'+unit if unit else '')+r'\s*;'))
        length, pitch = 10*scalar('fiberLength', 'cm'), 10*scalar('rodSize', 'cm')
        nr, nl = scalar('noRods'), scalar('noLayers')
        if nr != int(nr) or nl != int(nl):
            raise ValueError('Nonintegral replica count')
        radius = scalar('clad_Plastic_rMax', 'mm')
        core = scalar('core_Plastic_rMax', 'mm')
        if scalar('clad_S_rMax', 'mm') != radius or scalar('core_S_rMax', 'mm') != core:
            raise ValueError('Unsupported unequal source fiber radii')
        if 'fiberCoreQuartzLog = new G4LogicalVolume(fiberC,' not in g:
            raise ValueError('Quartz must use the source-defined Cherenkov solid')
        spacing = radius*float(one(r'double\s+R\s*=\s*clad_Plastic_rMax\s*\*\s*'+number+r'\s*\+\s*[0-9.]+\s*;'))
        spacing += float(one(r'double\s+R\s*=\s*clad_Plastic_rMax\s*\*\s*[0-9.]+\s*\+\s*'+number+r'\s*;'))
        if not (re.search(r'double\s+cx1\s*=\s*R\s*\*\s*cos\(30\.0\s*\*\s*deg\)', g) and
                re.search(r'double\s+cy1\s*=\s*R\s*\*\s*sin\(30\.0\s*\*\s*deg\)', g)):
            raise ValueError('Unsupported fiber-ring definition')
        symbols = dict(R=spacing, cx1=spacing*np.cos(np.pi/6), cy1=spacing*.5)
        def coordinate(token):
            token = token.strip()
            sign = -1 if token.startswith('-') else 1
            token = token.lstrip('+-').strip()
            if token in symbols:
                return float(sign*symbols[token])
            if re.fullmatch(r'[0-9]+\.?[0-9]*', token):
                return sign*float(token)
            raise ValueError('Unsupported fiber coordinate: '+token)
        fibers = []
        for xyz, logical, copy in re.findall(
                r'new G4PVPlacement\(0,\s*G4ThreeVector\(([^)]+)\),\s*(fiberPlasticLog|fiberQuartzLog|fiberSLog),'
                r'\s*"[^"]+",\s*holeLV,\s*false,\s*(\d+),\s*fCheckOverlaps\)', g):
            x, y, z = [coordinate(t) for t in xyz.split(',')]
            if z != 0:
                raise ValueError('Unsupported longitudinal fiber offset')
            kind = {'fiberPlasticLog':'plastic', 'fiberQuartzLog':'quartz', 'fiberSLog':'scintillation'}[logical]
            fibers.append((kind, int(copy), x, y))
        def index(name):
            match = re.search(r'G4double\s+RefractiveIndex_'+name+r'\[nEntries\]\s*=\s*\{([^}]+)\}', g)
            if match:
                values = [float(v.strip()) for v in match[1].split(',') if v.strip()]
                if not values or len(set(values)) != 1:
                    raise ValueError('Dispersive material needs an explicit group-index model')
                return values[0]
            return float(one(r'std::fill_n\(RefractiveIndex_'+name+r',\s*nEntries,\s*'+number+r'\)'))
        optics = (('plastic', index('PMMA'), index('FluorinatedPolymer')),
                  ('quartz', index('FusedSilica'), index('HardPolymer')),
                  ('scintillation', index('Polystyrene'), index('PMMA')))
        params = dict(re.findall(r'^#\$\$\$\s+(\w+)\s+(\S+)', macro, re.M))
        def parameter(name):
            if name not in params:
                raise ValueError('Run macro must explicitly define '+name)
            return float(params[name])
        gun = []
        for axis in 'xyz':
            lo, hi = parameter('gun_'+axis+'_min'), parameter('gun_'+axis+'_max')
            if lo != hi or parameter('pMomentum_'+axis+'_min') != (1 if axis == 'z' else 0) or parameter('pMomentum_'+axis+'_max') != (1 if axis == 'z' else 0):
                raise ValueError('Only a fixed world-z beam front is implemented')
            gun.append(10*lo)
        sample = [parameter('opSampleRod'), parameter('opSampleLayer')]
        if any(v != int(v) for v in sample):
            raise ValueError('Nonintegral sampling rod/layer')
        policy = re.sub(r'\s+', '', one(r'\bconst\s+bool\s+isSampleRod\s*=\s*([^;]+);', text=step))
        single = '(rodNumber==opSampleRod&&layerNumber==opSampleLayer)'
        full = '(opSampleRod==-1&&opSampleLayer==-1)||'+single
        if policy not in (single, full):
            raise ValueError('Unsupported upstream explicit-optical scope')
        if sample == [-1, -1] and policy != full:
            raise ValueError('This producer does not implement whole-detector optical tracking')
        return cls(revision, tuple(sorted((p, hashlib.sha256(s.encode()).hexdigest()) for p, s in sources.items())),
            length, pitch, int(nr), int(nl), 5*scalar('holeDiameter', 'cm'), core, radius,
            tuple(fibers), optics, (parameter('caloRotationX'), parameter('caloRotationY')),
            tuple(gun), tuple(int(v) for v in sample))

    def contract(self, material='plastic', **sensor_options):
        return dreamsim_front_contract(detector=self, material=material, **sensor_options)

    def check_contract(self, contract, material):
        expected, _, _ = self.contract(material)
        for name in ('geometry_id', 'length_mm', 'core_radius_mm', 'fiber_outer_radius_mm',
                     'n_phase', 'n_group', 'n_clad', 'front_x', 'front_y', 'front_z', 'front_offset_mm'):
            if getattr(contract, name) != getattr(expected, name):
                raise ValueError('DREAMSim source/material contract mismatch: '+name)


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


def localize_endpoints(points_world_mm, fiber, *, x_deg=2., y_deg=2., detector=None, material='plastic'):
    """Historical CPAD geometry by default; upstream requires its source profile."""
    if detector is not None:
        x_deg, y_deg = detector.rotation_deg
    calo = world_to_calorimeter(points_world_mm, x_deg=x_deg, y_deg=y_deg)
    f = np.asarray(fiber)
    if f.shape != (len(calo),) or f.dtype.kind not in "iu" or np.any((f < 0) | (f > 3)):
        raise ValueError("Invalid fiber identity")
    if detector is None:
        offsets, pitch, halfz, radius, nx, ny = FIBER_OFFSETS_MM, 4., 1000., .4, 90, 80
    else:
        by_id = {i: (x, y) for kind, i, x, y in detector.fibers if kind == material}
        if material not in ('plastic', 'quartz') or not np.isin(f, list(by_id)).all():
            raise ValueError('Fiber identity does not belong to selected material')
        offsets = np.zeros((4, 2))
        for i, xy in by_id.items(): offsets[i] = xy
        pitch, halfz, radius = detector.rod_pitch_mm, detector.length_mm/2, detector.outer_radius_mm
        nx, ny = detector.n_rods, detector.n_layers
    offset = offsets[f]
    rod = pitch*np.rint((calo[:, :2]-offset-pitch/2)/pitch)+pitch/2
    local = calo[:, :2]-offset-rod
    index = np.rint((rod-pitch/2)/pitch).astype(np.int64)
    if (np.any(abs(calo[:, 2]-halfz) > 1e-4) or np.any(np.linalg.norm(local, axis=1) > radius+1e-5) or
            np.any(index[:, 0] < -nx//2) or np.any(index[:, 0] >= nx//2) or
            np.any(index[:, 1] < -ny//2) or np.any(index[:, 1] >= ny//2)):
        raise ValueError("Downstream plane, fiber radius or rod limits mismatch")
    ids = ((index[:, 0]+2048) << 16) | ((index[:, 1]+2048) << 2) | f
    return local, rod+offset, index, ids, calo


def dreamsim_front_contract(*, gun_world_mm=None, x_deg=2., y_deg=2., detector=None, material='plastic', **sensor_optics):
    """Prompt world-z beam front and exact entrance time for pinned geometry.

    front time uses SENSOR CENTRES in the physical calo frame, an explicitly
    labelled sub-fiber transverse approximation (<~0.1 ps here). A general
    off-axis/slow shower secondary need not lie on this beta=1 plane.
    Returns (contract, entrance_time_ns, entrance_world_mm).
    """
    halfz = 1000.
    if detector is not None:
        if gun_world_mm is not None:
            raise ValueError('The run macro, not the reader, defines the upstream beam')
        if material not in ('plastic', 'quartz'):
            raise ValueError('Depth inference currently supports plastic or quartz Cherenkov cores only')
        x_deg, y_deg = detector.rotation_deg
        gun_world_mm, halfz = detector.gun_world_mm, detector.length_mm/2
        _, ncore, nclad = next(i for i in detector.indices if i[0] == material)
        values = dict(length_mm=detector.length_mm, core_radius_mm=detector.core_radius_mm,
            fiber_outer_radius_mm=detector.outer_radius_mm, n_phase=ncore, n_group=ncore, n_clad=nclad)
        if any(k in sensor_optics and sensor_optics[k] != v for k, v in values.items()):
            raise ValueError('DREAMSim defines geometry/optics; do not override them in dSiPM')
        sensor_optics.update(values)
    elif gun_world_mm is None:
        gun_world_mm = (25., -25., -1050.)
    r = rotation(x_deg, y_deg)
    gun = np.asarray(gun_world_mm, float)
    if gun.shape != (3,) or not np.isfinite(gun).all() or r[2, 2] <= 0:
        raise ValueError("Invalid world-z beam front")
    distance = (-halfz-(r @ gun)[2])/r[2, 2]
    if distance < 0:
        raise ValueError("Gun must precede entrance plane")
    face = gun+np.array([0., 0., distance])
    contract = DepthContract(front_z=float(r[2, 2]), front_x=float(r[0, 2]), front_y=float(r[1, 2]),
        front_offset_mm=float(-r[2, 2]*halfz-face[2]),
        geometry_id=(f'DREAMSim-source-{detector.sha256}-{material}-worldZ-front-v2' if detector else
                     f"DREAMSim-54b2447-RyRx-{x_deg:g}-{y_deg:g}-worldZ-front-v1"), **sensor_optics)
    if detector is None and (contract.length_mm != 2000. or contract.core_radius_mm != .39 or
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
