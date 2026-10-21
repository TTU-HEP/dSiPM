"""Observable-only, single-ended fiber-depth inference (development model).

Truth is accepted only by ``fit_calibration`` and the separate evaluation functions.
No estimator here is a validated detector resolution or a transport-bug fix.
Units throughout are millimetres and nanoseconds; z=0 is the entrance face.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import IntEnum
import hashlib
import json
from pathlib import Path

import numpy as np

C_MM_NS = 299.792458
ALGORITHM_VERSION = "fiber-depth-v2-core-radius"


def digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
                                     separators=(",", ":")).encode()).hexdigest()


@dataclass(frozen=True)
class DepthContract:
    """Compatibility contract; these defaults are assumptions, not measurements.

    front_path(z,x,y) = front_z*z + front_x*x + front_y*y + front_offset.
    This is the prompt beta=1 production front relative to the entrance clock.
    phase/group indices are deliberately independent. gate_ns is on the stored
    absolute arrival time, BEFORE instrument jitter and entrance subtraction.
    """
    length_mm: float = 2000.0
    core_radius_mm: float = 0.39
    fiber_outer_radius_mm: float = 0.4
    n_phase: float = 1.49
    n_clad: float = 1.42
    n_group: float = 1.49
    front_z: float = 1.0
    front_x: float = 0.0
    front_y: float = 0.0
    front_offset_mm: float = 0.0
    geometry_id: str = "ideal-axial-fiber-v1"
    pitch_um: float = 20.0
    pixel_origin_mm: float = -1.5
    pde: float = 0.2
    recovery_ns: float = 4.3
    gate_ns: tuple[float, float] = (0.0, 40.0)
    sensor_model: str = "independent-fiber-nonparalyzable-v1"
    workflow_sha256: str = ""  # Empty: standalone model, not canonical provenance.

    def __post_init__(self):
        vals = [v for v in asdict(self).values() if isinstance(v, (float, int))]
        if not np.isfinite(vals).all():
            raise ValueError("Nonfinite contract")
        if not (self.length_mm > 0 and 0 < self.core_radius_mm <= self.fiber_outer_radius_mm and
                self.n_phase > self.n_clad > 0 and self.n_group > abs(self.front_z)):
            raise ValueError("Invalid optical/front model")
        if not (self.pitch_um > 0 and 0 <= self.pde <= 1 and self.recovery_ns >= 0):
            raise ValueError("Invalid sensor model")
        if len(self.gate_ns) != 2 or not np.isfinite(self.gate_ns).all() or self.gate_ns[0] >= self.gate_ns[1]:
            raise ValueError("Invalid stored-time gate")
        if not self.geometry_id or self.sensor_model != "independent-fiber-nonparalyzable-v1":
            raise ValueError("Unspecified geometry or unsupported sensor model")
        if self.workflow_sha256 and (len(self.workflow_sha256) != 64 or
                                     any(c not in "0123456789abcdef" for c in self.workflow_sha256)):
            raise ValueError("Malformed sensor workflow SHA-256")

    @classmethod
    def from_workflow(cls, workflow=None, **overrides):
        """Read canonical PDE/pitch/holdoff and pin the complete JSON hash.

        workflow is a path (default: repository sensor_workflow.json). Shared
        timestamp export is a separate explicit sensor API, not automatic.
        """
        path = Path(workflow) if workflow is not None else Path(__file__).resolve().parents[1]/"sensor_workflow.json"
        raw = path.read_bytes()
        data = json.loads(raw)
        reference = data["reference_scenario"]
        if data["schema_version"] != 1 or reference["detection_model"] != "constant_effective_pde":
            raise ValueError("Unsupported sensor workflow")
        values = dict(pde=reference["effective_pde"], pitch_um=reference["spad_pitch_um"],
                      recovery_ns=reference["spad_recovery_ns"], workflow_sha256=hashlib.sha256(raw).hexdigest())
        geometry = data.get("geometry", {})
        for name in ("length_mm", "core_radius_mm", "fiber_outer_radius_mm", "n_phase", "n_clad", "n_group"):
            if name in geometry:
                values[name] = geometry[name]
        values.update(overrides)
        return cls(**values)

    @property
    def sha256(self):
        return digest(asdict(self))

    @property
    def radius_padding(self):
        # Half a square pixel's diagonal, in core-radius units.
        return self.pitch_um / 1000 / np.sqrt(2) / self.core_radius_mm

    def forward_time(self, z_mm, kappa, x_mm=0.0, y_mm=0.0, emission_delay_ns=0.0):
        """Forward arrival relative to entrance; emission delay is a diagnostic."""
        front = self.front_z * np.asarray(z_mm) + self.front_x*x_mm + self.front_y*y_mm + self.front_offset_mm
        return (front + self.n_group*np.asarray(kappa)*(self.length_mm-np.asarray(z_mm)))/C_MM_NS + emission_delay_ns

    def inverse(self, time_ns, kappa, x_mm=0.0, y_mm=0.0):
        nk = self.n_group * np.asarray(kappa)
        return (nk*self.length_mm + self.front_x*x_mm + self.front_y*y_mm + self.front_offset_mm
                - C_MM_NS*np.asarray(time_ns)) / (nk-self.front_z)

    def kappa_bound(self, radius):
        """Ideal-cylinder TIR bound; NOT a path measurement or loss correction."""
        r = np.minimum(np.asarray(radius, float), 1.0)
        sin2 = (1-(self.n_clad/self.n_phase)**2) / np.maximum(1-r*r, 1e-15)
        return np.where(sin2 < 1, 1/np.sqrt(np.maximum(1-sin2, 1e-15)), np.inf)


@dataclass(frozen=True)
class PhotonObservations:
    """Only measured quantities. event_key must identify source AND entry.

    photon_id is the original array index, not an index renumbered after cuts.
    Rod indices may be signed; fibers are separate for sensor holdoff, while
    neighbouring-light inference uses rods. Radius is a pixel-centre estimate.
    """
    event_key: np.ndarray
    photon_id: np.ndarray
    rod_x: np.ndarray
    rod_y: np.ndarray
    time_ns: np.ndarray
    radius: np.ndarray
    sensor_x_mm: np.ndarray
    sensor_y_mm: np.ndarray

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            value = np.asarray(getattr(self, name))
            if value.ndim != 1:
                raise ValueError(f"{name} must be one dimensional")
            if name == "event_key":
                value = value.astype(str)
                if np.any(value == ""):
                    raise ValueError("Empty event identity")
            elif name in ("photon_id", "rod_x", "rod_y"):
                if value.dtype.kind not in "iu":
                    raise ValueError(f"{name} must have an integer dtype")
                if value.dtype.kind == "u" and np.any(value > np.iinfo(np.int64).max):
                    raise ValueError("Integer identity overflow")
                value = value.astype(np.int64)
            else:
                value = value.astype(float)
                if not np.isfinite(value).all():
                    raise ValueError(f"Nonfinite {name}")
            value = value.copy()
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if len({len(getattr(self, k)) for k in self.__dataclass_fields__}) != 1:
            raise ValueError("Column length mismatch")
        if np.any(self.photon_id < 0) or np.any(self.radius < 0):
            raise ValueError("Negative photon identity or radius")
        if len(set(zip(self.event_key.tolist(), self.photon_id.tolist()))) != len(self.time_ns):
            raise ValueError("Duplicate event/photon identity")

    def subset(self, mask):
        return PhotonObservations(**{k: getattr(self, k)[mask] for k in self.__dataclass_fields__})


def channel_delay(obs: PhotonObservations):
    """Measured delay to first observed photon of the same event AND rod."""
    n = len(obs.time_ns)
    if not n:
        return np.empty(0)
    order = np.lexsort((obs.photon_id, obs.time_ns, obs.rod_y, obs.rod_x, obs.event_key))
    e, x, y = obs.event_key[order], obs.rod_x[order], obs.rod_y[order]
    first = np.r_[True, (e[1:] != e[:-1]) | (x[1:] != x[:-1]) | (y[1:] != y[:-1])]
    group = np.cumsum(first)-1
    result = np.empty(n)
    result[order] = obs.time_ns[order] - obs.time_ns[order][first][group]
    return result


@dataclass(frozen=True)
class Calibration:
    payload: dict

    def __post_init__(self):
        p = self.payload
        if p["algorithm"] != ALGORITHM_VERSION:
            raise ValueError("Calibration algorithm mismatch")
        DepthContract(**p["contract"])
        if p["contract_sha256"] != digest(p["contract"]):
            raise ValueError("Corrupt calibration contract")
        re, de, qs = [np.asarray(p[k], float) for k in ("radius_edges", "delay_edges_ns", "quantiles")]
        if (any(a.ndim != 1 or not np.isfinite(a).all() for a in (re, de, qs)) or
                len(re) < 2 or re[0] != 0 or re[-1] != 1 or np.any(np.diff(re) <= 0) or
                len(de) < 2 or de[0] != 0 or np.any(np.diff(de) <= 0) or
                len(qs) != 10 or not np.allclose(qs, np.linspace(.05, .95, 10), rtol=0, atol=1e-14)):
            raise ValueError("Malformed calibration bins/quantiles")
        nr, nt, nq = len(p["radius_edges"])-1, len(p["delay_edges_ns"])+1, len(p["quantiles"])
        for name, shape in [("radius_median", (nr,)), ("joint_median", (nr, nt)),
                            ("joint_quantiles", (nr, nt, nq))]:
            a = np.asarray(p[name])
            if a.shape != shape or not np.isfinite(a).all() or np.any(a < 1):
                raise ValueError(f"Malformed calibration {name}")
        if np.any(np.diff(np.asarray(p["joint_quantiles"]), axis=-1) < 0):
            raise ValueError("Unordered kappa quantiles")
        if not p["training_events"] or not all(isinstance(e, str) and e for e in p["training_events"]):
            raise ValueError("Empty or malformed calibration event identities")
        if p["min_count"] < 1 or not np.isfinite(p["training_sigma_ns"]) or p["training_sigma_ns"] < 0:
            raise ValueError("Malformed calibration metadata")
        if len(set(p["training_events"])) != len(p["training_events"]):
            raise ValueError("Duplicate calibration event")
        for name, shape in [("radius_counts", (nr,)), ("joint_counts", (nr, nt))]:
            a = np.asarray(p[name])
            if a.shape != shape or a.dtype.kind not in "iu" or np.any(a < 0) or int(a.sum()) != p["training_photons"]:
                raise ValueError("Calibration count closure mismatch")

    @property
    def sha256(self):
        return digest(self.payload)

    def save(self, path):
        """Immutable, content-checked artifact; never silently overwrite."""
        with Path(path).open("x") as f:
            json.dump({"sha256": self.sha256, "calibration": self.payload}, f, indent=2, allow_nan=False)
            f.write("\n")

    @classmethod
    def load(cls, path):
        p = json.loads(Path(path).read_text())
        if digest(p["calibration"]) != p["sha256"]:
            raise ValueError("Calibration checksum mismatch")
        return cls(p["calibration"])

    def bins(self, obs):
        ri = np.clip(np.searchsorted(self.payload["radius_edges"], obs.radius, side="right")-1,
                     0, len(self.payload["radius_edges"])-2)
        # Exact-zero delay has its OWN bin; (0, .01 ns), ... and overflow follow.
        dt = channel_delay(obs)
        ti = np.searchsorted(self.payload["delay_edges_ns"], dt, side="right")
        ti[dt == 0] = 0
        return ri, ti


def fit_calibration(obs, kappa_truth, contract, *, source_sha256, training_sigma_ns=0.0,
                    min_count=20, radius_edges=None, delay_edges_ns=None):
    """Fit truth-labelled kappa on CALIBRATION EVENTS ONLY.

    Equal-probability midpoint quantiles have equal initial mass. The true 50%
    median is stored separately (the 55% candidate is not called a median).
    Sparse bins fall back radius->global; counts/fallback levels are archived.
    source_sha256 identifies the immutable input manifest/cache, not a pathname.
    """
    k = np.asarray(kappa_truth, float)
    if k.shape != obs.time_ns.shape or not len(k) or not np.isfinite(k).all() or np.any(k < 1):
        raise ValueError("Need finite aligned calibration kappa >= 1")
    if len(source_sha256) != 64 or any(c not in "0123456789abcdef" for c in source_sha256):
        raise ValueError("Need SHA-256 of calibration input")
    if min_count < 1 or training_sigma_ns < 0 or not np.isfinite(training_sigma_ns):
        raise ValueError("Invalid calibration options")
    if np.any(obs.radius > 1+contract.radius_padding+1e-10):
        raise ValueError("Pixel radius incompatible with core/pitch")
    re = np.asarray(np.linspace(0, 1, 41) if radius_edges is None else radius_edges, float)
    de = np.asarray(np.logspace(-2, np.log10(60), 36) if delay_edges_ns is None else delay_edges_ns, float)
    if (re.ndim != 1 or len(re) < 2 or not np.isfinite(re).all() or re[0] != 0 or re[-1] != 1 or
            np.any(np.diff(re) <= 0) or de.ndim != 1 or not len(de) or
            not np.isfinite(de).all() or de[0] <= 0 or np.any(np.diff(de) <= 0)):
        raise ValueError("Invalid calibration bins")
    qs = np.linspace(.05, .95, 10)
    nr = len(re)-1
    # Add a zero edge to separate zero from arbitrarily small positive delays.
    edges = np.r_[0., de]
    nt = len(edges)+1
    ri = np.clip(np.searchsorted(re, obs.radius, side="right")-1, 0, nr-1)
    dt = channel_delay(obs)
    ti = np.searchsorted(edges, dt, side="right")
    ti[dt == 0] = 0
    med = np.full(nr, np.median(k)); qm = np.tile(np.quantile(k, qs), (nr, 1))
    rc = np.bincount(ri, minlength=nr)
    for r in range(nr):
        if rc[r] >= min_count:
            med[r], qm[r] = np.median(k[ri == r]), np.quantile(k[ri == r], qs)
    jm = np.repeat(med[:, None], nt, axis=1)
    jq = np.repeat(qm[:, None, :], nt, axis=1)
    counts = np.zeros((nr, nt), int)
    fallback = np.broadcast_to(np.where(rc >= min_count, 1, 2)[:, None], (nr, nt)).copy()
    for r, t in set(zip(ri.tolist(), ti.tolist())):
        mask = (ri == r) & (ti == t)
        counts[r, t] = np.sum(mask)
        if counts[r, t] >= min_count:
            jm[r, t], jq[r, t] = np.median(k[mask]), np.quantile(k[mask], qs)
            fallback[r, t] = 0
    return Calibration(dict(algorithm=ALGORITHM_VERSION, contract=asdict(contract),
        contract_sha256=contract.sha256, source_sha256=source_sha256,
        training_events=sorted(set(obs.event_key.tolist())), training_sigma_ns=training_sigma_ns,
        training_photons=len(k),
        min_count=min_count, radius_edges=re.tolist(), delay_edges_ns=edges.tolist(), quantiles=qs.tolist(),
        radius_median=med.tolist(), joint_median=jm.tolist(), joint_quantiles=jq.tolist(),
        radius_counts=rc.tolist(), joint_counts=counts.tolist(), fallback_level=fallback.tolist()))


class DepthStatus(IntEnum):
    SUPPORTED = 0
    CALIBRATED_FALLBACK = 1
    NO_PHYSICAL_CANDIDATE = 2


@dataclass(frozen=True)
class DepthResult:
    z_mm: np.ndarray
    status: np.ndarray
    candidate_z_mm: np.ndarray
    candidate_density: np.ndarray
    prior_z_mm: np.ndarray
    prior_seed: np.ndarray
    invalid_clean_seed_count: int
    calibration_sha256: str
    options: dict
    selected_kappa: np.ndarray
    photon_timing_sigma_z_mm: np.ndarray
    entrance_timing_sigma_z_mm: np.ndarray


def _kernel(sigma_bins):
    """One exact finite kernel, used for BOTH smoothing and self subtraction."""
    radius = int(4*sigma_bins + .5)
    offsets = np.arange(-radius, radius+1)
    weights = np.exp(-.5*(offsets/sigma_bins)**2)
    return weights / weights.sum()


def reconstruct(obs, calibration, contract, *, method="neighbour", sigma_t_ns=0.,
                entrance_sigma_ns=0., window=5, zbin_mm=20., clean_radius=.7,
                exclude_same_rod=False, allow_training=False, fixed_n_eff=None):
    """Fixed-calibration inference; no true z, kappa or emission time inputs.

    Compare radius / radius_delay / neighbour on the SAME accepted photons.
    Neighbour support is a heuristic density, NOT a calibrated posterior or
    uncertainty. It defaults to leave-one-photon-out; leave-one-rod-out is an
    explicit ablation. No physical candidate -> NaN plus failure status.
    Sparse/no-support fallback -> exact joint median, or nearest physical
    quantile in kappa if the median inverse is invalid (never a boundary clip).
    """
    p = calibration.payload
    if contract.sha256 != p["contract_sha256"]:
        raise ValueError("Calibration optical/geometry/sensor/gate contract mismatch")
    if not allow_training and set(obs.event_key).intersection(p["training_events"]):
        raise ValueError("Calibration/evaluation event overlap")
    if method not in ("fixed", "radius", "radius_delay", "neighbour"):
        raise ValueError("Unknown reconstruction method")
    if method == "fixed" and (fixed_n_eff is None or not np.isfinite(fixed_n_eff) or fixed_n_eff < contract.n_group):
        raise ValueError("Fixed-path baseline requires an explicit n_eff >= n_group")
    if (not isinstance(window, int) or window < 1 or window % 2 != 1 or
            not np.isfinite([sigma_t_ns, entrance_sigma_ns, zbin_mm, clean_radius]).all() or
            sigma_t_ns < 0 or entrance_sigma_ns < 0 or zbin_mm <= 0 or not 0 < clean_radius <= 1):
        raise ValueError("Invalid reconstruction options")
    if np.any(obs.radius > 1+contract.radius_padding+1e-10):
        raise ValueError("Radius outside pixel/core contract")
    n, length = len(obs.time_ns), contract.length_mm
    ri, ti = calibration.bins(obs)
    k1 = np.asarray(p["radius_median"])[ri]
    k2 = np.asarray(p["joint_median"])[ri, ti]
    bound = contract.kappa_bound(obs.radius+contract.radius_padding)
    kc = np.asarray(p["joint_quantiles"])[ri, ti]
    zc = contract.inverse(obs.time_ns[:, None], kc, obs.sensor_x_mm[:, None], obs.sensor_y_mm[:, None])
    zp = contract.inverse(obs.time_ns, k1, obs.sensor_x_mm, obs.sensor_y_mm)
    median_k = k1 if method == "radius" else k2
    if method == "fixed":
        median_k = np.full(n, fixed_n_eff/contract.n_group)
    selected_kappa = median_k.copy()
    z = contract.inverse(obs.time_ns, median_k, obs.sensor_x_mm, obs.sensor_y_mm)
    # Only sub-nanometre inversion roundoff is normalized at the endpoints.
    # This is not permission to clip any physically inconsistent estimate.
    physical = lambda v: np.isfinite(v) & (v >= -1e-9) & (v <= length+1e-9)
    def allowed(v, k, b):
        # A forward high-angle photon born right next to the face may arrive
        # before touching ANY wall. TIR is unnecessary for that population.
        # Its transverse path cannot exceed the longest chord ending at the
        # observed face radius (including pixel-centre uncertainty).
        transverse = np.maximum(length-v, 0)*np.sqrt(np.maximum(k*k-1, 0))
        chord = contract.core_radius_mm*(1+np.minimum(obs.radius+contract.radius_padding, 1))
        if v.ndim == 2:
            chord = chord[:, None]
        return physical(v) & ((k <= b+1e-12) | (transverse <= chord+1e-12))
    valid = allowed(zc, kc, bound[:, None])
    clean = obs.radius+contract.radius_padding < clean_radius
    seeds = clean & allowed(zp, k1, bound)  # NEVER clip invalid seeds to boundaries
    density = np.zeros_like(zc)
    status = np.full(n, int(DepthStatus.CALIBRATED_FALLBACK), dtype=np.int8)
    options = dict(method=method, sigma_t_ns=sigma_t_ns, entrance_sigma_ns=entrance_sigma_ns,
                   window=window, zbin_mm=zbin_mm, clean_radius=clean_radius,
                   exclude_same_rod=exclude_same_rod, allow_training=allow_training,
                   calibration_mode="fixed", fixed_n_eff=fixed_n_eff if method == "fixed" else None)
    if method == "neighbour":
        nz = int(np.ceil(length/zbin_mm))
        # Conservative front/timing scale; the entrance-clock term is common
        # within each event, not an independent photon fluctuation.
        width = max(30., np.hypot(sigma_t_ns, entrance_sigma_ns)*C_MM_NS/(contract.n_group-contract.front_z))
        kernel = _kernel(width/zbin_mm)
        half = len(kernel)//2
        for event in np.unique(obs.event_key):
            ids = np.flatnonzero(obs.event_key == event)
            hist = {}
            seed_bin = np.zeros(n, dtype=int)
            seed_ids = ids[seeds[ids]]
            seed_bin[seed_ids] = np.minimum((zp[seed_ids]/zbin_mm).astype(int), nz-1)
            for j in seed_ids:
                key = (int(obs.rod_x[j]), int(obs.rod_y[j]))
                if key not in hist:
                    hist[key] = np.zeros(nz)
                hist[key][seed_bin[j]] += 1
            smooth = {key: np.convolve(h, kernel, mode="full")[half:half+nz] for key, h in hist.items()}
            groups = {}
            for j in ids:
                groups.setdefault((int(obs.rod_x[j]), int(obs.rod_y[j])), []).append(j)
            w = window//2
            for (x, y), group in groups.items():
                js = np.asarray(group)
                h = np.zeros(nz)
                for dx in range(-w, w+1):
                    for dy in range(-w, w+1):
                        if exclude_same_rod and dx == dy == 0:
                            continue
                        h += smooth.get((x+dx, y+dy), 0.)
                bins = np.floor(np.clip(zc[js], 0, length)/zbin_mm).astype(int)
                bins = np.minimum(bins, nz-1)  # indexing only; validity is applied below
                d = h[bins]
                if not exclude_same_rod:
                    off = bins-seed_bin[js, None]
                    own = np.where(np.abs(off) <= half, kernel[np.clip(off+half, 0, len(kernel)-1)], 0.)
                    d -= seeds[js, None]*own
                # Tiny roundoff is NOT evidence for an isolated photon.
                d = np.where(d > 1e-12, d, 0.)
                density[js] = np.where(valid[js], d, 0.)
        physical_median = allowed(z, median_k, bound)
        no_median = ~physical_median & valid.any(axis=1)
        if np.any(no_median):
            distance = np.where(valid, abs(kc-median_k[:, None]), np.inf)
            j = np.argmin(distance[no_median], axis=1)
            z[no_median] = zc[no_median, j]
            selected_kappa[no_median] = kc[no_median, j]
        supported = density.max(axis=1) > 0 if n else np.empty(0, bool)
        if n:
            z[supported] = zc[supported, np.argmax(density[supported], axis=1)]
            selected_kappa[supported] = kc[supported, np.argmax(density[supported], axis=1)]
        status[supported] = DepthStatus.SUPPORTED
    # A supported/alternative candidate has already passed allowed(). A median
    # fallback must also satisfy its radius guidance OR direct-exit bound.
    failure = ~physical(z)
    if method not in ("neighbour", "fixed"):
        failure |= ~allowed(z, median_k, bound)
    elif method == "neighbour":
        failure |= (~physical_median & ~valid.any(axis=1))
    z[failure] = np.nan
    z[~failure] = np.clip(z[~failure], 0, length)
    status[failure] = DepthStatus.NO_PHYSICAL_CANDIDATE
    selected_kappa[failure] = np.nan
    derivative = C_MM_NS/(contract.n_group*selected_kappa-contract.front_z)
    # Conditional linear propagation only: NOT total depth errors/intervals.
    # It excludes path/emission ambiguity and candidate-selection uncertainty.
    return DepthResult(z, status, zc, density, zp, seeds, int(np.sum(clean & ~seeds)),
                       calibration.sha256, options, selected_kappa,
                       sigma_t_ns*derivative, entrance_sigma_ns*derivative)


# Truth-based evaluation. None of these truth arguments enter reconstruct().

def residual_summary(truth_mm, reconstructed_mm, *, tolerances_mm=(50., 100., 200.)):
    truth, z = np.asarray(truth_mm, float), np.asarray(reconstructed_mm, float)
    if truth.ndim != 1 or z.shape != truth.shape or not np.isfinite(truth).all():
        raise ValueError("Need finite aligned one-dimensional truth")
    if not np.isfinite(tolerances_mm).all() or np.any(np.asarray(tolerances_mm) <= 0):
        raise ValueError("Invalid tolerances")
    valid = np.isfinite(z)
    delta = z[valid]-truth[valid]
    q = np.percentile(delta, [16, 50, 84]) if len(delta) else [None]*3
    return dict(input=len(z), reconstructed=int(valid.sum()), unreconstructed=int((~valid).sum()),
        failure_fraction=float((~valid).mean()) if len(z) else None,
        conditional_bias_mm=None if not len(delta) else float(q[1]),
        conditional_68_halfwidth_mm=None if not len(delta) else float((q[2]-q[0])/2),
        inclusive_within={str(t): float(np.sum(abs(delta) <= t)/len(z)) if len(z) else None for t in tolerances_mm})


def evaluate(obs, truth_mm, result, contract, *, bootstrap_draws=200, seed=20261021,
             min_channel_photons=5, truth_edges_mm=None):
    truth = np.asarray(truth_mm, float)
    if truth.shape != obs.time_ns.shape or not np.isfinite(truth).all() or np.any((truth < 0) | (truth > contract.length_mm)):
        raise ValueError("Truth depth outside physical detector")
    if len(result.z_mm) != len(truth) or bootstrap_draws < 0 or min_channel_photons < 1:
        raise ValueError("Invalid evaluation options/lengths")
    if np.any(np.isfinite(result.z_mm) & ((result.z_mm < 0) | (result.z_mm > contract.length_mm))):
        raise ValueError("Unphysical depth incorrectly marked reconstructed")
    report = residual_summary(truth, result.z_mm)
    report["target"] = "creation depth of each accepted photon; failures retained"
    report["status_counts"] = {str(k): int(np.sum(result.status == k)) for k in [0, 1, 2]}
    if sum(report["status_counts"].values()) != len(truth):
        raise ValueError("Unknown reconstruction status")
    if not np.array_equal(result.status == 2, ~np.isfinite(result.z_mm)):
        raise ValueError("Failure-status/depth closure mismatch")
    report["excluded_invalid_prior_seeds"] = result.invalid_clean_seed_count
    edges = np.linspace(0, contract.length_mm, 11) if truth_edges_mm is None else np.asarray(truth_edges_mm, float)
    if not np.isfinite(edges).all() or len(edges) < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError("Invalid truth-depth bins")
    bins = np.searchsorted(edges, truth, side="right")-1
    bins[truth == edges[-1]] = len(edges)-2
    report["depth_bins"] = []
    for b in range(len(edges)-1):
        mask = bins == b
        report["depth_bins"].append(dict(low_mm=float(edges[b]), high_mm=float(edges[b+1]),
                                        **residual_summary(truth[mask], result.z_mm[mask])))
    report["truth_underflow"] = int(np.sum(truth < edges[0]))
    report["truth_overflow"] = int(np.sum(truth > edges[-1]))
    if sum(b["input"] for b in report["depth_bins"])+report["truth_underflow"]+report["truth_overflow"] != len(truth):
        raise ValueError("Truth-depth coverage failure")
    events = np.unique(obs.event_key)
    groups = {}
    for j, key in enumerate(zip(obs.event_key.tolist(), obs.rod_x.tolist(), obs.rod_y.tolist())):
        groups.setdefault(key, []).append(j)
    cm_truth, cm_reco = [], []
    for ids in groups.values():
        if len(ids) < min_channel_photons:
            continue
        ids = np.asarray(ids)
        rz = result.z_mm[ids]
        good = np.isfinite(rz)
        cm_truth.append(float(np.median(truth[ids])))
        cm_reco.append(float(np.median(rz[good])) if good.sum() >= min_channel_photons else np.nan)
    report["channel_medians"] = dict(target="median of accepted-photon creation depths; not pre-holdoff or energy truth",
        min_input_and_reconstructed_photons=min_channel_photons, **residual_summary(cm_truth, cm_reco))
    report["events"] = len(events)
    # Event-block bootstrap of PHOTON-WEIGHTED metrics; never resample individual photons.
    draws = []
    if len(events) >= 2 and bootstrap_draws:
        indices = [np.flatnonzero(obs.event_key == e) for e in events]
        rng = np.random.default_rng(seed)
        for _ in range(bootstrap_draws):
            ids = np.concatenate([indices[i] for i in rng.integers(0, len(events), len(events))])
            draws.append(residual_summary(truth[ids], result.z_mm[ids]))
    intervals = {}
    for metric in ("conditional_bias_mm", "conditional_68_halfwidth_mm", "failure_fraction"):
        values = [d[metric] for d in draws if d[metric] is not None]
        intervals[metric] = np.percentile(values, [2.5, 97.5]).tolist() if values else None
    for tolerance in ("50.0", "100.0", "200.0"):
        values = [d["inclusive_within"][tolerance] for d in draws]
        intervals["inclusive_within_"+tolerance+"mm"] = np.percentile(values, [2.5, 97.5]).tolist() if values else None
    report["event_block_bootstrap"] = dict(draws=len(draws), seed=seed, percentile_95=intervals,
        unavailable_reason="need at least two independent events" if len(events) < 2 else None)
    report["closure"] = report["input"] == report["reconstructed"]+report["unreconstructed"]
    return report
