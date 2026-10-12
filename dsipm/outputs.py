"""CSV, ROOT, and PNG outputs for device-comparison runs."""

from __future__ import annotations

from dataclasses import dataclass, field
import csv
import math
from pathlib import Path
from typing import Iterable

from .config import DeviceConfig
from .digitizer import DigitizedEvent
from .historical_outputs import HistoricalSlideOutputs


PLOT_STYLES = {"comparison", "slides", "both"}


EVENT_FIELDS = [
    "device",
    "event",
    "n_raw",
    "n_selected",
    "n_detected",
    "n_lost_pde",
    "n_accepted",
    "n_lost_deadtime",
    "n_active_spads",
    "n_active_timing_groups",
    "n_timing_records",
    "retained_fraction",
    "deadtime_loss_fraction",
    "mean_pde",
    "mean_wavelength_nm",
    "first_timestamp_ns",
    "last_timestamp_ns",
    "payload_bits_proxy",
]


@dataclass
class Accumulator:
    count: int = 0
    totals: dict[str, float] = field(default_factory=dict)
    accepted_squared: float = 0.0
    first_timestamp_count: int = 0
    last_timestamp_count: int = 0

    def add(self, result: DigitizedEvent) -> None:
        self.count += 1
        for name in (
            "n_selected",
            "n_detected",
            "n_accepted",
            "n_lost_deadtime",
            "n_active_spads",
            "n_active_timing_groups",
            "n_timing_records",
            "retained_fraction",
            "deadtime_loss_fraction",
            "mean_pde",
            "payload_bits_proxy",
        ):
            self.totals[name] = self.totals.get(name, 0.0) + float(
                getattr(result, name)
            )
        self.accepted_squared += float(result.n_accepted) ** 2
        if result.first_timestamp_ns is not None:
            self.totals["first_timestamp_ns"] = self.totals.get(
                "first_timestamp_ns", 0.0
            ) + float(result.first_timestamp_ns)
            self.first_timestamp_count += 1
        if result.last_timestamp_ns is not None:
            self.totals["last_timestamp_ns"] = self.totals.get(
                "last_timestamp_ns", 0.0
            ) + float(result.last_timestamp_ns)
            self.last_timestamp_count += 1

    def mean(self, name: str) -> float:
        if self.count == 0:
            return float("nan")
        return self.totals.get(name, 0.0) / self.count

    def accepted_std(self) -> float:
        if self.count == 0:
            return float("nan")
        variance = self.accepted_squared / self.count - self.mean("n_accepted") ** 2
        return math.sqrt(max(variance, 0.0))

    def aggregate_ratio(self, numerator: str, denominator: str) -> float:
        total_denominator = self.totals.get(denominator, 0.0)
        if total_denominator == 0:
            return 0.0
        return self.totals.get(numerator, 0.0) / total_denominator


class ComparisonOutputs:
    """Stream event rows while accumulating compact plotting products."""

    def __init__(
        self,
        output_dir: str | Path,
        devices: list[DeviceConfig],
        *,
        plot_style: str = "comparison",
    ):
        import ROOT

        ROOT.gROOT.SetBatch(True)
        if plot_style not in PLOT_STYLES:
            available = ", ".join(sorted(PLOT_STYLES))
            raise ValueError(f"Unknown plot style {plot_style!r}. Choose from {available}")
        self.ROOT = ROOT
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.devices = devices
        self.plot_style = plot_style
        self.accumulators = {device.name: Accumulator() for device in devices}

        self._event_handle = (self.output_dir / "event_summary.csv").open(
            "w", newline=""
        )
        self._event_writer = csv.DictWriter(self._event_handle, fieldnames=EVENT_FIELDS)
        self._event_writer.writeheader()

        self.root_path = self.output_dir / "comparison.root"
        self.root_file = ROOT.TFile(str(self.root_path), "RECREATE")
        self.histograms: dict[str, dict[str, object]] = {}
        for device in devices:
            self.histograms[device.name] = self._make_histograms(device)
        self.historical_outputs = (
            HistoricalSlideOutputs(self.root_file, self.output_dir, devices)
            if plot_style in {"slides", "both"}
            else None
        )

    def _make_histograms(self, device: DeviceConfig) -> dict[str, object]:
        ROOT = self.ROOT
        prefix = device.name
        histograms = {
            "selected": ROOT.TH1D(
                f"{prefix}_selected_per_event",
                f"{device.label};Selected photons/event;Events",
                200,
                0,
                20000,
            ),
            "accepted": ROOT.TH1D(
                f"{prefix}_accepted_per_event",
                f"{device.label};Accepted photons/event;Events",
                200,
                0,
                20000,
            ),
            "retained": ROOT.TH1D(
                f"{prefix}_retained_fraction",
                f"{device.label};Accepted / selected;Events",
                101,
                0,
                1.01,
            ),
            "deadtime": ROOT.TH1D(
                f"{prefix}_deadtime_loss_fraction",
                f"{device.label};Dead-time loss / detected;Events",
                101,
                0,
                1.01,
            ),
            "records": ROOT.TH1D(
                f"{prefix}_timing_records_per_event",
                f"{device.label};Timing records/event;Events",
                200,
                0,
                20000,
            ),
            "payload": ROOT.TH1D(
                f"{prefix}_payload_bits_proxy",
                f"{device.label};Payload proxy [bits/event];Events",
                200,
                0,
                10000000,
            ),
            "detected_time": ROOT.TH1D(
                f"{prefix}_detected_time",
                f"{device.label};Photon time [ns];Detected photons",
                700,
                5,
                40,
            ),
            "accepted_time": ROOT.TH1D(
                f"{prefix}_accepted_time",
                f"{device.label};Recorded time [ns];Accepted photons",
                700,
                5,
                40,
            ),
            "multiplicity": ROOT.TH1D(
                f"{prefix}_photons_per_hit_spad",
                f"{device.label};Detected photons/hit SPAD;Hit SPADs",
                30,
                0,
                30,
            ),
        }
        for histogram in histograms.values():
            histogram.SetDirectory(self.root_file)
        return histograms

    def add(self, result: DigitizedEvent) -> None:
        self._event_writer.writerow(result.as_row())
        self.accumulators[result.device_name].add(result)
        histograms = self.histograms[result.device_name]
        histograms["selected"].Fill(result.n_selected)
        histograms["accepted"].Fill(result.n_accepted)
        histograms["retained"].Fill(result.retained_fraction)
        histograms["deadtime"].Fill(result.deadtime_loss_fraction)
        histograms["records"].Fill(result.n_timing_records)
        histograms["payload"].Fill(result.payload_bits_proxy)
        for value in result.detected_time_ns:
            histograms["detected_time"].Fill(float(value))
        for value in result.accepted_time_ns:
            histograms["accepted_time"].Fill(float(value))
        for value in result.cell_multiplicities:
            histograms["multiplicity"].Fill(float(value))
        if self.historical_outputs is not None:
            self.historical_outputs.add(result)

    def _summary_rows(self) -> list[dict[str, object]]:
        rows = []
        for device in self.devices:
            accumulator = self.accumulators[device.name]
            first_mean = (
                accumulator.totals.get("first_timestamp_ns", 0.0)
                / accumulator.first_timestamp_count
                if accumulator.first_timestamp_count
                else ""
            )
            last_mean = (
                accumulator.totals.get("last_timestamp_ns", 0.0)
                / accumulator.last_timestamp_count
                if accumulator.last_timestamp_count
                else ""
            )
            rows.append(
                {
                    "device": device.name,
                    "label": device.label,
                    "category": device.category,
                    "maturity": device.maturity,
                    "spad_pitch_um": device.spad_pitch_um,
                    "timing_group": (
                        f"{device.timing_group_spads_x}x"
                        f"{device.timing_group_spads_y} SPADs"
                    ),
                    "dead_time_ns": ""
                    if device.dead_time_ns is None
                    else device.dead_time_ns,
                    "pde_model": device.pde_model,
                    "pde_basis": device.pde_basis,
                    "timestamp_policy": device.timestamp_policy,
                    "events": accumulator.count,
                    "mean_selected": accumulator.mean("n_selected"),
                    "mean_detected": accumulator.mean("n_detected"),
                    "mean_accepted": accumulator.mean("n_accepted"),
                    "std_accepted": accumulator.accepted_std(),
                    "mean_retained_fraction": accumulator.mean("retained_fraction"),
                    "mean_deadtime_loss_fraction": accumulator.mean(
                        "deadtime_loss_fraction"
                    ),
                    "aggregate_retained_fraction": accumulator.aggregate_ratio(
                        "n_accepted", "n_selected"
                    ),
                    "aggregate_deadtime_loss_fraction": accumulator.aggregate_ratio(
                        "n_lost_deadtime", "n_detected"
                    ),
                    "mean_active_spads": accumulator.mean("n_active_spads"),
                    "mean_active_timing_groups": accumulator.mean(
                        "n_active_timing_groups"
                    ),
                    "mean_timing_records": accumulator.mean("n_timing_records"),
                    "mean_first_timestamp_ns": first_mean,
                    "mean_last_timestamp_ns": last_mean,
                    "mean_payload_bits_proxy": accumulator.mean(
                        "payload_bits_proxy"
                    ),
                    "source_url": device.source_url,
                    "notes": device.notes,
                }
            )
        return rows

    def _write_summary(self, rows: list[dict[str, object]]) -> None:
        if not rows:
            return
        with (self.output_dir / "device_summary.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def _draw_bar(
        self,
        rows: list[dict[str, object]],
        field: str,
        title: str,
        y_title: str,
        filename: str,
        *,
        log_y: bool = False,
        maximum: float | None = None,
    ) -> None:
        ROOT = self.ROOT
        canvas = ROOT.TCanvas(f"canvas_{field}", title, 1400, 900)
        canvas.SetBottomMargin(0.42)
        canvas.SetLeftMargin(0.11)
        canvas.SetGridy()
        canvas.SetLogy(log_y)
        histogram = ROOT.TH1D(
            f"summary_{field}",
            f"{title};;{y_title}",
            len(rows),
            0,
            len(rows),
        )
        histogram.SetDirectory(self.root_file)
        for index, row in enumerate(rows, start=1):
            histogram.SetBinContent(index, float(row[field]))
            histogram.GetXaxis().SetBinLabel(index, str(row["label"]))
        histogram.SetFillColor(ROOT.kAzure - 9)
        histogram.SetLineColor(ROOT.kAzure + 2)
        histogram.SetLineWidth(2)
        histogram.SetBarWidth(0.78)
        histogram.SetBarOffset(0.11)
        histogram.SetStats(0)
        if maximum is not None:
            histogram.SetMaximum(maximum)
        if log_y:
            positive = [float(row[field]) for row in rows if float(row[field]) > 0]
            if positive:
                histogram.SetMinimum(min(positive) * 0.5)
        histogram.LabelsOption("v", "X")
        histogram.GetXaxis().SetLabelSize(0.035)
        histogram.GetXaxis().SetLabelOffset(0.012)
        histogram.GetXaxis().SetTitleOffset(2.8)
        histogram.Draw("BAR2")
        canvas.SaveAs(str(self.output_dir / filename))

    def _draw_overlay(
        self,
        histogram_key: str,
        title: str,
        x_title: str,
        filename: str,
        *,
        normalize: bool,
        log_y: bool,
    ) -> None:
        ROOT = self.ROOT
        canvas = ROOT.TCanvas(f"canvas_{histogram_key}", title, 1100, 750)
        canvas.SetLeftMargin(0.11)
        canvas.SetRightMargin(0.28)
        canvas.SetGridx()
        canvas.SetGridy()
        canvas.SetLogy(log_y)
        legend = ROOT.TLegend(0.74, 0.15, 0.99, 0.90)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        colors = [
            ROOT.kBlue + 1,
            ROOT.kRed + 1,
            ROOT.kGreen + 2,
            ROOT.kMagenta + 1,
            ROOT.kOrange + 7,
            ROOT.kCyan + 2,
            ROOT.kViolet + 2,
            ROOT.kSpring + 5,
            ROOT.kPink + 6,
            ROOT.kAzure + 7,
            ROOT.kGray + 2,
        ]
        drawn = False
        clones = []
        for index, device in enumerate(self.devices):
            original = self.histograms[device.name][histogram_key]
            clone = original.Clone(f"plot_{device.name}_{histogram_key}")
            clone.SetDirectory(0)
            if normalize and clone.Integral() > 0:
                clone.Scale(1.0 / clone.Integral())
            clone.SetLineColor(colors[index % len(colors)])
            clone.SetLineWidth(2)
            clone.SetStats(0)
            clone.SetTitle(f"{title};{x_title};{'Normalized entries' if normalize else 'Entries'}")
            if log_y:
                clone.SetMinimum(1e-8 if normalize else 0.5)
            clone.Draw("HIST" if not drawn else "HIST SAME")
            legend.AddEntry(clone, device.label, "l")
            clones.append(clone)
            drawn = True
        if drawn:
            legend.Draw()
            canvas.SaveAs(str(self.output_dir / filename))

    def close(self) -> list[dict[str, object]]:
        self._event_handle.close()
        rows = self._summary_rows()
        self._write_summary(rows)

        if self.plot_style in {"comparison", "both"}:
            self._draw_bar(
                rows,
                "aggregate_retained_fraction",
                "Photon retention after PDE and recovery",
                "Mean accepted / selected",
                "retained_fraction.png",
                maximum=1.05,
            )
            self._draw_bar(
                rows,
                "aggregate_deadtime_loss_fraction",
                "Dead-time loss among detected photons",
                "Mean lost / detected",
                "deadtime_loss_fraction.png",
                maximum=1.05,
            )
            self._draw_bar(
                rows,
                "mean_timing_records",
                "Timing information exported",
                "Mean records/event",
                "timing_records.png",
                log_y=True,
            )
            self._draw_bar(
                rows,
                "mean_payload_bits_proxy",
                "Transparent payload-size proxy",
                "Mean bits/event",
                "payload_bits_proxy.png",
                log_y=True,
            )
            self._draw_overlay(
                "accepted_time",
                "Recorded photon-time distributions",
                "Recorded time [ns]",
                "accepted_time_overlay.png",
                normalize=True,
                log_y=True,
            )
            self._draw_overlay(
                "multiplicity",
                "Incident multiplicity in hit SPADs",
                "Detected photons per hit SPAD",
                "photons_per_spad_overlay.png",
                normalize=True,
                log_y=True,
            )

        if self.historical_outputs is not None:
            self.historical_outputs.close(rows)

        self.root_file.Write()
        self.root_file.Close()
        return rows

    def __enter__(self) -> "ComparisonOutputs":
        return self

    def __exit__(self, exception_type, exception, traceback) -> None:
        if exception_type is None:
            self.close()
        else:
            self._event_handle.close()
            self.root_file.Close()


def print_summary(rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        return
    print()
    print(
        f"{'Device':26s} {'accepted':>11s} {'retained':>10s} "
        f"{'DT loss':>10s} {'records':>10s} {'bits':>12s}"
    )
    print("-" * 85)
    for row in rows:
        print(
            f"{str(row['label'])[:26]:26s} "
            f"{float(row['mean_accepted']):11.2f} "
            f"{float(row['aggregate_retained_fraction']):10.4f} "
            f"{float(row['aggregate_deadtime_loss_fraction']):10.4f} "
            f"{float(row['mean_timing_records']):10.2f} "
            f"{float(row['mean_payload_bits_proxy']):12.1f}"
        )
