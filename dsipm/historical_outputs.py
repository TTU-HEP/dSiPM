"""Slide-matched ROOT products for Chris Madrid's 2025 occupancy study."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .config import DeviceConfig
from .digitizer import DigitizedEvent


HISTORICAL_BINNING = {
    "sensor_physical_width_mm": 3.0,
    "sensor_display_xy_mm": [-0.5, 0.5],
    "xy_bins": "round(3000 um / configured channel pitch)",
    "coordinate_scale": 3.0,
    "radius_mm": [0.0, 0.5],
    "radius_bins": 60,
    "time_ns": [5.0, 40.0],
    "time_bins": 700,
    "time_bin_ps": 50.0,
    "production_z_mm": [0.0, 2000.0],
    "production_z_bins": 500,
    "production_z_bin_mm": 4.0,
    "photons_per_channel": [0.0, 30.0],
    "photons_per_channel_bins": 30,
}


def _pitch_tag(device: DeviceConfig) -> str:
    pitch = device.spad_pitch_um
    value = str(int(pitch)) if float(pitch).is_integer() else f"{pitch:g}"
    value = value.replace(".", "p")
    return f"{value}x{value}"


def _pitch_label(device: DeviceConfig) -> str:
    pitch = device.spad_pitch_um
    value = str(int(pitch)) if float(pitch).is_integer() else f"{pitch:g}"
    return f"{value}#times{value} #mum^{{2}}"


class HistoricalSlideOutputs:
    """Accumulate the same observables and binning used in the 2025 slides.

    The original scripts stored a third axis that the slide plots immediately
    projected away.  A single aggregate bin preserves those ROOT names and the
    plotted projections while avoiding several gigabytes of redundant storage.
    """

    def __init__(self, root_file, output_dir: str | Path, devices: list[DeviceConfig]):
        import ROOT

        self.ROOT = ROOT
        self.root_file = root_file
        self.output_dir = Path(output_dir)
        self.devices = [
            device
            for device in devices
            if device.category == "historical idealization"
        ]
        self.map_devices = [
            device for device in self.devices if device.spad_pitch_um >= 10.0
        ]
        self.device_by_name = {device.name: device for device in self.map_devices}
        self.histograms = {
            device.name: self._make_histograms(device) for device in self.map_devices
        }
        with (self.output_dir / "historical_binning.json").open("w") as handle:
            json.dump(HISTORICAL_BINNING, handle, indent=2)
            handle.write("\n")

    def _make_histograms(self, device: DeviceConfig) -> dict[str, object]:
        ROOT = self.ROOT
        tag = _pitch_tag(device)
        pitch_mm = device.spad_pitch_um / 1000.0
        xy_bins = max(1, int(round(3.0 / pitch_mm)))

        definitions = {
            "xyt_all": ROOT.TH3D(
                f"nPhotons_xyt_all_{tag}",
                "Incident photons; y [mm]; x [mm]; t [ns]",
                xy_bins,
                -0.5,
                0.5,
                xy_bins,
                -0.5,
                0.5,
                1,
                5.0,
                40.0,
            ),
            "rte_all": ROOT.TH3D(
                f"nPhotons_rte_all_{tag}",
                "Incident photons; r [mm]; t [ns]; aggregate event bin",
                60,
                0.0,
                0.5,
                700,
                5.0,
                40.0,
                1,
                0.0,
                1.0,
            ),
            "rze_all": ROOT.TH3D(
                f"nPhotons_rze_all_{tag}",
                "Incident photons; r [mm]; z [mm]; aggregate event bin",
                60,
                0.0,
                0.5,
                500,
                0.0,
                2000.0,
                1,
                0.0,
                1.0,
            ),
            "tze_all": ROOT.TH3D(
                f"nPhotons_tze_all_{tag}",
                "Incident photons; t [ns]; z [mm]; aggregate event bin",
                700,
                5.0,
                40.0,
                500,
                0.0,
                2000.0,
                1,
                0.0,
                1.0,
            ),
            "multiplicity": ROOT.TH1D(
                f"nPhotonsPerChannel_{tag}",
                "Photons per channel;NPhotons per channel;Channels",
                30,
                0.0,
                30.0,
            ),
        }
        for stem, title in (
            ("xyt", "Accepted photons; y [mm]; x [mm]; t [ns]"),
            ("rte", "Accepted photons; r [mm]; t [ns]; aggregate event bin"),
            ("rze", "Accepted photons; r [mm]; z [mm]; aggregate event bin"),
            ("tze", "Accepted photons; t [ns]; z [mm]; aggregate event bin"),
        ):
            source = definitions[f"{stem}_all"]
            accepted = source.Clone(f"nPhotons_{stem}_oneHit_{tag}")
            accepted.Reset()
            accepted.SetTitle(title)
            definitions[f"{stem}_oneHit"] = accepted

        for histogram in definitions.values():
            histogram.SetDirectory(self.root_file)
        for suffix in ("all", "oneHit"):
            for stem in ("xyt", "rte", "rze", "tze"):
                source = definitions[f"{stem}_{suffix}"]
                accumulator = ROOT.TH2D(
                    f"aggregate_{source.GetName()}",
                    source.GetTitle(),
                    source.GetXaxis().GetNbins(),
                    source.GetXaxis().GetXmin(),
                    source.GetXaxis().GetXmax(),
                    source.GetYaxis().GetNbins(),
                    source.GetYaxis().GetXmin(),
                    source.GetYaxis().GetXmax(),
                )
                accumulator.SetDirectory(0)
                definitions[f"{stem}_{suffix}_2d"] = accumulator
        return definitions

    @staticmethod
    def _fill_2d(histogram, x: np.ndarray, y: np.ndarray) -> None:
        if x.size == 0:
            return
        x_values = np.ascontiguousarray(x, dtype=np.float64)
        y_values = np.ascontiguousarray(y, dtype=np.float64)
        weights = np.ones(x_values.size, dtype=np.float64)
        histogram.FillN(int(x_values.size), x_values, y_values, weights)

    def _fill_maps(
        self,
        histograms: dict[str, object],
        suffix: str,
        x_mm: np.ndarray,
        y_mm: np.ndarray,
        time_ns: np.ndarray,
        production_z_mm: np.ndarray,
    ) -> None:
        radius_mm = np.hypot(x_mm, y_mm)
        self._fill_2d(histograms[f"xyt_{suffix}_2d"], y_mm, x_mm)
        self._fill_2d(
            histograms[f"rte_{suffix}_2d"], radius_mm, time_ns
        )
        self._fill_2d(
            histograms[f"rze_{suffix}_2d"],
            radius_mm,
            production_z_mm,
        )
        self._fill_2d(
            histograms[f"tze_{suffix}_2d"],
            time_ns,
            production_z_mm,
        )

    def _fill_multiplicity(
        self,
        histogram,
        device: DeviceConfig,
        x_mm: np.ndarray,
        y_mm: np.ndarray,
    ) -> None:
        pitch_mm = device.spad_pitch_um / 1000.0
        xy_bins = max(1, int(round(3.0 / pitch_mm)))
        in_sensor = (
            (x_mm >= -0.5)
            & (x_mm < 0.5)
            & (y_mm >= -0.5)
            & (y_mm < 0.5)
        )
        sensor_x = x_mm[in_sensor]
        sensor_y = y_mm[in_sensor]
        if sensor_x.size:
            ix = np.floor((sensor_x + 0.5) * xy_bins).astype(np.int64)
            iy = np.floor((sensor_y + 0.5) * xy_bins).astype(np.int64)
            cell_ids = iy * xy_bins + ix
            _, counts = np.unique(cell_ids, return_counts=True)
        else:
            counts = np.empty(0, dtype=np.int64)

        zero_channels = xy_bins * xy_bins - int(counts.size)
        if zero_channels:
            histogram.Fill(0.0, float(zero_channels))
        for count in counts:
            histogram.Fill(float(count))

    def add(self, result: DigitizedEvent) -> None:
        if result.device_name not in self.histograms:
            return
        device = self.device_by_name[result.device_name]
        histograms = self.histograms[result.device_name]
        self._fill_maps(
            histograms,
            "all",
            result.detected_x_mm,
            result.detected_y_mm,
            result.detected_time_ns,
            result.detected_production_z_mm,
        )
        self._fill_maps(
            histograms,
            "oneHit",
            result.accepted_x_mm,
            result.accepted_y_mm,
            result.accepted_time_true_ns,
            result.accepted_production_z_mm,
        )
        self._fill_multiplicity(
            histograms["multiplicity"],
            device,
            result.detected_x_mm,
            result.detected_y_mm,
        )

    def _project(self, histogram, option: str, name: str):
        projection = histogram.Project3D(option)
        projection.SetName(name)
        projection.SetDirectory(0)
        projection.SetStats(0)
        return projection

    def _copy_accumulators_to_root_histograms(self) -> None:
        ROOT = self.ROOT
        if not hasattr(ROOT, "dSiPMCopyTH2ToTH3"):
            ROOT.gInterpreter.Declare(
                r"""
                void dSiPMCopyTH2ToTH3(const TH2D* source, TH3D* target) {
                    for (int ix = 0; ix <= source->GetNbinsX() + 1; ++ix) {
                        for (int iy = 0; iy <= source->GetNbinsY() + 1; ++iy) {
                            const double value = source->GetBinContent(ix, iy);
                            if (value != 0.0) {
                                target->SetBinContent(ix, iy, 1, value);
                            }
                        }
                    }
                    target->SetEntries(source->GetEntries());
                }
                """
            )
        for histograms in self.histograms.values():
            for suffix in ("all", "oneHit"):
                for stem in ("xyt", "rte", "rze", "tze"):
                    ROOT.dSiPMCopyTH2ToTH3(
                        histograms[f"{stem}_{suffix}_2d"],
                        histograms[f"{stem}_{suffix}"],
                    )

    def _map_triplet(self, device: DeviceConfig, *, efficiency: bool) -> None:
        ROOT = self.ROOT
        histograms = self.histograms[device.name]
        pitch = int(device.spad_pitch_um)
        suffix = "deadtime_efficiency" if efficiency else "photon_distributions"
        title = (
            f"Full-event one-hit efficiency: {_pitch_label(device)}"
            if efficiency
            else f"Incident photon distributions: {_pitch_label(device)}"
        )
        canvas = ROOT.TCanvas(f"canvas_{suffix}_{pitch}", title, 1800, 600)
        canvas.Divide(3, 1)
        specifications = (
            ("xyt", "yx", "Sensor plane; y [mm]; x [mm]"),
            ("rte", "xy", "Arrival time vs radius; t [ns]; r [mm]"),
            ("tze", "xy", "Arrival time vs production depth; z [mm]; t [ns]"),
        )
        retained = []
        for pad_index, (stem, option, plot_title) in enumerate(specifications, 1):
            denominator = self._project(
                histograms[f"{stem}_all"],
                option,
                f"plot_{stem}_all_{pitch}_{suffix}",
            )
            plot = denominator
            retained.append(denominator)
            if efficiency:
                numerator = self._project(
                    histograms[f"{stem}_oneHit"],
                    option,
                    f"plot_{stem}_oneHit_{pitch}_{suffix}",
                )
                plot = numerator.Clone(f"plot_{stem}_efficiency_{pitch}")
                plot.SetDirectory(0)
                plot.Divide(numerator, denominator, 1.0, 1.0, "B")
                plot.SetMinimum(0.0)
                plot.SetMaximum(1.0)
                retained.extend([numerator, plot])
            plot.SetTitle(plot_title)
            canvas.cd(pad_index)
            ROOT.gPad.SetLeftMargin(0.12)
            ROOT.gPad.SetRightMargin(0.17)
            ROOT.gPad.SetBottomMargin(0.12)
            ROOT.gPad.SetTicks(1, 1)
            plot.Draw("COLZ")
        canvas.SaveAs(str(self.output_dir / f"{suffix}_{pitch}um.png"))

    def _multiplicity_overlay(self) -> None:
        ROOT = self.ROOT
        canvas = ROOT.TCanvas(
            "canvas_historical_photons_per_channel",
            "Number of photons per channel",
            1100,
            750,
        )
        canvas.SetLogy()
        canvas.SetTicks(1, 1)
        canvas.SetLeftMargin(0.12)
        canvas.SetRightMargin(0.25)
        legend = ROOT.TLegend(0.76, 0.50, 0.99, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, ROOT.kOrange + 7]
        devices = [device for device in self.devices if device.spad_pitch_um >= 10.0]
        devices.sort(key=lambda device: device.spad_pitch_um)
        retained = []
        for index, device in enumerate(devices):
            source = self.histograms[device.name]["multiplicity"]
            plot = source.Clone(f"plot_multiplicity_{_pitch_tag(device)}")
            plot.SetDirectory(0)
            integral = plot.Integral(0, plot.GetNbinsX() + 1)
            if integral:
                plot.Scale(1.0 / integral)
            plot.SetLineColor(colors[index % len(colors)])
            plot.SetLineWidth(2)
            plot.SetStats(0)
            plot.SetTitle(
                "Number of photons per channel;NPhotons per channel;Normalized"
            )
            plot.GetXaxis().SetRangeUser(0.0, 15.0)
            plot.SetMinimum(1.0e-8)
            plot.SetMaximum(1.5)
            plot.Draw("HIST" if index == 0 else "HIST SAME")
            legend.AddEntry(plot, _pitch_label(device), "l")
            retained.append(plot)
        if retained:
            legend.Draw()
            canvas.SaveAs(str(self.output_dir / "photons_per_channel_overlay.png"))

    def _efficiency_vs_pitch(self, summary_rows: list[dict[str, object]]) -> None:
        ROOT = self.ROOT
        historical_names = {device.name for device in self.devices}
        points = sorted(
            (
                float(row["spad_pitch_um"]),
                1.0 - float(row["aggregate_deadtime_loss_fraction"]),
            )
            for row in summary_rows
            if row["device"] in historical_names
        )
        if not points:
            return
        x = np.asarray([point[0] for point in points], dtype=np.float64)
        y = np.asarray([point[1] for point in points], dtype=np.float64)
        graph = ROOT.TGraph(int(x.size), x, y)
        graph.SetName("historical_efficiency_vs_pitch")
        graph.SetTitle(
            "Full-event one-hit efficiency;Channel pitch [#mum];Accepted / incident"
        )
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1.3)
        graph.SetLineWidth(2)
        graph.SetMinimum(0.0)
        graph.SetMaximum(1.05)
        canvas = ROOT.TCanvas(
            "canvas_historical_efficiency_vs_pitch",
            "Efficiency versus channel pitch",
            900,
            650,
        )
        canvas.SetGridy()
        canvas.SetTicks(1, 1)
        canvas.SetLeftMargin(0.12)
        canvas.SetBottomMargin(0.12)
        graph.Draw("ALP")
        canvas.SaveAs(str(self.output_dir / "deadtime_efficiency_vs_pitch.png"))
        self.root_file.cd()
        graph.Write()

    def close(self, summary_rows: list[dict[str, object]]) -> None:
        if not self.devices:
            return
        self._copy_accumulators_to_root_histograms()
        self.ROOT.gStyle.SetOptStat(0)
        self.ROOT.gStyle.SetPalette(self.ROOT.kViridis)
        self._multiplicity_overlay()
        self._efficiency_vs_pitch(summary_rows)
        for target_pitch in (25.0, 100.0):
            matches = [
                device
                for device in self.devices
                if math.isclose(device.spad_pitch_um, target_pitch)
            ]
            if matches:
                self._map_triplet(matches[0], efficiency=False)
                self._map_triplet(matches[0], efficiency=True)
