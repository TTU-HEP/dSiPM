"""CMS-inspired ROOT helpers for FLASH/HG-DREAM figures.

The defaults are intentionally presentation friendly: the typography and
label hierarchy follow the public CMS plots used in Chris Madrid's recent
talks, while the wider canvases fit naturally into a 16:9 slide.  Individual
analysis scripts should use :func:`style_pad`, :func:`style_axes`, and
:func:`make_legend` instead of reimplementing margins and font sizes.
"""

from __future__ import annotations

import shutil
import subprocess


def set_flash_style(ROOT) -> dict[str, int]:
    """Apply a clean TDR-like style and return the project color palette."""

    ROOT.gROOT.SetBatch(True)
    style = ROOT.TStyle("FLASHStyle", "FLASH CMS-inspired style")
    style.SetCanvasBorderMode(0)
    style.SetCanvasColor(ROOT.kWhite)
    style.SetCanvasDefH(720)
    style.SetCanvasDefW(900)
    style.SetPadBorderMode(0)
    style.SetPadColor(ROOT.kWhite)
    style.SetPadGridX(False)
    style.SetPadGridY(False)
    style.SetFrameBorderMode(0)
    style.SetFrameFillColor(ROOT.kWhite)
    style.SetFrameLineWidth(2)
    style.SetOptStat(0)
    style.SetOptTitle(0)
    style.SetTitleFont(42, "XYZ")
    style.SetTitleSize(0.050, "XYZ")
    style.SetTitleXOffset(1.08)
    style.SetTitleYOffset(1.35)
    style.SetLabelFont(42, "XYZ")
    style.SetLabelSize(0.041, "XYZ")
    style.SetLabelOffset(0.008, "XYZ")
    style.SetAxisColor(ROOT.kBlack, "XYZ")
    style.SetStripDecimals(True)
    style.SetTickLength(0.028, "XYZ")
    style.SetNdivisions(510, "XYZ")
    style.SetPadTickX(1)
    style.SetPadTickY(1)
    style.SetPadTopMargin(0.105)
    style.SetPadBottomMargin(0.150)
    style.SetPadLeftMargin(0.155)
    style.SetPadRightMargin(0.055)
    style.SetLegendBorderSize(0)
    style.SetLegendFillColor(ROOT.kWhite)
    style.SetLegendFillStyle(0)
    style.SetLegendFont(42)
    style.SetHistLineWidth(3)
    style.SetLineWidth(2)
    style.SetMarkerStyle(20)
    style.SetMarkerSize(1.0)
    style.SetEndErrorSize(4)
    style.SetErrorX(0.0)
    ROOT.gROOT.SetStyle("FLASHStyle")
    ROOT.gROOT.ForceStyle()
    return {
        "black": ROOT.kBlack,
        "red": ROOT.TColor.GetColor("#BD1F01"),
        "blue": ROOT.TColor.GetColor("#3F90DA"),
        "green": ROOT.TColor.GetColor("#00A6A6"),
        "orange": ROOT.TColor.GetColor("#FFA90E"),
        "purple": ROOT.TColor.GetColor("#832DB6"),
        "gray": ROOT.TColor.GetColor("#5F6368"),
        "light_gray": ROOT.TColor.GetColor("#D9D9D9"),
        "yellow": ROOT.TColor.GetColor("#E5AE38"),
    }


def make_canvas(
    ROOT,
    name: str,
    *,
    width: int = 900,
    height: int = 720,
    right: float | None = None,
):
    canvas = ROOT.TCanvas(name, name, width, height)
    canvas.SetTicks(1, 1)
    if right is not None:
        canvas.SetRightMargin(right)
    return canvas


def style_pad(
    pad,
    *,
    left: float = 0.155,
    right: float = 0.045,
    top: float = 0.125,
    bottom: float = 0.155,
    grid_y: bool = False,
) -> None:
    """Apply consistent margins and ticks to a ROOT pad."""

    pad.SetLeftMargin(left)
    pad.SetRightMargin(right)
    pad.SetTopMargin(top)
    pad.SetBottomMargin(bottom)
    pad.SetTicks(1, 1)
    pad.SetGridx(False)
    pad.SetGridy(grid_y)


def style_axes(histogram, *, title_size: float = 0.050, label_size: float = 0.041) -> None:
    """Style every axis owned by a histogram or graph."""

    for axis_name in ("X", "Y", "Z"):
        getter = getattr(histogram, f"Get{axis_name}axis", None)
        if getter is None:
            continue
        axis = getter()
        if not axis:
            continue
        axis.SetTitleFont(42)
        axis.SetLabelFont(42)
        axis.SetTitleSize(title_size)
        axis.SetLabelSize(label_size)
        axis.SetLabelOffset(0.008)
    histogram.GetXaxis().SetTitleOffset(1.10)
    histogram.GetYaxis().SetTitleOffset(1.38)
    get_z_axis = getattr(histogram, "GetZaxis", None)
    if get_z_axis is not None and get_z_axis():
        get_z_axis().SetTitleOffset(1.30)


def make_legend(
    ROOT,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    *,
    text_size: float = 0.038,
    columns: int = 1,
):
    """Return a transparent, borderless CMS-style legend."""

    legend = ROOT.TLegend(x1, y1, x2, y2)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(text_size)
    legend.SetNColumns(columns)
    return legend


def draw_flash_label(
    ROOT,
    pad,
    *,
    extra: str = "Simulation",
    right: str = "",
    x: float | None = None,
    y: float = 0.955,
) -> list[object]:
    """Draw a CMS-style experiment/qualifier line in normalized coordinates."""

    pad.cd()
    left = pad.GetLeftMargin() if x is None else x
    keep: list[object] = []
    flash = ROOT.TLatex()
    flash.SetNDC(True)
    flash.SetTextFont(62)
    flash.SetTextSize(0.044)
    flash.DrawLatex(left, y, "FLASH")
    keep.append(flash)
    if extra:
        qualifier = ROOT.TLatex()
        qualifier.SetNDC(True)
        qualifier.SetTextFont(52)
        qualifier.SetTextSize(0.033)
        qualifier.DrawLatex(left + 0.150, y, extra)
        keep.append(qualifier)
    if right:
        label = ROOT.TLatex()
        label.SetNDC(True)
        label.SetTextFont(42)
        label.SetTextAlign(31)
        label.SetTextSize(0.031)
        label.DrawLatex(1.0 - pad.GetRightMargin(), y, right)
        keep.append(label)
    return keep


def draw_panel_title(
    ROOT,
    pad,
    title: str,
    *,
    subtitle: str = "",
    x: float | None = None,
    y: float = 0.955,
) -> list[object]:
    """Draw a bold panel heading with an optional compact subtitle."""

    pad.cd()
    left = pad.GetLeftMargin() if x is None else x
    keep: list[object] = []
    heading = ROOT.TLatex()
    heading.SetNDC(True)
    heading.SetTextFont(62)
    heading.SetTextSize(0.043)
    heading.DrawLatex(left, y, title)
    keep.append(heading)
    if subtitle:
        detail = ROOT.TLatex()
        detail.SetNDC(True)
        detail.SetTextFont(42)
        detail.SetTextColor(ROOT.TColor.GetColor("#5F6368"))
        detail.SetTextSize(0.030)
        detail.DrawLatex(left, y - 0.055, subtitle)
        keep.append(detail)
    return keep


def save_canvas(canvas, stem) -> list:
    pdf_output = stem.with_suffix(".pdf")
    png_output = stem.with_suffix(".png")
    outputs = [pdf_output, png_output]
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs(str(png_output))

    # ROOT's direct PDF backend uses the global paper size rather than the
    # canvas aspect ratio.  Wide, slide-ready canvases therefore land at the
    # top of portrait pages.  Export through SVG to preserve vector content
    # and obtain a PDF page whose geometry matches the canvas.
    converter = shutil.which("rsvg-convert")
    if converter:
        svg_output = stem.with_suffix(".tmp.svg")
        temporary_pdf = stem.with_suffix(".tmp.pdf")
        try:
            canvas.SaveAs(str(svg_output))
            subprocess.run(
                [
                    converter,
                    "-f",
                    "pdf1.5",
                    "-o",
                    str(temporary_pdf),
                    str(svg_output),
                ],
                check=True,
            )
            temporary_pdf.replace(pdf_output)
        finally:
            svg_output.unlink(missing_ok=True)
            temporary_pdf.unlink(missing_ok=True)
    else:
        canvas.SaveAs(str(pdf_output))
    return outputs
