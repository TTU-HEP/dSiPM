import ROOT
import re

def is_th2(hist):
    return isinstance(hist, ROOT.TH2)

def find_histogram_pairs(file):
    """Find (numerator, denominator) histogram pairs based on naming pattern."""
    hist_map = {key.GetName(): file.Get(key.GetName()) for key in file.GetListOfKeys() if is_th2(file.Get(key.GetName()))}

    pairs = []
    for name in hist_map:
        match = re.match(r"(.+)_oneHit_(.+)", name)
        if match:
            base, res = match.groups()
            denom_name = f"{base}_{res}"
            if denom_name in hist_map:
                num = hist_map[name].Clone()
                den = hist_map[denom_name].Clone()
                pairs.append((name, denom_name, num, den))
    return pairs

def compute_efficiencies(pairs):
    results = []
    for num_name, den_name, num_hist, den_hist in pairs:
        eff_name = num_name.replace("oneHit", "eff")
        eff_hist = num_hist.Clone(eff_name)
        eff_hist.Divide(num_hist, den_hist, 1.0, 1.0, "B")  # Binomial errors
        results.append((eff_name, eff_hist))
    return results

def main(root_file_path):
    ROOT.gROOT.SetBatch(True)
    file = ROOT.TFile.Open(root_file_path)

    # Find histogram pairs and compute efficiency
    pairs = find_histogram_pairs(file)
    efficiency_hists = compute_efficiencies(pairs)

    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)

    # Save to file and optionally draw
    out_file = ROOT.TFile("efficiencies.root", "RECREATE")
    for name, hist in efficiency_hists:
        hist.Write()
        hist.Draw("COLZ")
        c1.SaveAs(f"output/{name}.png")
    out_file.Close()
    print(f"Saved {len(efficiency_hists)} efficiency histograms to efficiencies.root")

if __name__ == "__main__":
    main("outfile.root")
