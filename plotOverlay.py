import ROOT
import re

# File and histogram base name
file_path = "outfile.root"
histogram_suffix = "nPhotonsPerChannel"

# Open ROOT file
f = ROOT.TFile.Open(file_path)
if not f or f.IsZombie():
    print("Error opening file!")
    exit()

# Get list of keys and filter histograms
hist_names = [key.GetName() for key in f.GetListOfKeys()]
matching_hists = [name for name in hist_names if name.startswith(histogram_suffix)]

def extract_sort_key(name):
    match = re.match(histogram_suffix + r"(\d+)x(\d+)_", name)
    if match:
        return int(match.group(1))  # Change to group(2) if sorting by second y
    return float('inf')  # Push non-matching names to the end

# Sort histograms by the extracted number
matching_hists.sort(key=extract_sort_key)

# Prepare canvas and legend
canvas = ROOT.TCanvas("canvas", "Overlayed Histograms", 800, 600)
canvas.SetLogy()
canvas.SetTickx()
canvas.SetTicky()
legend = ROOT.TLegend(0.65, 0.5, 0.88, 0.88)
legend.SetBorderSize(0)
legend.SetFillStyle(0)

# Colors for plotting
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange+7, ROOT.kCyan+2]

# Draw histograms
for i, hist_name in enumerate(matching_hists):
    hist = f.Get(hist_name)
    if not hist:
        continue

    if hist.Integral() > 0:
        hist.Scale(1.0/hist.Integral())

    hist.SetLineColor(colors[i % len(colors)])
    hist.SetLineWidth(2)
    hist.SetStats(0)

    # Set axis titles only once
    if i == 0:
        hist.GetXaxis().SetTitle("NPhotons per channel")
        hist.GetYaxis().SetTitle("Normalized")
        hist.SetMinimum(1e-6)
        hist.GetXaxis().SetRange(1,15)
        hist.Draw("HIST")
    else:
        hist.Draw("HIST SAME")

    # Add legend entry with prefix only
    prefix = hist_name.replace(histogram_suffix + "_", "") + (" #mum^{2}")
    legend.AddEntry(hist, prefix, "l")


f.Get(matching_hists[0]).Draw("AXIS SAME")
legend.Draw()
#canvas.SetGrid()
canvas.Update()

# Save output if needed
canvas.SaveAs("output/overlayed_histograms.png")
