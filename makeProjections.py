import ROOT

def doProjections(h3d):
    base_name = h3d.GetName()

    hist_xy = h3d.Project3D("xy")
    hist_xy.SetName(base_name + "_projXY")

    hist_xz = h3d.Project3D("xz")
    hist_xz.SetName(base_name + "_projXZ")

    hist_yz = h3d.Project3D("yz")
    hist_yz.SetName(base_name + "_projYZ")    
    return hist_xy, hist_xz, hist_yz

def makeEff(num, den, c):
    eff = num.Clone(num.GetName()+"eff")
    eff.Divide(num, den, 1.0, 1.0, "B")  # Binomial errors
    eff.Write()
    eff.SetStats(0)
    eff.SetMaximum(1.0)
    eff.Draw("COLZ")
    c.SaveAs("output/" + eff.GetName() + ".png")  
    return eff

def main(input_filename):
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile.Open(input_filename)
    histogram_suffix = "nPhotons_xyt_all_"

    # Get list of keys and filter histograms
    hist_names = [key.GetName() for key in f.GetListOfKeys()]
    matching_hists = [name.replace(histogram_suffix,"") for name in hist_names if name.startswith(histogram_suffix)]

    # Save and plot
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)

    out_file = ROOT.TFile("projections.root", "RECREATE")

    for i, hist_name in enumerate(matching_hists):
        histNum = f.Get("nPhotons_xyt_oneHit_" + hist_name)
        histDen = f.Get(   "nPhotons_xyt_all_" + hist_name)

        #Make 2D histograms
        hNum_xy, hNum_xz, hNum_yz = doProjections(histNum)
        hDen_xy, hDen_xz, hDen_yz = doProjections(histDen)

        #Make eff histograms
        eff_xy = makeEff(hNum_xy, hDen_xy, c1)
        eff_xz = makeEff(hNum_xz, hDen_xz, c1)
        eff_yz = makeEff(hNum_yz, hDen_yz, c1)

        #Save things
        for hist in [hNum_xy, hNum_xz, hNum_yz, hDen_xy, hDen_xz, hDen_yz]:
            hist.Write()
            hist.SetStats(0)
            hist.Draw("COLZ")
            c1.SaveAs("output/" + hist.GetName() + ".png")

    out_file.Close()

if __name__ == "__main__":
    main("outfile.root")
