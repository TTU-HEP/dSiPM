import ROOT

def doProjections(h3d, histType):
    base_name = h3d.GetName()

    hist_xy = h3d.Project3D("xy")
    hist_xy.SetName(base_name + "_projXY")

    if histType == "rte":
        hist_xy.GetXaxis().SetRangeUser(5, 20)
        return hist_xy, None, None

    hist_xz = h3d.Project3D("xz")
    hist_xz.SetName(base_name + "_projXZ")

    hist_yz = h3d.Project3D("yz")
    hist_yz.SetName(base_name + "_projYZ")    
    return hist_xy, hist_xz, hist_yz

def compute_2d_efficiency_manual(num_hist, den_hist, name="eff_hist"):
    if not (num_hist.GetNbinsX() == den_hist.GetNbinsX() and
            num_hist.GetNbinsY() == den_hist.GetNbinsY()):
        raise ValueError("Histograms must have the same binning.")

    eff_hist = num_hist.Clone(name)
    eff_hist.Reset()

    for ix in range(1, num_hist.GetNbinsX() + 1):
        for iy in range(1, num_hist.GetNbinsY() + 1):
            num = num_hist.GetBinContent(ix, iy)
            den = den_hist.GetBinContent(ix, iy)

            if den > 0:
                eff = num / den
                if num == 0: 
                    eff = 0.001
                eff_hist.SetBinContent(ix, iy, eff)
                # binomial error
                err = (eff * (1 - eff) / den)**0.5 if num > 0 else 0.0
                eff_hist.SetBinError(ix, iy, err)
            else:
                eff_hist.SetBinContent(ix, iy, 0.0)
                eff_hist.SetBinError(ix, iy, 0.0)

    return eff_hist

def makeEff(num, den, c):
    if num == None:
        return None
    #eff = num.Clone(num.GetName()+"eff")
    #eff.Divide(num, den, 1.0, 1.0, "B")  # Binomial errors

    eff = compute_2d_efficiency_manual(num, den, num.GetName()+"eff")

    eff.Write()
    eff.SetStats(0)
    eff.SetMaximum(1.0)
    eff.Draw("COLZ")
    c.SaveAs("output/" + eff.GetName() + ".png")  
    return eff

def makePlots(input_filename, histType):
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile.Open(input_filename)
    histogram_suffix = "nPhotons_{}_all_".format(histType)

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
        histNum = f.Get("nPhotons_{}_oneHit_{}".format(histType, hist_name))
        histDen = f.Get(   "nPhotons_{}_all_{}".format(histType, hist_name))

        #Make 2D histograms
        hNum_xy, hNum_xz, hNum_yz = doProjections(histNum, histType)
        hDen_xy, hDen_xz, hDen_yz = doProjections(histDen, histType)

        #Make eff histograms
        eff_xy = makeEff(hNum_xy, hDen_xy, c1)
        eff_xz = makeEff(hNum_xz, hDen_xz, c1)
        eff_yz = makeEff(hNum_yz, hDen_yz, c1)

        #Save things
        for hist in [hNum_xy, hNum_xz, hNum_yz, hDen_xy, hDen_xz, hDen_yz]:
            if hist == None: continue
            hist.Write()
            hist.SetStats(0)
            hist.Draw("COLZ")
            c1.SaveAs("output/" + hist.GetName() + ".png")

    out_file.Close()

def main(input_filename):
    #makePlots(input_filename, "rte")
    makePlots(input_filename, "xyt")

if __name__ == "__main__":
    main("outfile.root")
