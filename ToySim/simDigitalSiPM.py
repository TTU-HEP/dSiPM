import random
import math
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)

def makeBranch(tree, name, t):
    b = None
    if t == "vector":
        b = ROOT.vector('double')()
        tree.Branch(name, b)
    elif t == "D":
        b = np.zeros(1, dtype=float)
        tree.Branch(name, b, name+"/"+t)
    else:
        print("Branch type not supported yet. Add it here")
    return b
    
def saveHisto(c, hDic, name, drawOpt=""):
    hDic[name].Draw(drawOpt)
    c.SaveAs(name+".png")    

def rand(m=0,M=100):
    a = random.uniform(m,M)
    return a
    
def randgauss(m,w):
    return random.gauss(m,w)
    
def noise1D(histo):
    for b in range(histo.GetXaxis().GetNbins()+1):
        histo.SetBinContent(b, randgauss(0, 0.5))
    return histo

def hitOneArray(h, t):
    binNum = h.GetMaximumBin()
    b = h.FindBin(t[2], t[3])
    return binNum==b

def main():
    # Some useful hardcoded stuff
    nEvents = 1
    outfile = "dSIPMs.root"
    fiberWidth = 0.2

    # Define TTrees and TBranches
    root_file = ROOT.TFile(outfile, "RECREATE")
    tree = ROOT.TTree("tree","tree")
    amp_   = makeBranch(tree, "amp", "vector")
    time_   = makeBranch(tree, "time", "vector")
    
    # Define histos
    histos = {}

    nChannels = [
        ("500x500", 250), 
        ("100x100", 100), 
        ("50x50",    50), 
        ("25x25",    25), 
        ("10x10",    10), 
        ("5x5",       5), 
        ("1x1",       1),
    ]

    for c in nChannels: histos["nPhotons_xy_{}".format(c[0])] = ROOT.TH2D("nPhotons_xy_{}".format(c[0]),"nPhotons_xy; x [mm]; y [mm]; nPhotons", c[1],-1.5,1.5, c[1],-1.5,1.5)
    for c in nChannels: histos["nPhotons_time_{}".format(c[0])] = ROOT.TH1D("nPhotons_time_{}".format(c[0]),"nPhotons_time; time [ns]; nPhotons", 500,0.0,25.0)
    for c in nChannels: histos["signal_time_{}".format(c[0])] = ROOT.TH1D("signal_time_{}".format(c[0]),"signal_time; time [ns]; Amplitude [mV]", 500,0.0,25.0)

    # Loop over events
    for i in range(nEvents):
        #############################################
        # Generate photon distribution
        #############################################
        truthPhoton = [] #nphotons, time, x, y
        nSignal = 100
        for i in range(nSignal):
            truthPhoton.append((round(rand(1,10)), rand(0,25), randgauss(0.0, fiberWidth), randgauss(0.0,fiberWidth)))
            truthPhoton.append((round(rand(1,10)), rand(0,25), randgauss(1.0, fiberWidth), randgauss(0.0,fiberWidth)))
            truthPhoton.append((round(rand(1,10)), rand(0,25), randgauss(0.0, fiberWidth), randgauss(1.0,fiberWidth)))
            truthPhoton.append((round(rand(1,10)), rand(0,25), randgauss(1.0, fiberWidth), randgauss(1.0,fiberWidth)))
            truthPhoton.append((round(rand(1,10)), rand(0,25), randgauss(-1.0, fiberWidth), randgauss(0.0,fiberWidth)))
            truthPhoton.append((round(rand(1,10)), rand(0,25), randgauss(0.0, fiberWidth), randgauss(-1.0,fiberWidth)))
            truthPhoton.append((round(rand(1,10)), rand(0,25), randgauss(-1.0, fiberWidth), randgauss(-1.0,fiberWidth)))
            truthPhoton.append((round(rand(1,10)), rand(0,25), randgauss(-1.0, fiberWidth), randgauss(1.0,fiberWidth)))
            truthPhoton.append((round(rand(1,10)), rand(0,25), randgauss(1.0, fiberWidth), randgauss(-1.0,fiberWidth)))

        #############################################
        # Calculate adc per channel
        #############################################
        
        # Fill histos
        #for c in nChannels: histos["nPhotons_xy_{}".format(c[0])] = noise1D(histos["nPhotons_xy_{}".format(c[0])])
        
        for t in truthPhoton:
            for c in nChannels:
                histos["nPhotons_xy_{}".format(c[0])].Fill(t[2], t[3], t[0])

        for t in truthPhoton:
            for c in nChannels:
                if(hitOneArray(histos["nPhotons_xy_{}".format(c[0])], t)): histos["nPhotons_time_{}".format(c[0])].Fill(t[1], t[0])

        for c in nChannels: 
            histos["signal_time_{}".format(c[0])] = noise1D(histos["signal_time_{}".format(c[0])])
            histos["signal_time_{}".format(c[0])].Add(histos["nPhotons_time_{}".format(c[0])])

        # # Save to TTree
        # xtruth_[0] = xtruth
        # for x in ampTri1:
        #     amp1_.push_back(x)
        # for x in ampTri2:
        #     amp2_.push_back(x)

        # tree.Fill()
        # amp1_.clear()
        
    # Draw histos
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    for c in nChannels: saveHisto(c1, histos, "nPhotons_xy_{}".format(c[0]))
    for c in nChannels: saveHisto(c1, histos, "nPhotons_time_{}".format(c[0]),"hist")
    for c in nChannels: saveHisto(c1, histos, "signal_time_{}".format(c[0]),"hist")

    # Write everything out
    root_file.Write()
    tree.Write()
    root_file.Close()
        
if __name__ == '__main__':
    main()
