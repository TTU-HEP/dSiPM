import random
import math
import numpy as np
import ROOT
import copy
import os
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()

class Photons:
    def __init__(self, event):
        self.isCoreC = np.array(event.OP_isCoreC)
        self.pos_final_x = np.array(event.OP_pos_final_x)
        self.pos_final_y = np.array(event.OP_pos_final_y)
        self.pos_final_z = np.array(event.OP_pos_final_z)
        self.time_final = np.array(event.OP_time_final)

    def nPhotons(self):
        return len(self.pos_final_x)

def saveHisto(c, hDic, name, drawOpt="", doLog=False):
    hDic[name].Draw(drawOpt)
    hDic[name].SetStats(0)
    if doLog:
        ROOT.gPad.SetLogy()
        #ROOT.gPad.SetLogx()
    else:
        ROOT.gPad.SetLogy(0)
    c.SaveAs("output/"+name+".png")    

def rand(m=0,M=100):
    a = random.uniform(m,M)
    return a
    
def randgauss(m,w):
    return random.gauss(m,w)
    
def noise1D(histo):
    for b in range(1, histo.GetXaxis().GetNbins()+1):
        histo.SetBinContent(b, randgauss(0, 0.5))
    return histo

def hitOneArray(h, x, y):
    binNum = h.GetMaximumBin()
    b = h.FindBin(x, y)
    return binNum==b

def getSiPMResponse(time):
    l = ROOT.TF1("landau","landau(0)+gaus(3)",0, 25)
    l.SetParameter(0, 0.5)
    l.SetParameter(1, time)
    l.SetParameter(2, 0.2)
    l.SetParameter(3, 1)
    l.SetParameter(4, time)
    l.SetParameter(5, 0.25)
    l.SetNpx(500)

    h = ROOT.TH1D(l.GetHistogram())
    h.SetName("{}".format(time))
    if h.Integral()>0:
        h.Scale(1.0/(0.9*h.GetMaximum()))
    return h

def main():
    # Some useful hardcoded stuff
    #input_file = ROOT.TFile("mc_testjob_run001_003_Test_20evt_pi+_100.0_100.0.root", "READ")
    input_file = ROOT.TFile("mc_testjob_run001_003_Test_20evt_e+_100.0_100.0.root", "READ")
    tree = input_file.Get("tree")
    root_file = ROOT.TFile("outfile.root", "RECREATE")
    os.makedirs("output/", exist_ok=True)

    # Define histos
    histos = {}

    nChannels = [
        ("1x1",     3000), 
        #("2x2",     1500), 
        #("5x5",      600), 
        ("10x10",    300), 
        #("20x20",    150), 
        ("25x25",    120), 
        #("30x30",    100), 
        #("40x40",     75), 
        ("50x50",     60), 
        #("60x60",     50), 
        ("75x75",     40), 
        ("100x100",   30), 
        #("120x120",   25), 
        #("125x125",   24), 
        #("200x200",   15),
        #("600x600",    5),
        ("3000x3000",  1),
    ]

    for c in nChannels: histos["nPhotons_xy_{}".format(c[0])] = ROOT.TH2D("nPhotons_xy_{}".format(c[0]),"nPhotons_xy; x [cm]; y [cm]; nPhotons", c[1],-3.45,-3.15, c[1],3.5,3.8)
    for c in nChannels: histos["dummy_nPhotons_xy_{}".format(c[0])] = ROOT.TH2D("dummy_nPhotons_xy_{}".format(c[0]),"nPhotons_xy; x [cm]; y [cm]; nPhotons", c[1],-3.45,-3.15, c[1],3.5,3.8)
    for c in nChannels: histos["nPhotonsPerChannel_{}".format(c[0])] = ROOT.TH1D("nPhotonsPerChannel_{}".format(c[0]),"nPhotonsChannel; nPhotons", 20, 0, 20)
    for c in nChannels: histos["nPhotonsPerChannel_time_{}".format(c[0])] = ROOT.TH2D("nPhotonsPerChannel_time_{}".format(c[0]),"nPhotonsChannel; nPhotons; time", 500,0.0,25.0, 20,0,20)
    for c in nChannels: histos["nPhotons_time_{}".format(c[0])] = ROOT.TH1D("nPhotons_time_{}".format(c[0]),"nPhotons_time; time [ns]; nPhotons", 500,0.0,25.0)
    for c in nChannels: histos["nPhotons_time_smear_{}".format(c[0])] = ROOT.TH1D("nPhotons_time_smear_{}".format(c[0]),"nPhotons_time_smear; time [ns]; nPhotons Smear", 500,0.0,25.0)
    for c in nChannels: histos["signal_time_{}".format(c[0])] = ROOT.TH1D("signal_time_{}".format(c[0]),"signal_time; time [ns]; Amplitude [mV]", 500,0.0,25.0)

    # Loop over events
    nEvents = 0
    for event in tree:
        nEvents += 1
        if nEvents < 7:
            continue

        if nEvents >= 8: break
        # break

        gammas = Photons(event)
        print("Total number of photons:", gammas.nPhotons())

        #Loop over photons in the event
        for i, _ in np.ndenumerate(gammas.pos_final_x):
            x = gammas.pos_final_x[i]
            y = gammas.pos_final_y[i]
            z = gammas.pos_final_z[i]
            t = gammas.time_final[i]
            isCoreC = bool(gammas.isCoreC[i])
            for c in nChannels:
                if z>0 and isCoreC:
                    histos["nPhotons_xy_{}".format(c[0])].Fill(x, y)
                    histos["dummy_nPhotons_xy_{}".format(c[0])].Fill(x, y)

        for i, _ in np.ndenumerate(gammas.pos_final_x):
            x = gammas.pos_final_x[i]
            y = gammas.pos_final_y[i]
            z = gammas.pos_final_z[i]
            t = gammas.time_final[i]
            isCoreC = bool(gammas.isCoreC[i])
            for c in nChannels:
                if z>0 and isCoreC:
                    if(hitOneArray(histos["nPhotons_xy_{}".format(c[0])], x, y)): 
                        histos["nPhotons_time_{}".format(c[0])].Fill(t)
        #                 rSiPM = getSiPMResponse(t)
        #                 histos["nPhotons_time_smear_{}".format(c[0])].Add(rSiPM)

        for c in nChannels: 
            # #Make realistic waveform plots
            # histos["signal_time_{}".format(c[0])] = noise1D(histos["signal_time_{}".format(c[0])])
            # histos["signal_time_{}".format(c[0])].Add(histos["nPhotons_time_{}".format(c[0])])
            # temp = copy.deepcopy(histos["nPhotons_time_smear_{}".format(c[0])])
            # temp.Scale(8)
            # histos["signal_time_{}".format(c[0])].Add(temp)

            # Fill the nPhotons per channel summary histogram
            h = histos["dummy_nPhotons_xy_{}".format(c[0])]
            #print(h.GetEntries())
            for bx in range(1, h.GetXaxis().GetNbins()+1):
                for by in range(1, h.GetYaxis().GetNbins()+1):
                    nPhotonsPerChannel = h.GetBinContent(bx, by)
                    #h.SetBinContent(bx, by, 0)
                    #if nPhotonsPerChannel == 0: continue
                    histos["nPhotonsPerChannel_{}".format(c[0])].Fill(nPhotonsPerChannel)
            histos["dummy_nPhotons_xy_{}".format(c[0])].Reset()
            
    # Draw and save histos
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)
    for c in nChannels: saveHisto(c1, histos, "nPhotons_xy_{}".format(c[0]), "colz")
    for c in nChannels: saveHisto(c1, histos, "nPhotonsPerChannel_{}".format(c[0]),"hist", True)
    for c in nChannels: saveHisto(c1, histos, "nPhotons_time_{}".format(c[0]),"hist")
    #for c in nChannels: saveHisto(c1, histos, "nPhotons_time_smear_{}".format(c[0]),"hist")
    #for c in nChannels: saveHisto(c1, histos, "signal_time_{}".format(c[0]),"hist")

    # Write everything out
    root_file.Write()
    tree.Write()
    root_file.Close()
        
if __name__ == '__main__':
    main()
