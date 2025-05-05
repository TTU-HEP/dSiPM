import random
import math
import numpy as np
import ROOT
import copy
import os
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()

xShift = [ 3.2908,  3.2205,  3.3597,  3.2898]
yShift = [-3.6878, -3.7287, -3.7284, -3.6075]

class Photons:
    def __init__(self, event):
        #Sort by arrival time of the photon
        self.time_final = np.array(event.OP_time_final)
        self.array4Sorting = np.argsort(self.time_final)
        self.time_final = self.time_final[self.array4Sorting]

        #Get needed variables
        self.productionFiber = self.getArrayFromEvent(event.OP_productionFiber)
        self.isCoreC = self.getArrayFromEvent(event.OP_isCoreC)
        self.pos_final_x = self.getArrayFromEvent(event.OP_pos_final_x)
        self.pos_final_y = self.getArrayFromEvent(event.OP_pos_final_y)
        self.pos_final_z = self.getArrayFromEvent(event.OP_pos_final_z)
        self.pos_produced_z = self.getArrayFromEvent(event.OP_pos_produced_z)
        self.w = np.full(self.nPhotons(), 1.0)

    def getArrayFromEvent(self, var):
        return np.array(var)[self.array4Sorting]

    def nPhotons(self):
        return len(self.pos_final_x)
    def x(self, i):
        return self.pos_final_x[i] + xShift[self.productionFiber[i]]
    def y(self, i):
        return self.pos_final_y[i] + yShift[self.productionFiber[i]]
    def z(self, i):
        return self.pos_produced_z[i]
    def zEnd(self, i):
        return self.pos_final_z[i]
    def t(self, i):
        return self.time_final[i]

class SiPMInfo:
    def __init__(self, channelSize, nBins):
        self.channelSize = channelSize #microns
        self.nBins = nBins
        self.name = "{0}x{0}".format(int(self.channelSize))

def saveHisto(c, hDic, name, drawOpt="", doLog=False):
    hDic[name].Draw(drawOpt)
    hDic[name].SetStats(0)
    if doLog:
        ROOT.gPad.SetLogy()
        #ROOT.gPad.SetLogx()
    else:
        ROOT.gPad.SetLogy(0)
    c.SaveAs("output/"+name+".png")    

def main():
    # Some useful hardcoded stuff
    input_file = ROOT.TFile("mc_dreamsim_e+_0_run100_0_Test_50evt_e+_100_101.root", "READ")
    tree = input_file.Get("tree")
    root_file = ROOT.TFile("outfile.root", "RECREATE")
    os.makedirs("output/", exist_ok=True)

    # Define histos
    histos = {}

    nChannels = [
        #SiPMInfo(   1, 3000), 
        #SiPMInfo(   2, 1500), 
        #SiPMInfo(   5,  600), 
        SiPMInfo(  10,  300), 
        #SiPMInfo(  20,  150), 
        SiPMInfo(  25,  120), 
        #SiPMInfo(  30,  100), 
        #SiPMInfo(  40,   75), 
        SiPMInfo(  50,   60), 
        #SiPMInfo(  60,   50), 
        SiPMInfo(  75,   40), 
        SiPMInfo( 100,   30), 
        #SiPMInfo( 120,   25), 
        #SiPMInfo( 125,   24), 
        SiPMInfo( 200,   15),
        SiPMInfo( 300,   10),
        SiPMInfo( 600,    5),
        #SiPMInfo(3000,    1),
    ]

    for c in nChannels: histos[        "nPhotons_xyt_all_{}".format(c.name)] = ROOT.TH3D(        "nPhotons_xyt_all_{}".format(c.name),"nPhotons_xyt; x [mm]; y [mm]; t [ns]; nPhotons", c.nBins    ,-0.5,0.5, c.nBins,-0.5,0.5, 700,5.0, 40.0)
    for c in nChannels: histos[        "nPhotons_rte_all_{}".format(c.name)] = ROOT.TH3D(        "nPhotons_rte_all_{}".format(c.name),"nPhotons_rt;  r [mm]; t [ns]; events; nPhotons", 60,0.0,0.5,    700,5.0, 40.0, 100,  0,  100)
    for c in nChannels: histos[        "nPhotons_rze_all_{}".format(c.name)] = ROOT.TH3D(        "nPhotons_rze_all_{}".format(c.name),"nPhotons_rz;  r [mm]; z [mm]; events; nPhotons", 60,0.0,0.5,    500,0.0,2000.0, 100,  0,  100)
    for c in nChannels: histos[        "nPhotons_tze_all_{}".format(c.name)] = ROOT.TH3D(        "nPhotons_tze_all_{}".format(c.name),"nPhotons_rz;  t [ns]; z [mm]; events; nPhotons", 700,5.0, 40.0, 500,0.0,2000.0, 100,  0,  100)
    for c in nChannels: histos[     "nPhotons_xyt_oneHit_{}".format(c.name)] = ROOT.TH3D(     "nPhotons_xyt_oneHit_{}".format(c.name),"nPhotons_xyt; x [mm]; y [mm]; t [ns]; nPhotons", c.nBins    ,-0.5,0.5, c.nBins,-0.5,0.5, 700,5.0, 40.0)
    for c in nChannels: histos[     "nPhotons_rte_oneHit_{}".format(c.name)] = ROOT.TH3D(     "nPhotons_rte_oneHit_{}".format(c.name),"nPhotons_rt;  r [mm]; t [ns]; events; nPhotons", 60,0.0,0.5,    700,5.0, 40.0, 100,  0,  100)
    for c in nChannels: histos[     "nPhotons_rze_oneHit_{}".format(c.name)] = ROOT.TH3D(     "nPhotons_rze_oneHit_{}".format(c.name),"nPhotons_rz;  r [mm]; z [mm]; events; nPhotons", 60,0.0,0.5,    500,0.0,2000.0, 100,  0,  100)
    for c in nChannels: histos[     "nPhotons_tze_oneHit_{}".format(c.name)] = ROOT.TH3D(     "nPhotons_tze_oneHit_{}".format(c.name),"nPhotons_rz;  t [ns]; z [mm]; events; nPhotons", 700,5.0, 40.0, 500,0.0,2000.0, 100,  0,  100)
    for c in nChannels: histos["temp_nPhotons_xyt_oneHit_{}".format(c.name)] = ROOT.TH2D("temp_nPhotons_xyt_oneHit_{}".format(c.name),"nPhotons_xyt; x [mm]; y [mm]; t [ns]; nPhotons", c.nBins    ,-0.5,0.5, c.nBins,-0.5,0.5)

    for c in nChannels: histos[       "dummy_nPhotons_xy_{}".format(c.name)] = ROOT.TH2D(       "dummy_nPhotons_xy_{}".format(c.name),"nPhotons_xy; x [cm]; y [cm]; nPhotons", c.nBins,-0.5,0.5, c.nBins,-0.5,0.5)
    for c in nChannels: histos[      "nPhotonsPerChannel_{}".format(c.name)] = ROOT.TH1D(      "nPhotonsPerChannel_{}".format(c.name),"nPhotonsChannel; nPhotons", 30, 0, 30)

    # Loop over events
    nEvents = 0
    for event in tree:
        nEvents += 1

        g = Photons(event)
        print("Total number of photons:", g.nPhotons())

        #Loop over photons in the event
        for i, _ in np.ndenumerate(g.pos_final_x):
            for c in nChannels:
                if g.zEnd(i)>0 and bool(g.isCoreC[i]) and g.t(i)>0.0 and g.t(i)<40.0:
                    x, y, z, t, w = 10*g.x(i), 10*g.y(i), 20*g.z(i) + 2000, g.t(i), g.w[i]
                    r = math.sqrt(x**2 + y**2)
                    histos["nPhotons_xyt_all_{}".format(c.name)].Fill(x, y, t, w)          
                    histos["nPhotons_rte_all_{}".format(c.name)].Fill(r, t, nEvents, w)          
                    histos["nPhotons_rze_all_{}".format(c.name)].Fill(r, z, nEvents, w)          
                    histos["nPhotons_tze_all_{}".format(c.name)].Fill(t, z, nEvents, w)          
                    histos["dummy_nPhotons_xy_{}".format(c.name)].Fill(x, y, w)

                    #Check if bin was hit
                    b = histos["temp_nPhotons_xyt_oneHit_{}".format(c.name)].FindBin(x, y)
                    nPhotonsInBin = histos["temp_nPhotons_xyt_oneHit_{}".format(c.name)].GetBinContent(b)
                    if nPhotonsInBin < 1:
                        histos["nPhotons_xyt_oneHit_{}".format(c.name)].Fill(x, y, t, w)
                        histos["nPhotons_rte_oneHit_{}".format(c.name)].Fill(r, t, nEvents, w)
                        histos["nPhotons_rze_oneHit_{}".format(c.name)].Fill(r, z, nEvents, w)
                        histos["nPhotons_tze_oneHit_{}".format(c.name)].Fill(t, z, nEvents, w)
                        histos["temp_nPhotons_xyt_oneHit_{}".format(c.name)].Fill(x, y, w)

        #Clear temp histos
        for c in nChannels: 
            # Fill the nPhotons per channel summary histogram
            h = histos["dummy_nPhotons_xy_{}".format(c.name)]
            for bx in range(1, h.GetXaxis().GetNbins()+1):
                for by in range(1, h.GetYaxis().GetNbins()+1):
                    nPhotonsPerChannel = h.GetBinContent(bx, by)
                    histos["nPhotonsPerChannel_{}".format(c.name)].Fill(nPhotonsPerChannel)
            histos["dummy_nPhotons_xy_{}".format(c.name)].Reset()
            histos["temp_nPhotons_xyt_oneHit_{}".format(c.name)].Reset()

    # Draw and save histos
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)
    # for c in nChannels: saveHisto(c1, histos, "nPhotons_xy_{}".format(c.name), "colz")

    # Write everything out
    root_file.Write()
    #tree.Write()
    root_file.Close()
        
if __name__ == '__main__':
    main()
