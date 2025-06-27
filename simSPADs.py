import sys
import random
import math
import numpy as np
import ROOT
import copy
import os

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()

# xShift = [ 3.2908,  3.2205,  3.3597,  3.2898]
# yShift = [-3.6878, -3.7287, -3.7284, -3.6075]

xShift = [ 3.7308,  3.5768,  3.8008,  3.6468]
yShift = [-3.7293, -3.6878, -3.6878, -3.6488]

# (distance_limit_from_center, shift_amount_toward_center)
shrink_rules = [
    (0.1,  0.00),
    (0.5,  0.23),
    (0.9,  0.46),
    (1.3,  0.69),
    (1.7,  0.92),
    (2.1,  1.15),
    (2.5,  1.38),
    (2.9,  1.61),
    (3.3,  1.84),
    (3.7,  2.07),
    (4.1,  2.30),
    (4.5,  2.53),
    (4.9,  2.76),
    (5.3,  2.99),
    (5.7,  3.22),
    (6.1,  3.45),
    (6.5,  3.68),
    (6.9,  3.91),
    (7.3,  4.14),
    (7.7,  4.37),
    (8.1,  4.60),
    (8.5,  4.83),
    (8.9,  5.06),
    (9.3,  5.29),
    (9.7,  5.52),
    (10.1, 5.75),
    (10.5, 5.98),
    (10.9, 6.21),
    (11.3, 6.44),
    (11.7, 6.67),
    (12.1, 6.90),
    (12.5, 7.13),
    (12.9, 7.36),
    (13.3, 7.59),
    (13.7, 7.82),
    (14.1, 8.05),
    (14.5, 8.28),
    (14.9, 8.51),
    (15.3, 8.74),
    (15.7, 8.97),
]

def shrink_toward_center(val: float) -> float:
    """
    Move 'val' toward zero according to the graduated shrink_rules table.
    """
    abs_val = abs(val)
    for limit, shift in shrink_rules:
        if abs_val <= limit:
            return val - shift * np.sign(val)
    # Anything farther than the last limit: apply a constant max-shift (2.0 here)
    return val - 2.0 * np.sign(val)

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
        raw_x = self.pos_final_x[i] + xShift[self.productionFiber[i]]
        return shrink_toward_center(raw_x)

    def y(self, i):
        raw_y = self.pos_final_y[i] + yShift[self.productionFiber[i]]
        return shrink_toward_center(raw_y)

    def z(self, i):
        return self.pos_produced_z[i]
    def zEnd(self, i):
        return self.pos_final_z[i]
    def t(self, i):
        return self.time_final[i]
    def fiberNumber(self, i):
        return self.productionFiber[i]

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

def getNBins(l,h,s):
    return int((h - l)/s)

def main():
    # Some useful hardcoded stuff
    if len(sys.argv) < 2:
       print("❌ Error: please provide a ROOT file as argument.")
       sys.exit(1)
    input_file_path = sys.argv[1]
    input_file = ROOT.TFile(input_file_path, "READ")
    tree = input_file.Get("tree")
    root_file = ROOT.TFile("outfile.root", "RECREATE")
    os.makedirs("output/", exist_ok=True)

    # Define histos
    histos = {}

    xBinL = -40.0
    xBinH =  40.0

    print(getNBins(xBinL,xBinH,0.01))

    nChannels = [
        #SiPMInfo(   1, getNBins(xBinL,xBinH,0.001), 
        #SiPMInfo(  10, getNBins(xBinL,xBinH,0.010)), 
        SiPMInfo(  20,  getNBins(xBinL,xBinH,0.020)), 
        SiPMInfo(  30,  getNBins(xBinL,xBinH,0.030)), 
        SiPMInfo(  40,  getNBins(xBinL,xBinH,0.040)), 
        SiPMInfo(  50,  getNBins(xBinL,xBinH,0.050)), 
        SiPMInfo(  60,  getNBins(xBinL,xBinH,0.060)), 
        SiPMInfo(  70,  getNBins(xBinL,xBinH,0.070)), 
        SiPMInfo( 100,  getNBins(xBinL,xBinH,0.100)), 
        #SiPMInfo( 120,   25), 
        #SiPMInfo( 125,   24), 
        #SiPMInfo( 200,   15),
        #SiPMInfo( 300,   10),
        #SiPMInfo( 600,    5),
        #SiPMInfo(3000,    1),
    ]

    for c in nChannels: histos[        "nPhotons_xyt_all_{}".format(c.name)] = ROOT.TH3D(        "nPhotons_xyt_all_{}".format(c.name),"nPhotons_xyt; y [mm]; x [mm]; t [ns]; nPhotons", c.nBins    ,xBinL,xBinH, c.nBins,xBinL,xBinH, 1,5.0, 40.0)
    for c in nChannels: histos[        "nPhotons_rte_all_{}".format(c.name)] = ROOT.TH3D(        "nPhotons_rte_all_{}".format(c.name),"nPhotons_rt;  r [mm]; t [ns]; events; nPhotons", 60,0.0,0.5,    700,5.0, 40.0, 100,  0,  100)
    for c in nChannels: histos[        "nPhotons_rze_all_{}".format(c.name)] = ROOT.TH3D(        "nPhotons_rze_all_{}".format(c.name),"nPhotons_rz;  r [mm]; z [mm]; events; nPhotons", 60,0.0,0.5,    500,0.0,2000.0, 100,  0,  100)
    for c in nChannels: histos[        "nPhotons_tze_all_{}".format(c.name)] = ROOT.TH3D(        "nPhotons_tze_all_{}".format(c.name),"nPhotons_rz;  t [ns]; z [mm]; events; nPhotons", 700,5.0, 40.0, 500,0.0,2000.0, 100,  0,  100)
    for c in nChannels: histos[     "nPhotons_xyt_oneHit_{}".format(c.name)] = ROOT.TH3D(     "nPhotons_xyt_oneHit_{}".format(c.name),"nPhotons_xyt; y [mm]; x [mm]; t [ns]; nPhotons", c.nBins    ,xBinL,xBinH, c.nBins,xBinL,xBinH, 1,5.0, 40.0)
    for c in nChannels: histos[     "nPhotons_rte_oneHit_{}".format(c.name)] = ROOT.TH3D(     "nPhotons_rte_oneHit_{}".format(c.name),"nPhotons_rt;  r [mm]; t [ns]; events; nPhotons", 60,0.0,0.5,    700,5.0, 40.0, 100,  0,  100)
    for c in nChannels: histos[     "nPhotons_rze_oneHit_{}".format(c.name)] = ROOT.TH3D(     "nPhotons_rze_oneHit_{}".format(c.name),"nPhotons_rz;  r [mm]; z [mm]; events; nPhotons", 60,0.0,0.5,    500,0.0,2000.0, 100,  0,  100)
    for c in nChannels: histos[     "nPhotons_tze_oneHit_{}".format(c.name)] = ROOT.TH3D(     "nPhotons_tze_oneHit_{}".format(c.name),"nPhotons_rz;  t [ns]; z [mm]; events; nPhotons", 700,5.0, 40.0, 500,0.0,2000.0, 100,  0,  100)
    for c in nChannels: histos["temp_nPhotons_xyt_oneHit_{}".format(c.name)] = ROOT.TH2D("temp_nPhotons_xyt_oneHit_{}".format(c.name),"nPhotons_xyt; y [mm]; x [mm]; t [ns]; nPhotons", c.nBins    ,xBinL,xBinH, c.nBins,xBinL,xBinH)

    for c in nChannels: histos[       "dummy_nPhotons_xy_{}".format(c.name)] = ROOT.TH2D(       "dummy_nPhotons_xy_{}".format(c.name),"nPhotons_xy; y [cm]; x [cm]; nPhotons", c.nBins,xBinL,xBinH, c.nBins,xBinL,xBinH)
    for c in nChannels: histos[      "nPhotonsPerChannel_{}".format(c.name)] = ROOT.TH1D(      "nPhotonsPerChannel_{}".format(c.name),"nPhotonsChannel; nPhotons", 30, 0, 30)

    # Loop over events
    nEvents = -1
    for event in tree:
        nEvents += 1

        g = Photons(event)
        if nEvents % 5 == 0:
            print("Event number: {0} Total number of photons: {1}".format(nEvents, g.nPhotons()))

        #Loop over photons in the event
        nP = 0
        for i, _ in np.ndenumerate(g.pos_final_x):
            x, y, z, t, w = 10*g.x(i), 10*g.y(i), 20*g.z(i) + 2000, g.t(i), g.w[i]
            r = math.sqrt(x**2 + y**2)
            isCFiber = bool(g.isCoreC[i])
            #isCFiber = bool(g.isCoreC[i])
            isGoodPhoton = g.zEnd(i)>0 and g.t(i)>0.0 and g.t(i)<40.0
            if isCFiber and isGoodPhoton:
                nP += 1
                for c in nChannels:
                    histos["nPhotons_xyt_all_{}".format(c.name)].Fill(y, x, t, w)          
                    histos["nPhotons_rte_all_{}".format(c.name)].Fill(r, t, nEvents, w)          
                    histos["nPhotons_rze_all_{}".format(c.name)].Fill(r, z, nEvents, w)          
                    histos["nPhotons_tze_all_{}".format(c.name)].Fill(t, z, nEvents, w)          
                    histos["dummy_nPhotons_xy_{}".format(c.name)].Fill(y, x, w)

                    #Check if bin was hit
                    b = histos["temp_nPhotons_xyt_oneHit_{}".format(c.name)].FindBin(y, x)
                    nPhotonsInBin = histos["temp_nPhotons_xyt_oneHit_{}".format(c.name)].GetBinContent(b)
                    if nPhotonsInBin < 1:
                        histos["nPhotons_xyt_oneHit_{}".format(c.name)].Fill(y, x, t, w)
                        histos["nPhotons_rte_oneHit_{}".format(c.name)].Fill(r, t, nEvents, w)
                        histos["nPhotons_rze_oneHit_{}".format(c.name)].Fill(r, z, nEvents, w)
                        histos["nPhotons_tze_oneHit_{}".format(c.name)].Fill(t, z, nEvents, w)
                        histos["temp_nPhotons_xyt_oneHit_{}".format(c.name)].Fill(y, x, w)
        #print(nP)
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
