#!/usr/bin/env python
from math import *
import re
import os, os.path
from array import array
import logging
from collections import OrderedDict
import numbaDefines
from math import sqrt

## safe batch mode
import sys
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gInterpreter.ProcessLine(".O3")
ROOT.EnableImplicitMT()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",  type=str, default='horace-photos/minnlo', help="File to read histos from")
parser.add_argument("-b", "--boson", type=str, default='Z', help="boson")
parser.add_argument("-w", "--withISR", type=bool, default=True, help="Consider ISR photons")
args = parser.parse_args()

ROOT.gStyle.SetOptStat(0)

def makeRatio(boson, name1, name2):
    print('making plot for %s: %s/%s' % (boson, name1, name2))

    hist = hists[name1].Clone('weights')
    hist.Divide(hists[name1], hists[name2])
    
    ROOT.gStyle.SetPalette(ROOT.kThermometer)
    ROOT.gStyle.SetNumberContours(100)
    
    weightHistOrig = ROOT.TH1D('weightHistOrig', ';weight;bins', 60, 0, 3)
    weightHistSmooth = ROOT.TH1D('weightHistSmooth', ';weight;bins', 60, 0, 3)
    
    # smoothing
    for i in range(hist.GetNbinsX()+2):
        for j in range(hist.GetNbinsY()+2):
            val = hist.GetBinContent(i, j)
            if val == 0.:
                hist.SetBinContent(i, j, 1.)
                continue
            errWeight = 2. # approach 1 when error is 50%
            relErr = min(1., errWeight * hist.GetBinError(i, j) / val)
            ipolVal = (1.-relErr) * val + relErr
            hist.SetBinContent(i, j, ipolVal)
            weightHistOrig.Fill(val)
            weightHistSmooth.Fill(ipolVal)

    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.12)
    c.SetLeftMargin(0.15)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.cd()

    hist.GetZaxis().SetRangeUser(0., 2.)
    hist.SetContour(100)
    hist.Draw('colz')

    c.Print('plots/ewPhotonKinematics/%s/ewPKRatio_%s_%s.pdf'  % (boson, name1, name2))
    c.Print('plots/ewPhotonKinematics/%s/ewPKRatio_%s_%s.png'  % (boson, name1, name2))
    # c.Print('plots/ewPhotonKinematics/%s/ewPKRatio_%s_%s.root' % (args.boson, name1, name2))
    rootFile = ROOT.TFile('plots/ewPhotonKinematics/%s/ewPKRatio_%s_%s.root' % (boson, name1, name2), 'RECREATE')
    hist.Write()
    rootFile.Close()
    
    # Plot error histogram
    ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator)
    hist_err = hist.Clone('weights_err')
    for i in range(0, hist_err.GetNbinsX()+2):
        for j in range(0, hist_err.GetNbinsY()+2):
            hist_err.SetBinContent(i, j, hist.GetBinError(i,j))
    c.SetLogz()
    hist_err.GetZaxis().SetRangeUser(1e-3, 10)
    hist_err.Draw('colz')
    
    c.Print('plots/ewPhotonKinematics/%s/ewPKRatio_%s_%s_err.pdf'  % (boson, name1, name2))
    c.Print('plots/ewPhotonKinematics/%s/ewPKRatio_%s_%s_err.png'  % (boson, name1, name2))
    c.Print('plots/ewPhotonKinematics/%s/ewPKRatio_%s_%s_err.root' % (boson, name1, name2))
    
    # Plot weight histogram
    c.SetRightMargin(0.05)
    weightHistOrig.SetLineColor(ROOT.kRed+1)
    weightHistOrig.SetLineStyle(7)
    weightHistSmooth.SetLineColor(ROOT.kBlue+1)
    
    weightHistOrig.Draw()
    weightHistSmooth.Draw('same')
    
    legend = ROOT.TLegend(0.5,0.6,0.95,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    legend.AddEntry(weightHistOrig, 'Original weights', "l")
    legend.AddEntry(weightHistSmooth, 'Smoothed weights', "l")
    legend.Draw()
    
    c.Print('plots/ewPhotonKinematics/%s/ewPKRatio_%s_%s_1d.png'  % (boson, name1, name2))
    c.Print('plots/ewPhotonKinematics/%s/ewPKRatio_%s_%s_1d.root' % (boson, name1, name2))
    c.Print('plots/ewPhotonKinematics/%s/ewPKRatio_%s_%s_1d.pdf'  % (boson, name1, name2))

def getHist(boson, name):
    tfile = ROOT.TFile('plots/ewPhotonKinematics/%s/ewPK_%s.root' % (boson, name))
    canvas = tfile.Get('c')
    hist = canvas.GetPrimitive(name)

    # merge no/low-emission bins
    maxBinLE = hist.GetYaxis().FindFixBin(-2.)
    
    # integralLE = hist.Integral(11, maxBinLE, 11, maxBinLE)
    for i in range(0, hist.GetNbinsX()+2):
        integralLE = 0.
        integralLEerr = 0.
        for j in range(0, maxBinLE):
            integralLE += hist.GetBinContent(i, j)
            integralLEerr += hist.GetBinError(i, j)
        for j in range(0, maxBinLE):
            hist.SetBinContent(i, j, integralLE)
            hist.SetBinError(i, j, integralLEerr)

    hist.Scale(1./hist.Integral())
    return hist

if args.boson == 'all':
    bosons = ['Z', 'Wplus', 'Wminus']
else:
    bosons = [args.boson]

hists = {}

for boson in bosons:
    if args.input == 'all':
        allHists = ['minnlo', 'horace-photos', 'horace-exp', 'horace-exp-old']
        for name in allHists:
            hists[name] = getHist(boson, name + '_withISR')
        makeRatio(boson, 'minnlo', 'horace-photos')
        makeRatio(boson, 'horace-exp', 'horace-photos')
        makeRatio(boson, 'horace-exp-old', 'horace-photos')
    else:
        name1,name2 = tuple(args.input.split('/'))
        for name in [name1, name2]:
            hists[name] = getHist(boson, name + '_withISR')
        makeRatio(boson, name1, name2)