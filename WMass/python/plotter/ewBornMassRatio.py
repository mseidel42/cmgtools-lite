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

def makeRatio(boson, name1, name2, obs='rap'):
    print('making plot for %s: %s/%s' % (boson, name1, name2))
    
    legend = ROOT.TLegend(0.5,0.6,0.95,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    
    tfile1 = ROOT.TFile('plots/ewBornMass/%s/mass_rap_%s.root' % (boson, name1))
    canvas1 = tfile1.Get('c')
    hist1 = canvas1.GetPrimitive(name1+'_'+obs)

    hist1.Scale(1./hist1.Integral())
    
    colors = [0, ROOT.kRed+1, ROOT.kOrange+1, ROOT.kGreen+1, ROOT.kCyan+1, ROOT.kAzure+1,  ROOT.kBlue+1, ROOT.kViolet+1, ROOT.kMagenta+1]
    
    projs1 = []
    ybins = [1, 3, 5, 7]
    for ybin in ybins:
        low = hist1.GetYaxis().GetBinLowEdge(ybin)
        high = hist1.GetYaxis().GetBinUpEdge(ybin+1)
        if low >= 4.:
            break
        proj = hist1.ProjectionX(name1+'_'+obs+str(ybin), ybin, ybin+1)
        proj.SetLineColor(colors[ybin])
        proj.GetYaxis().SetRangeUser(0.5, 1.3)
        proj.GetYaxis().SetTitle('%s / %s' % (name1, name2))
        proj.Scale(1./proj.Integral())
        projs1.append(proj)
        legend.AddEntry(proj, '%.2f < |Y| < %.2f' % (low, high), "l")
    
    tfile2 = ROOT.TFile('plots/ewBornMass/%s/mass_rap_%s.root' % (boson, name2))
    canvas2 = tfile2.Get('c')
    hist2 = canvas2.GetPrimitive(name2+'_'+obs)

    hist2.Scale(1./hist2.Integral())
    
    projs2 = []
    for ybin in ybins:
        proj = hist2.ProjectionX(name2+'_'+obs+str(ybin), ybin, ybin)
        proj.SetLineColor(1+ybin)
        proj.Scale(1./proj.Integral())
        projs2.append(proj)

    print(projs1)
    
    hist = hist1.Clone('ratio')
    hist.Divide(hist1, hist2)
    
    ROOT.gStyle.SetPalette(ROOT.kThermometer)
    ROOT.gStyle.SetNumberContours(100)

    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.12)
    c.SetLeftMargin(0.15)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.cd()

    hist.GetZaxis().SetRangeUser(0.5, 1.5)
    hist.SetContour(100)
    hist.Draw('colz')

    c.Print('plots/ewBornMass/%s/mass_rap_ratio_%s_%s.pdf'  % (boson, name1, name2))
    c.Print('plots/ewBornMass/%s/mass_rap_ratio_%s_%s.png'  % (boson, name1, name2))
    # c.Print('plots/ewBornMass/%s/mass_rap_ratio__%s_%s.root' % (args.boson, name1, name2))
    rootFile = ROOT.TFile('plots/ewBornMass/%s/mass_rap_ratio_%s_%s.root' % (boson, name1, name2), 'RECREATE')
    hist.Write()
    rootFile.Close()
    
    c.SetRightMargin(0.05)
    
    for i in range(len(projs1)):
        projs1[i].Divide(projs2[i])
        projs1[i].GetXaxis().SetRangeUser(87, 95)
        projs1[i].GetYaxis().SetRangeUser(0.95, 1.10)
        if i == 0:
            projs1[i].Draw()
        else:
            projs1[i].Draw('same')
    legend.Draw()
    
    c.Print('plots/ewBornMass/%s/mass_rap_ratio1d_%s_%s.pdf'  % (boson, name1, name2))
    c.Print('plots/ewBornMass/%s/mass_rap_ratio1d_%s_%s.png'  % (boson, name1, name2))
    
    legend = ROOT.TLegend(0.5,0.6,0.95,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    
    hist1rap = hist1.ProjectionY(name1+'_rap')
    hist2rap = hist2.ProjectionY(name2+'_rap')
    hist1rap.SetLineColor(ROOT.kRed+1)
    hist2rap.SetLineColor(ROOT.kBlue+1)
    legend.AddEntry(hist1rap, name1, "l")
    legend.AddEntry(hist2rap, name2, "l")
    
    hist1rap.Draw()
    hist2rap.Draw('same')
    legend.Draw()
    
    c.Print('plots/ewBornMass/%s/rap_%s_%s.pdf'  % (boson, name1, name2))
    c.Print('plots/ewBornMass/%s/rap_%s_%s.png'  % (boson, name1, name2))


if args.boson == 'all':
    bosons = ['Z', 'Wplus', 'Wminus']
else:
    bosons = [args.boson]

hists = {}

for boson in bosons:
    if args.input == 'all':
        makeRatio(boson, 'minnlo', 'horace-photos')
        makeRatio(boson, 'horace-exp', 'horace-photos')
        makeRatio(boson, 'horace-exp-old', 'horace-photos')
    else:
        name1,name2 = tuple(args.input.split('/'))
        makeRatio(boson, name1, name2)