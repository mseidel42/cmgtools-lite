#!/usr/bin/env python
from math import *
import re
import os, os.path
from array import array
import logging
from collections import OrderedDict
import numbaDefines
import numpy as np
from glob import glob
from tqdm import tqdm

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gInterpreter.ProcessLine(".O3")
ROOT.EnableImplicitMT()

from copy import *
from cutsFile import *
from fakeRate import *
from mcCorrections import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",  type=str, default='minnlo', help="File to read histos from")
parser.add_argument("-b", "--boson", type=str, default='Z', help="boson")
parser.add_argument("-m", "--maxFiles", type=int, default=50, help="Maximum number of files")
parser.add_argument('--snapshot', action='store_true')
parser.add_argument('--binning', action='store_true')
parser.add_argument('-z', '--deleteZombies', action='store_true')
args = parser.parse_args()

def makeBinning(mass = 91.1535, width = 2.4932, initialStep = 0.1):
    maxVal = ROOT.Math.breitwigner_pdf(mass, width, mass)
    bins = [mass]
    currentMass = mass
    while currentMass - mass < 100:
        binSize = maxVal / ROOT.Math.breitwigner_pdf(currentMass, width, mass) * initialStep
        currentMass += binSize
        bins.append(currentMass)
        lowMass = 2*mass - currentMass
        if lowMass - binSize > 0:
            bins.insert(0, lowMass)
    bins.insert(0, 0.)
    return bins

def normalizeAndDivideByBinWidth(hist):
    hist.Scale(1./hist.Integral())
    for i in range(1, hist.GetNbinsX()+1):
        for j in range(1, hist.GetNbinsY()+1):
            hist.SetBinContent(i, j, hist.GetBinContent(i,j) / hist.GetXaxis().GetBinWidth(i))
            hist.SetBinError  (i, j, hist.GetBinError(i,j)   / hist.GetXaxis().GetBinWidth(i))
    return hist

if args.binning:
    print(makeBinning())
    sys.exit()

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/functions.cc")

# ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kViridis)
ROOT.TColor.InvertPalette()

def makePlot(name):
    print('making plot for ' + name)

    chain = ROOT.TChain('Events')

    path = paths[args.boson][name]
    files = glob(path+'/*.root')
    nFilesAdded = 0
    print('Adding files')
    for f in tqdm(files):
        if args.maxFiles > 0 and nFilesAdded > args.maxFiles:
            break
        if '.root' in f:
            markedZombie = False
            if args.deleteZombies and not 'powhegMiNNLO' in path:
                try:
                    thisFile = ROOT.TFile(f)
                    thisFile.Close()
                except OSError:
                    os.remove(f)
                    markedZombie = True
            if not markedZombie:
                chain.Add(f)
                nFilesAdded += 1
    print('Added %i files' % nFilesAdded)

    rdf = ROOT.RDataFrame(chain)
    rdf = rdf.Define('genPrefsrLeps', 'Numba::prefsrLeptons(GenPart_status, GenPart_statusFlags, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt)')
    rdf = rdf.Define('ptVgen', 'transversemomentum(GenPart_pt[genPrefsrLeps],GenPart_phi[genPrefsrLeps])')
    rdf = rdf.Define('rapVgen', 'abs(rapidity(GenPart_pt[genPrefsrLeps],GenPart_eta[genPrefsrLeps],GenPart_phi[genPrefsrLeps],GenPart_mass[genPrefsrLeps]))')
    rdf = rdf.Define('massVgen', 'invariantmass(GenPart_pt[genPrefsrLeps],GenPart_eta[genPrefsrLeps],GenPart_phi[genPrefsrLeps],GenPart_mass[genPrefsrLeps])')

    if args.snapshot:
        rdf.Snapshot('Events', 'snapshot_%s.root' % name, ['ptVgen', 'ewSel', 'nGenPart', 'GenPart_pdgId', 'GenPart_status', 'GenPart_statusFlags', 'GenPart_genPartIdxMother', 'GenPart_pt', 'GenPart_eta', 'GenPart_phi', 'ptVgen', 'rapVgen', 'massVgen'])

    if args.boson == 'Z':
        massBins = makeBinning(mass = 91.1535, width = 2.4932, initialStep=0.10)
    else:
        massBins = makeBinning(mass = 80.3815, width = 2.0904, initialStep=0.10)
    ptBins = np.logspace(-1,3,9)
    ptModel = ROOT.RDF.TH2DModel(name+'_pt', ';Boson mass [GeV];Boson p_{T} [GeV]', len(massBins)-1, array('d', massBins), len(ptBins)-1, array('d', ptBins))
    ptHist = rdf.Histo2D(ptModel, 'massVgen', 'ptVgen')
    
    rapBins = np.linspace(0, 5, 11)
    rapModel = ROOT.RDF.TH2DModel(name+'_rap', ';Boson mass [GeV];Boson rapidity', len(massBins)-1, array('d', massBins), len(rapBins)-1, array('d', rapBins))
    rapHist = rdf.Histo2D(rapModel, 'massVgen', 'rapVgen')

    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.12)
    c.SetLeftMargin(0.15)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.SetLogz()
    c.cd()
    if not os.path.isdir('plots/ewBornMass/%s/' % args.boson):
        os.makedirs('plots/ewBornMass/%s/' % args.boson)

    c.SetLogy(True)
    ptHist.Draw('colz')
    c.Update()
    
    st = ptHist.FindObject('stats')
    st.SetX1NDC(0.60)
    st.SetX2NDC(0.87)
    c.Modified()
    c.Update()

    c.Print('plots/ewBornMass/%s/mass_pt_%s.pdf'  % (args.boson, name))
    c.Print('plots/ewBornMass/%s/mass_pt_%s.png'  % (args.boson, name))
    c.Print('plots/ewBornMass/%s/mass_pt_%s.root' % (args.boson, name))
    
    c.SetLogy(False)
    rapHist.Draw('colz')
    c.Update()
    
    st = rapHist.FindObject('stats')
    st.SetX1NDC(0.60)
    st.SetX2NDC(0.87)
    c.Modified()
    c.Update()

    c.Print('plots/ewBornMass/%s/mass_rap_%s.pdf'  % (args.boson, name))
    c.Print('plots/ewBornMass/%s/mass_rap_%s.png'  % (args.boson, name))
    c.Print('plots/ewBornMass/%s/mass_rap_%s.root' % (args.boson, name))

baseNano = '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoAOD/'
baseNanoGen = '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/'
paths = {}
paths['Z'] = {
    'minnlo': baseNano+'DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MC*VFP/210319_*/*/',
    'horace-photos': baseNanoGen+'ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photos-isr-pythia/',
}
paths['Wplus'] = {
    'minnlo': baseNano+'WplusJetsToMuNu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/MC*VFPWeightFix/210701_234251/*/',
    'horace-photos': baseNanoGen+'WplusToMuNu_TuneCP5_13TeV-horace-born-fsr-photos-isr-pythia',
}
paths['Wminus'] = {
    'minnlo': baseNano+'WminusJetsToMuNu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/MC*VFPWeightFix/210701_234131/*/',
    'horace-photos': baseNanoGen+'WminusToMuNu_TuneCP5_13TeV-horace-born-fsr-photos-isr-pythia',
}


if args.input == 'all':
    for pathname in paths[args.boson]:
        makePlot(pathname)
else:
    makePlot(args.input)