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
parser.add_argument("-m", "--maxFiles", type=int, default=500, help="Maximum number of files")
parser.add_argument("-w", "--withISR", type=bool, default=True, help="Consider ISR photons")
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

    if args.withISR:
        name = name + "_withISR"

    rdf = ROOT.RDataFrame(chain)
    # rdf = ROOT.RDataFrame('Events', '/afs/cern.ch/work/m/mseidel/WMass/MC/CMSSW_10_6_19_patch2/src/Configuration/WMassNanoGen/NanoGen.root')
    rdf = rdf.Define('genPrefsrLeps', 'Numba::prefsrLeptons(GenPart_status, GenPart_statusFlags, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt)')
    rdf = rdf.Define('ptVgen', 'transversemomentum(GenPart_pt[genPrefsrLeps],GenPart_phi[genPrefsrLeps])')
    rdf = rdf.Define('rapVgen', 'rapidity(GenPart_pt[genPrefsrLeps],GenPart_phi[genPrefsrLeps],GenPart_phi[genPrefsrLeps],GenPart_mass[genPrefsrLeps])')
    rdf = rdf.Define('massVgen', 'invariantmass(GenPart_pt[genPrefsrLeps],GenPart_phi[genPrefsrLeps],GenPart_phi[genPrefsrLeps],GenPart_mass[genPrefsrLeps])')
    rdf = rdf.Define('ewSel', 'Numba::ewPhotonKinematicsSel(GenPart_status, GenPart_statusFlags, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, %s)' % str(args.withISR).lower())
    rdf = rdf.Define('sij', 'log10(invMLepPhotons(GenPart_pt[ewSel==1],GenPart_eta[ewSel==1],GenPart_phi[ewSel==1],0.105658,GenPart_pt[ewSel==3],GenPart_eta[ewSel==3],GenPart_phi[ewSel==3]))')
    rdf = rdf.Define('sik', '(invMLepPhotons(GenPart_pt[ewSel==1|ewSel==2],GenPart_eta[ewSel==1|ewSel==2],GenPart_phi[ewSel==1|ewSel==2],0.105658, GenPart_pt[ewSel==99],GenPart_eta[ewSel==99],GenPart_phi[ewSel==99]))')
    rdf = rdf.Define('sjk', 'log10(invMLepPhotons(GenPart_pt[ewSel==2],GenPart_eta[ewSel==2],GenPart_phi[ewSel==2],0.105658,GenPart_pt[ewSel==3],GenPart_eta[ewSel==3],GenPart_phi[ewSel==3]))')
    rdf = rdf.Define('sijk', '(invMLepPhotons(GenPart_pt[ewSel==1|ewSel==2],GenPart_eta[ewSel==1|ewSel==2],GenPart_phi[ewSel==1|ewSel==2],0.105658, GenPart_pt[ewSel==3],GenPart_eta[ewSel==3],GenPart_phi[ewSel==3]))')
    rdf = rdf.Define('logMassDiff', 'log10(sijk-sik+1e-5)')
    rdf = rdf.Define('pt2ijk', 'pt2ijk(GenPart_pt[ewSel==1],GenPart_eta[ewSel==1],GenPart_phi[ewSel==1],0.105658, GenPart_pt[ewSel==3],GenPart_eta[ewSel==3],GenPart_phi[ewSel==3],0., GenPart_pt[ewSel==2],GenPart_eta[ewSel==2],GenPart_phi[ewSel==2],0.105658)')

    if args.snapshot:
        rdf.Snapshot('Events', 'snapshot_%s.root' % name, ['ptVgen', 'ewSel', 'nGenPart', 'GenPart_pdgId', 'GenPart_status', 'GenPart_statusFlags', 'GenPart_genPartIdxMother', 'GenPart_pt', 'GenPart_eta', 'GenPart_phi', 'sij', 'sik', 'sjk', 'sijk', 'pt2ijk'])

    # hist = rdf.Histo2D((name, ';log_{10} s(#mu- #sum #gamma);log_{10} s(#mu+ #sum #gamma)', 100, -1.5, 3.5, 100, -1.5, 3.5), 'sij', 'sjk')

    if args.boson == 'Z':
        massBins = makeBinning(mass = 91.1535, width = 2.4932, initialStep=0.010)
    else:
        massBins = makeBinning(mass = 80.3815, width = 2.0904, initialStep=0.010)
    massDiffBins = np.linspace(-5, 5, 101)
    hmodel = ROOT.RDF.TH2DModel(name, ';M(#mu+ #mu-) [GeV];log_{10} #{}{ M(#mu+ #mu- #sum #gamma) - M(#mu+ #mu-) + 10^{-5} } [GeV]', len(massBins)-1, array('d', massBins), len(massDiffBins)-1, array('d', massDiffBins))
    hist = rdf.Histo2D(hmodel, 'sik', 'logMassDiff')

    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.12)
    c.SetLeftMargin(0.15)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.SetLogz()
    c.cd()

    hist.Draw('colz')
    c.Update()
    
    st = hist.FindObject('stats')
    st.SetX1NDC(0.60)
    st.SetX2NDC(0.87)
    c.Modified()
    c.Update()

    c.Print('plots/ewPhotonKinematics/%s/ewPK_%s.pdf'  % (args.boson, name))
    c.Print('plots/ewPhotonKinematics/%s/ewPK_%s.png'  % (args.boson, name))
    c.Print('plots/ewPhotonKinematics/%s/ewPK_%s.root' % (args.boson, name))

baseNano = '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoAOD/'
baseNanoGen = '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoGen/'
paths = {}
paths['Z'] = {
    'minnlo': baseNano+'DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV8MC*VFP/210319_*/*/',
    'horace-photos': baseNanoGen+'ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photos-isr-pythia/',
    'horace-pythia': baseNanoGen+'ZToMuMu_TuneCP5_13TeV-horace-born-fsr-pythia-isr-pythia/',
    'horace-exp': baseNanoGen+'ZToMuMu_TuneCP5_13TeV-horace-exp-fsr-off-isr-off/',
    'horace-exp-old': baseNanoGen+'ZToMuMu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia/',
    'horace-alpha-old': baseNanoGen+'ZToMuMu_TuneCP5_13TeV-horace-alpha-old-fsr-off-isr-pythia/',
    'horace-photoslow': baseNanoGen+'ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photoslow-isr-pythia/',
    'horace-photosnopair': baseNanoGen+'ZToMuMu_TuneCP5_13TeV-horace-born-fsr-photosnopair-isr-pythia/',
}
paths['Wplus'] = {
    'minnlo': baseNano+'WplusJetsToMuNu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/MC*VFPWeightFix/210701_234251/*/',
    # 'horace-photosnopair': baseNanoGen+'WplusToMuNu_TuneCP5_13TeV-horace-born-fsr-photosnopair-isr-pythia',
    'horace-photos': baseNanoGen+'WplusToMuNu_TuneCP5_13TeV-horace-born-fsr-photos-isr-pythia',
    'horace-exp': baseNanoGen+'WplusToMuNu_TuneCP5_13TeV-horace-exp-fsr-off-isr-off',
    'horace-exp-old': baseNanoGen+'WplusToMuNu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia',
}
paths['Wminus'] = {
    'minnlo': baseNano+'WminusJetsToMuNu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/MC*VFPWeightFix/210701_234131/*/',
    'horace-photos': baseNanoGen+'WminusToMuNu_TuneCP5_13TeV-horace-born-fsr-photos-isr-pythia',
    'horace-exp': baseNanoGen+'WminusToMuNu_TuneCP5_13TeV-horace-exp-fsr-off-isr-off',
    'horace-exp-old': baseNanoGen+'WminusToMuNu_TuneCP5_13TeV-horace-exp-old-fsr-off-isr-pythia',
}


if args.input == 'all':
    for pathname in paths[args.boson]:
        makePlot(pathname)
else:
    makePlot(args.input)