import ROOT
import os, array, copy
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True
from math import *
from ctypes import c_float

from CMGTools.WMass.postprocessing.framework.datamodel import Collection 
from CMGTools.WMass.postprocessing.framework.eventloop import Module
from PhysicsTools.HeppyCore.utils.deltar import deltaR

class lepIsoEAProducer(Module):
    def __init__(self,EAfile,rho='rho'):
        self.rho = rho
        self.EAinputfile = EAfile
        if "/EffectiveAreas_cc.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("$CMSSW_RELEASE_BASE/lib/$SCRAM_ARCH/libFWCoreParameterSet.so")
            ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/postprocessing/helpers/EffectiveAreas.cc+" % os.environ['CMSSW_BASE'])

    def beginJob(self):
        self._worker         = ROOT.EffectiveAreas(self.EAinputfile)
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("LepGood_relIso04EA", "F", lenVar="nLepGood")
        self.out.branch("LepGood_customId", "I", lenVar="nLepGood")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def passIDnoIso2016(self,lep,wp='loose'):
        if wp not in ['loose','medium']: 
            print "WARNING! Only loose or medium implemented here"
            return True
        wpthr = {'loose':1, 'medium':2}
        if abs(lep.pdgId)==11:
            if abs(lep.etaSc)<1.479:
                return lep.tightId >= wpthr[wp] and abs(lep.dxy) < 0.05 and abs(lep.dz) < 0.1 and lep.lostHits <= 1 and lep.convVeto == 1
            else:
                return lep.tightId >= wpthr[wp] and abs(lep.dxy) < 0.1 and abs(lep.dz) < 0.2 and lep.lostHits <= 1 and lep.convVeto == 1
        else:
            return True # not yet implemented in friends
    def passCustomIso2016(self,lep,isoEA):
        if abs(lep.pdgId)==11:
            return isoEA < (0.2 if abs(lep.etaSc)<1.479 else 0.0821)
        else: return True
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        leps = Collection(event, "LepGood")
        iso = []
        customId = []
        for l in leps:
            if abs(l.pdgId)!=11: # implemented only for electrons
                iso.append(-1000) 
                customId.append(1)
            else:
                eA = self._worker.getEffectiveArea(abs(l.eta))
                rho = getattr(event,self.rho)
                # approximation, since we don't have the three components of the isolation
                # should be chad + max(0.0, nhad + pho - rho*eA)
                iso.append((l.relIso04*l.pt - rho*eA)/l.pt) 
                # print "eta=%f,pt=%f,eA=%f,rho=%f,reliso=%f,relisocorr=%f" % (l.eta,l.pt,eA,rho,l.relIso04,(l.relIso04*l.pt - rho*eA)/l.pt)
                wp = 'loose' if abs(l.etaSc)<1.479 else 'medium'
                customId.append(self.passCustomIso2016(l,iso[-1]) and self.passIDnoIso2016(l,wp))
        self.out.fillBranch("LepGood_relIso04EA", iso)
        self.out.fillBranch("LepGood_customId", customId)
        return True

class lepAwayJetProducer(Module):
    def __init__(self, lepSel = lambda lep: True, jetSel = lambda jet : True, pairSel = lambda lep, jet : True):
        self.lepSel = lepSel
        self.jetSel = jetSel
        self.pairSel = pairSel
        self.jetSort = lambda jet: jet.pt 
        self.vars = ("pt","eta","phi")
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for V in self.vars:
            self.out.branch("LepGood_awayJet_"+V, "F", lenVar="nLepGood")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        leps = filter(self.lepSel, Collection(event, "LepGood"))
        jets = filter(self.jetSel, Collection(event, "Jet"))
        jets.sort(key = self.jetSort, reverse=True)
        ret = {}
        for V in self.vars: ret[V] = []
        for lep in leps: lep.awayJet = None
        for lep in leps:
            for jet in jets:
                if self.pairSel(lep,jet):
                    lep.awayJet = jet
                    break
            for V in self.vars: 
                ret[V].append(getattr(lep.awayJet,V) if lep.awayJet else 0)
        for V in self.vars:
            self.out.fillBranch("LepGood_awayJet_"+V,ret[V])
        return True

class lepCalibratedEnergyProducer(Module):
    def __init__(self,correctionFile,seed=0,synchronization=False):
        self.corrFile = correctionFile
        self.seed = seed
        self.synchronization = synchronization
        if "/EnergyScaleCorrection_class_cc.so" not in ROOT.gSystem.GetLibraries():
            print 'did not find EnergyScaleCorrection shared object, rebuilding'
            ROOT.gSystem.Load("$CMSSW_RELEASE_BASE/lib/$SCRAM_ARCH/libFWCoreParameterSet.so")
            ROOT.gROOT.ProcessLine(".L %s/src/EgammaAnalysis/ElectronTools/src/EnergyScaleCorrection_class.cc+" % os.environ['CMSSW_BASE'])

        print 'trying to load all the rochester correction libraries'
        if "/RoccoR_cc.so" not in ROOT.gSystem.GetLibraries():
            print 'did not find legacy data roccor. rebuilding...'
            ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/postprocessing/helpers/RoccoR.cc+" % os.environ['CMSSW_BASE'])
            print 'done building legacy data roccor.'

        ##ROOT.gSystem.Load("/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/helpers/RoccoR_cc.so")
        ##ROOT.gSystem.Load("/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/helpers/RoccoR80XMC_cc.so")
        ##ROOT.gSystem.Load("/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/helpers/rochcor201680XMC_cc.so")

    def beginJob(self):
        print 'begin job of lepCalProducer'

        print 'initializing all the rochester correction objects'
        ## rochester corrections for data and mc on legacy rereco
        self._rochester  = ROOT.RoccoR()
        self._rochester.init('/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/data/rochesterCorrections/RoccoR2016.txt')
    
        ## not needed at the moment nLayersHistFile = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/data/rochesterCorrections/nLayersInner_goodmu.root",'read')
        ## not needed at the moment self._nLayersHist = copy.deepcopy(nLayersHistFile.Get("nLayersInner"))

        print 'initializing the egamma calibrator'
        ROOT.gSystem.Load(os.environ['CMSSW_BASE']+"/src/EgammaAnalysis/ElectronTools/src/EnergyScaleCorrection_class_cc.so")
        self._worker = ROOT.EnergyScaleCorrection_class(self.corrFile,self.seed)
        self.rng = ROOT.TRandom3()
        self.rng.SetSeed(self.seed)
        f_resCorr = ROOT.TFile.Open("%s/src/CMGTools/WMass/python/postprocessing/data/leptonScale/el/smoothPtScale_electrons_residualPtCorr.root" % os.environ['CMSSW_BASE'])
        self.h_resCorr = f_resCorr.Get("histSmooth").Clone("dm_diff")
        self.h_resCorr.SetDirectory(None)

    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("LepGood_calPt_step1", "F", lenVar="nLepGood")
        self.out.branch("LepGood_calPt", "F", lenVar="nLepGood")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def gauss(self):
        if self.synchronization: return 1.0
        else: return self.rng.Gaus()
    def residualScale(self,pt,eta,isData):
        if not isData: 
            return 1
        etabin = max(1, min(self.h_resCorr.GetNbinsX(), self.h_resCorr.GetXaxis().FindFixBin(abs(eta))))
        ptbin = max(1, min(self.h_resCorr.GetNbinsY(), self.h_resCorr.GetYaxis().FindFixBin(pt)))
        MZ0 = 91.1876
        scale = 1. - self.h_resCorr.GetBinContent(etabin,ptbin)/MZ0/sqrt(2.)
        if scale<0:
            print "WARNING in residualScale() function: scale < 0 --> returning 0."
            return 0
        else:
            return scale
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        leps = Collection(event, "LepGood")
        calPt_step1 = []
        calPt = []
        for l in leps:
            if abs(l.pdgId)==13: # implemented only for electrons
                if event.isData:
                    scale = self._rochester.kScaleDT(l.charge, l.pt, l.eta, l.phi)
                    calPt_step1.append(scale)
                    calPt.append(l.pt * scale)
                else:
                    ## ATTENTION. THE NTRK IS RANDOMLY CHOSEN FROM A HISTOGRAM DISTRIBUTED LIKE THE SIGNAL MC
                    ## tmp_nlayers = int(self._nLayersHist.GetRandom())
                    scale = self._rochester.kSpreadMC(l.charge, l.pt, l.eta, l.phi, l.mcPt, 0, 0)
                    calPt_step1.append(scale)
                    calPt.append(l.pt * scale )
            else:
                if event.isData:
                    scale = self._worker.ScaleCorrection(event.run,abs(l.etaSc)<1.479,l.r9,abs(l.eta),l.pt)
                    calPt_step1.append(l.pt * scale)
                    calPt.append(l.pt * scale * self.residualScale(l.pt,l.eta,event.isData))
                else:
                    smear = self._worker.getSmearingSigma(event.run,abs(l.etaSc)<1.479,l.r9,abs(l.eta),l.pt,0.,0.)
                    corr = 1.0 + smear * self.gauss()
                    calPt_step1.append(l.pt * corr)
                    calPt.append(l.pt * corr)
        self.out.fillBranch("LepGood_calPt_step1", calPt_step1)
        self.out.fillBranch("LepGood_calPt", calPt)
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

eleRelIsoEA = lambda : lepIsoEAProducer("%s/src/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt" % os.environ['CMSSW_BASE'])
lepQCDAwayJet = lambda : lepAwayJetProducer(jetSel = lambda jet : jet.pt > 30 and abs(jet.eta) < 2.4,
                                            pairSel =lambda lep, jet: deltaR(lep.eta,lep.phi, jet.eta, jet.phi) > 0.7)
eleCalibrated = lambda : lepCalibratedEnergyProducer("CMGTools/WMass/python/postprocessing/data/leptonScale/el/Legacy2016_07Aug2017_FineEtaR9_ele")
