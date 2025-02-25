#!/bin/env python
from math import *
from PhysicsTools.Heppy.analyzers.core.autovars import NTupleObjectType  
from PhysicsTools.Heppy.analyzers.objects.autophobj import  *
from PhysicsTools.HeppyCore.utils.deltar import deltaR

from CMGTools.TTHAnalysis.signedSip import *
from CMGTools.TTHAnalysis.tools.functionsTTH import _ttH_idEmu_cuts_E2_obj,_soft_MuonId_2016ICHEP,_medium_MuonId_2016ICHEP
from CMGTools.TTHAnalysis.tools.functionsRAX import _susy2lss_idEmu_cuts_obj,_susy2lss_idIsoEmu_cuts_obj

##------------------------------------------  
## LEPTON
##------------------------------------------  

leptonTypeWMass = NTupleObjectType("leptonWMass", baseObjectTypes = [ leptonType ], variables = [
    # Lepton MVA-id related variables
    NTupleVariable("jetDR",      lambda lepton : deltaR(lepton.eta(),lepton.phi(),lepton.jet.eta(),lepton.jet.phi()) if hasattr(lepton,'jet') else -1, help="deltaR(lepton, nearest jet)"),
    NTupleVariable("r9",      lambda lepton : lepton.full5x5_r9() if abs(lepton.pdgId()) == 11 else -99, help="SuperCluster 5x5 r9 variable, only for electrons; -99 for muons"),
    NTupleVariable("chIso04", lambda lepton : lepton.chargedHadronIsoR(0.4), help="Lepton charged hadron isolation in a cone DeltaR=0.4"),
    NTupleVariable("nhIso04", lambda lepton : lepton.neutralHadronIsoR(0.4), help="Lepton neutral hadron isolation in a cone DeltaR=0.4"),
    NTupleVariable("phIso04", lambda lepton : lepton.photonIsoR(0.4), help="Lepton photon isolation in a cone DeltaR=0.4"),
    #2016 muon Id
    NTupleVariable("softMuonId2016", lambda lepton: _soft_MuonId_2016ICHEP(lepton), help="Soft muon ID retuned for ICHEP 2016"),
    NTupleVariable("mediumMuonID2016", lambda lepton: _medium_MuonId_2016ICHEP(lepton), help="Medium muon ID retuned for ICHEP 2016"),
    # More
    NTupleVariable("tightChargeFix",  lambda lepton : ( lepton.isGsfCtfScPixChargeConsistent() + lepton.isGsfScPixChargeConsistent() ) if abs(lepton.pdgId()) == 11 else 2*(lepton.muonBestTrack().ptError()/lepton.muonBestTrack().pt() < 0.2), int, help="Tight charge criteria: for electrons, 2 if isGsfCtfScPixChargeConsistent, 1 if only isGsfScPixChargeConsistent, 0 otherwise; for muons, 2 if ptError/pt < 0.20, 0 otherwise (using the muon best track)"),
    NTupleVariable("muonTrackType",  lambda lepton : 1 if abs(lepton.pdgId()) == 11 else lepton.muonBestTrackType(), int, help="Muon best track type"),
    NTupleVariable("chargeConsistency",  lambda lepton : ( lepton.isGsfCtfScPixChargeConsistent() + lepton.isGsfScPixChargeConsistent() ) if abs(lepton.pdgId()) == 11 else abs(lepton.muonBestTrack().charge() + lepton.innerTrack().charge() + lepton.tunePMuonBestTrack().charge() + ( lepton.globalTrack().charge() + lepton.outerTrack().charge() if lepton.isGlobalMuon() else 0) ), int, help="Tight charge criteria: for electrons, 2 if isGsfCtfScPixChargeConsistent, 1 if only isGsfScPixChargeConsistent, 0 otherwise; for muons, absolute value of the sum of all the charges (5 for global-muons, 3 for global muons)"),
    NTupleVariable("ptErrTk",  lambda lepton : ( lepton.gsfTrack().ptError() ) if abs(lepton.pdgId()) == 11 else (lepton.muonBestTrack().ptError()), help="pt error, for the gsf track or muon best track"),
    # variables for isolated electron trigger matching cuts
    NTupleVariable("ecalPFClusterIso", lambda lepton :  lepton.ecalPFClusterIso() if abs(lepton.pdgId())==11 else -999, help="Electron ecalPFClusterIso"),
    NTupleVariable("hcalPFClusterIso", lambda lepton :  lepton.hcalPFClusterIso() if abs(lepton.pdgId())==11 else -999, help="Electron hcalPFClusterIso"),
    NTupleVariable("dr03TkSumPt", lambda lepton: lepton.dr03TkSumPt() if abs(lepton.pdgId())==11 else -999, help="Electron dr03TkSumPt isolation"),
    NTupleVariable("trackIso", lambda lepton :  lepton.trackIso() if abs(lepton.pdgId())==11 else -999, help="Electron trackIso (in cone of 0.4)"),
    NTupleVariable("trackIso03", lambda lepton :  lepton.isolationR03().sumPt/lepton.pt() if abs(lepton.pdgId())==13 else -999, help="Muon track isolation in cone 0.3"),
    NTupleVariable("etaSc", lambda x : x.superCluster().eta() if abs(x.pdgId())==11 else -100, help="Electron supercluster pseudorapidity"),
    NTupleVariable("energySc", lambda x : x.superCluster().energy() if abs(x.pdgId())==11 else -100, help="Electron supercluster pseudorapidity"),
    NTupleVariable("e5x5", lambda x: x.full5x5_e5x5() if (abs(x.pdgId())==11 and hasattr(x,"full5x5_e5x5")) else -999, help="Electron full5x5_e5x5"),
    NTupleVariable("sigmaIetaIeta", lambda x: x.full5x5_sigmaIetaIeta() if (abs(x.pdgId())==11 and hasattr(x,"full5x5_sigmaIetaIeta")) else -999, help="Electron full5x5_sigmaIetaIeta"),
    NTupleVariable("hcalOverEcal", lambda x: x.full5x5_hcalOverEcal() if (abs(x.pdgId())==11 and hasattr(x,"full5x5_hcalOverEcal")) else -999, help="Electron full5x5_hcalOverEcal"),
    NTupleVariable("eSuperClusterOverP", lambda x: x.eSuperClusterOverP() if (abs(x.pdgId())==11 and hasattr(x,"eSuperClusterOverP")) else -999, help="Electron eSuperClusterOverP"),
    NTupleVariable("matchedTrgObjElePt" , lambda x: x.matchedTrgObjwmassEle.pt() if  x.matchedTrgObjwmassEle else -999.        , help="Matched trigger object (cone dR<0.3) pT to Ele27_WPTight"),
    NTupleVariable("matchedTrgObjEleDR" , lambda x: deltaR(x, x.matchedTrgObjwmassEle) if  x.matchedTrgObjwmassEle else -999.  , help="Matched trigger object (cone dR<0.3) dR to Ele27_WPTight"),
    NTupleVariable("matchedTrgObjMuPt"  , lambda x: x.matchedTrgObjwmassMu.pt() if  x.matchedTrgObjwmassMu else -999.          , help="Matched trigger object (cone dR<0.3) pT to IsoMu24"),
    NTupleVariable("matchedTrgObjMuDR"  , lambda x: deltaR(x, x.matchedTrgObjwmassMu) if  x.matchedTrgObjwmassMu else -999.    , help="Matched trigger object (cone dR<0.3) dR to IsoMu24"),
    NTupleVariable("matchedTrgObjTkMuPt", lambda x: x.matchedTrgObjwmassTkMu.pt() if  x.matchedTrgObjwmassTkMu else -999.      , help="Matched trigger object (cone dR<0.3) pT to IsoTkMu24"),
    NTupleVariable("matchedTrgObjTkMuDR", lambda x: deltaR(x, x.matchedTrgObjwmassTkMu) if  x.matchedTrgObjwmassTkMu else -999., help="Matched trigger object (cone dR<0.3) dR to IsoTkMu24"),
    NTupleVariable("matchedTrgObjMu50Pt", lambda x: x.matchedTrgObjwmassMu50.pt() if  x.matchedTrgObjwmassMu50 else -999.      , help="Matched trigger object (cone dR<0.3) pT to Mu50"),
    NTupleVariable("matchedTrgObjMu50DR", lambda x: deltaR(x, x.matchedTrgObjwmassMu50) if  x.matchedTrgObjwmassMu50 else -999., help="Matched trigger object (cone dR<0.3) dR to Mu50"),
    NTupleVariable("nLayersInner", lambda lepton: lepton.innerTrack().hitPattern().trackerLayersWithMeasurement() if abs(lepton.pdgId()) == 13 else -999., help="Number of layers with measurements in inner track hit pattern for muons."),
])


##------------------------------------------  
## LHE weights with a reduced precision
##------------------------------------------  
lightWeightsInfoType = NTupleObjectType("LightWeightsInfo", mcOnly=True, variables = [
#    NTupleVariable("id",   lambda x : x.id, int),
    NTupleVariable("wgt",   lambda x : x.wgt, storageType="H"),
])

