##########################################################
##       CONFIGURATION FOR TTH MULTILEPTON TREES       ##
##########################################################
import PhysicsTools.HeppyCore.framework.config as cfg
import re


#-------- LOAD ALL ANALYZERS -----------

from CMGTools.WMass.analyzers.dmCore_modules_cff import *
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption

#-------- SET OPTIONS AND REDEFINE CONFIGURATIONS -----------

runData = getHeppyOption("runData",False)
runDataQCD = getHeppyOption("runDataQCD",False)
runFRMC = getHeppyOption("runFRMC",False)
scaleProdToLumi = float(getHeppyOption("scaleProdToLumi",-1)) # produce rough equivalent of X /pb for MC datasets
removeJetReCalibration = getHeppyOption("removeJetReCalibration",False)
removeJecUncertainty = getHeppyOption("removeJecUncertainty",False)
doLepCorr = getHeppyOption("doLepCorr",False)
doPhotonCorr = getHeppyOption("doPhotonCorr",False)
doMETpreprocessor = getHeppyOption("doMETpreprocessor",False)
skipT1METCorr = getHeppyOption("skipT1METCorr",False)
forcedSplitFactor = getHeppyOption("splitFactor",-1)
forcedFineSplitFactor = getHeppyOption("fineSplitFactor",-1)
isTest = getHeppyOption("test",None) != None and not re.match("^\d+$",getHeppyOption("test"))
selectedEvents=getHeppyOption("selectEvents","")

# save PDF information and do not skim. Do only for needed MC samples
runOnSignal = True
doTriggerMatching = True
keepLHEweights = True
signalZ = False
diLeptonSkim = False
useBasicRECOLeptons = False
doRecoilVariables = True

# Lepton Skimming
ttHLepSkim.minLeptons = 1
ttHLepSkim.maxLeptons = 999
#ttHLepSkim.idCut  = ""
ttHLepSkim.ptCuts = [23]

if diLeptonSkim:
    ttHLepSkim.minLeptons = 2
    ttHLepSkim.maxLeptons = 999
    #ttHLepSkim.idCut  = ""
    ttHLepSkim.ptCuts = [23, 23]

if doTriggerMatching:
    print 'adding the trigger match analyzer'
    dmCoreSequence.insert(dmCoreSequence.index(triggerFlagsAna)+1, triggerMatchAnaEle )
    dmCoreSequence.insert(dmCoreSequence.index(triggerFlagsAna)+1, triggerMatchAnaMu  )
    dmCoreSequence.insert(dmCoreSequence.index(triggerFlagsAna)+1, triggerMatchAnaTkMu)
    dmCoreSequence.insert(dmCoreSequence.index(triggerFlagsAna)+1, triggerMatchAnaMu50)

# Run miniIso
lepAna.doMiniIsolation = True
lepAna.packedCandidates = 'packedPFCandidates'
lepAna.miniIsolationPUCorr = 'rhoArea'
lepAna.miniIsolationVetoLeptons = None # use 'inclusive' to veto inclusive leptons and their footprint in all isolation cones
lepAna.doIsolationScan = False

# Lepton Preselection
lepAna.loose_electron_id = "POG_Cuts_ID_SPRING16_25ns_v1_HLT"
isolation = None

if useBasicRECOLeptons:
    lepAna.inclusive_muon_id  = None
    lepAna.inclusive_muon_pt  = 23
    lepAna.inclusive_muon_eta = 2.4
    lepAna.inclusive_muon_dxy = 1000.
    lepAna.inclusive_muon_dz  = 1000.
    # loose muon selection
    lepAna.loose_muon_id     = None
    lepAna.loose_muon_pt     = 23
    lepAna.loose_muon_eta    = 2.4
    lepAna.loose_muon_dxy    = 1000.
    lepAna.loose_muon_dz     = 1000.
    lepAna.loose_muon_relIso = 1000.
    # inclusive very loose electron selection
    lepAna.inclusive_electron_id  = None
    lepAna.inclusive_electron_pt  = 25
    lepAna.inclusive_electron_eta = 2.5
    lepAna.inclusive_electron_dxy = 1000.
    lepAna.inclusive_electron_dz  = 1000.
    lepAna.inclusive_electron_lostHits = 100
    # loose electron selection
    lepAna.loose_electron_id     = None #"", #POG_MVA_ID_NonTrig_full5x5",
    lepAna.loose_electron_pt     = 25
    lepAna.loose_electron_eta    = 2.5
    lepAna.loose_electron_dxy    = 1000.
    lepAna.loose_electron_dz     = 1000.
    lepAna.loose_electron_relIso = 1000.
    lepAna.loose_electron_lostHits = 100

    ## ttHLepSkim.minLeptons = 2
    ## ttHLepSkim.maxLeptons = 999
    ## #ttHLepSkim.idCut  = ""
    ## ttHLepSkim.ptCuts = [23, 23]


jetAna.lepSelCut = lambda lep : False # no cleaning of jets with leptons
jetAnaScaleDown.lepSelCut = lambda lep : False # no cleaning of jets with leptons
jetAnaScaleUp.lepSelCut = lambda lep : False # no cleaning of jets with leptons
jetAna.copyJetsByValue = True # do not remove this
metAna.copyMETsByValue = True # do not remove this
jetAna.doQG = True
if not removeJecUncertainty:
    jetAna.addJECShifts = True
    jetAna.jetPtOrUpOrDnSelection = True
    jetAnaScaleDown.copyJetsByValue = True # do not remove this
    jetAnaScaleDown.doQG = False
    metAnaScaleDown.copyMETsByValue = True # do not remove this
    jetAnaScaleUp.copyJetsByValue = True # do not remove this
    jetAnaScaleUp.doQG = False
    metAnaScaleUp.copyMETsByValue = True # do not remove this
    dmCoreSequence.insert(dmCoreSequence.index(jetAna)+1, jetAnaScaleDown)
    dmCoreSequence.insert(dmCoreSequence.index(jetAna)+1, jetAnaScaleUp)
    dmCoreSequence.insert(dmCoreSequence.index(metAna)+1, metAnaScaleDown)
    dmCoreSequence.insert(dmCoreSequence.index(metAna)+1, metAnaScaleUp)


if isolation == "miniIso": 
    lepAna.loose_muon_isoCut     = lambda muon : muon.miniRelIso < 0.4 and muon.sip3D() < 8
    lepAna.loose_electron_isoCut = lambda elec : elec.miniRelIso < 0.4 and elec.sip3D() < 8
elif isolation == None:
    lepAna.loose_muon_isoCut     = lambda muon : True
    lepAna.loose_electron_isoCut = lambda elec : True
elif isolation == "absIso04":
    lepAna.loose_muon_isoCut     = lambda muon : muon.RelIsoMIV04*muon.pt() < 10 and muon.sip3D() < 8
    lepAna.loose_electron_isoCut = lambda elec : elec.RelIsoMIV04*elec.pt() < 10 and elec.sip3D() < 8
else:
    # nothing to do, will use normal relIso03
    pass

# lepton scale / resolution corrections
def doECalElectronCorrections(sync=False,era="25ns"):
    global lepAna, monoJetCtrlLepSkim
    lepAna.doElectronScaleCorrections = {
        'data' : 'EgammaAnalysis/ElectronTools/data/ScalesSmearings/Winter_2016_reReco_v1_ele',
        'GBRForest': ('$CMSSW_BASE/src/CMGTools/RootTools/data/egamma_epComb_GBRForest_76X.root',
                      'gedelectron_p4combination_'+era),
        'isSync': sync
    }
def doECalPhotonCorrections(sync=False):
    global photonAna, gammaJetCtrlSkimmer
    photonAna.doPhotonScaleCorrections = {
        'data' : 'EgammaAnalysis/ElectronTools/data/ScalesSmearings/80X_ichepV2_2016_pho',
        'isSync': sync
    }
def doKalmanMuonCorrections(sync=False,smear="basic"):
    global lepAna
    lepAna.doMuonScaleCorrections = ( 'Kalman', {
        'MC': 'MC_80X_13TeV',
        'Data': 'DATA_80X_13TeV',
        'isSync': sync,
        'smearMode':smear
    })

if doLepCorr: 
    doECalElectronCorrections(era="25ns")
    doKalmanMuonCorrections()
if doPhotonCorr:
    doECalPhotonCorrections()

# Switch off slow photon MC matching
photonAna.do_mc_match = False

# switch off un-needed analyzers
dmCoreSequence.remove(photonAna)

#-------- ADDITIONAL ANALYZERS -----------

## Event Analyzer for susy multi-lepton (at the moment, it's the TTH one)
from CMGTools.TTHAnalysis.analyzers.ttHLepEventAnalyzer import ttHLepEventAnalyzer
ttHEventAna = cfg.Analyzer(
    ttHLepEventAnalyzer, name="ttHLepEventAnalyzer",
    minJets25 = 0,
    )

from CMGTools.WMass.analyzers.treeProducerWMass import * 

# Spring16 electron MVA - follow instructions on pull request for correct area setup
leptonTypeSusy.addVariables([
        NTupleVariable("mvaIdSpring16HZZ",   lambda lepton : lepton.mvaRun2("Spring16HZZ") if abs(lepton.pdgId()) == 11 else 1, help="EGamma POG MVA ID, Spring16, HZZ; 1 for muons"),
        NTupleVariable("mvaIdSpring16GP",   lambda lepton : lepton.mvaRun2("Spring16GP") if abs(lepton.pdgId()) == 11 else 1, help="EGamma POG MVA ID, Spring16, GeneralPurpose; 1 for muons"),
        ])

if not removeJecUncertainty:
    wmass_globalObjects.update({
            "met_jecUp" : NTupleObject("met_jecUp", metType, help="PF E_{T}^{miss}, after type 1 corrections (JEC plus 1sigma)"),
            "met_jecDown" : NTupleObject("met_jecDown", metType, help="PF E_{T}^{miss}, after type 1 corrections (JEC minus 1sigma)"),
            })

if runOnSignal: wmass_globalVariables += pdfsVariables

## Tree Producer
treeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='treeProducerWMass',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
     PDFWeights = PDFWeights,
     globalVariables = wmass_globalVariables,
     globalObjects = wmass_globalObjects,
     collections = wmass_collections,
)


## histo counter
dmCoreSequence.insert(dmCoreSequence.index(skimAnalyzer),
                        histoCounter)

# HBHE new filter
from CMGTools.TTHAnalysis.analyzers.hbheAnalyzer import hbheAnalyzer
hbheAna = cfg.Analyzer(
    hbheAnalyzer, name="hbheAnalyzer", IgnoreTS4TS5ifJetInLowBVRegion=False
    )
dmCoreSequence.insert(dmCoreSequence.index(ttHCoreEventAna),hbheAna)
treeProducer.globalVariables.append(NTupleVariable("hbheFilterNew50ns", lambda ev: ev.hbheFilterNew50ns, int, help="new HBHE filter for 50 ns"))
treeProducer.globalVariables.append(NTupleVariable("hbheFilterNew25ns", lambda ev: ev.hbheFilterNew25ns, int, help="new HBHE filter for 25 ns"))
treeProducer.globalVariables.append(NTupleVariable("hbheFilterIso", lambda ev: ev.hbheFilterIso, int, help="HBHE iso-based noise filter"))
treeProducer.globalVariables.append(NTupleVariable("Flag_badChargedHadronFilter", lambda ev: ev.badChargedHadron, help="bad charged hadron filter decision"))
treeProducer.globalVariables.append(NTupleVariable("Flag_badMuonFilter", lambda ev: ev.badMuon, help="bad muon filter decision"))

#additional MET quantities
metAna.doTkMet = True
treeProducer.globalVariables.append(NTupleVariable("met_trkPt", lambda ev : ev.tkMet.pt() if  hasattr(ev,'tkMet') else  0, help="tkmet p_{T}"))
treeProducer.globalVariables.append(NTupleVariable("met_trkPhi", lambda ev : ev.tkMet.phi() if  hasattr(ev,'tkMet') else  0, help="tkmet phi"))

if not skipT1METCorr:
    if doMETpreprocessor: 
        print "WARNING: you're running the MET preprocessor and also Type1 MET corrections. This is probably not intended."
    jetAna.calculateType1METCorrection = True
    metAna.recalibrate = "type1"
    jetAnaScaleUp.calculateType1METCorrection = True
    metAnaScaleUp.recalibrate = "type1"
    jetAnaScaleDown.calculateType1METCorrection = True
    metAnaScaleDown.recalibrate = "type1"

# add puppi met
puppiMetAna=metAna.clone(
    name='pupimetAnalyzer',
    metCollection='slimmedMETsPuppi',
    noPUMetCollection='slimmedMETsPuppi',
    doTkMet=False,
    includeTkMetCHS=False,
    includeTkMetPVTight=False,
    doMetNoPU=False,
    storePuppiExtra=False,
    collectionPostFix='puppi')

                                         
#-------- SAMPLES AND TRIGGERS -----------


from CMGTools.RootTools.samples.triggers_13TeV_DATA2016 import *
triggerFlagsAna.triggerBits = {
    'DoubleMu' : triggers_mumu_iso,
    'DoubleMuSS' : triggers_mumu_ss,
    'DoubleMuNoIso' : triggers_mumu_noniso,
    'DoubleEl' : triggers_ee,
    'MuEG'     : triggers_mue,
    'SingleMu' : triggers_1mu_iso,
    'SingleMuNoIso' : triggers_1mu_noniso,
    'SingleEl'     : triggers_1e,
}
triggerFlagsAna.unrollbits = True
triggerFlagsAna.saveIsUnprescaled = False
triggerFlagsAna.checkL1Prescale = False

## do some trigger matching for the trigger efficiency studies

from CMGTools.RootTools.samples.samples_13TeV_RunIISummer16MiniAODv2 import *
from CMGTools.RootTools.samples.samples_13TeV_DATA2016 import *
from CMGTools.HToZZ4L.tools.configTools import printSummary, configureSplittingFromTime, cropToLumi, prescaleComponents, insertEventSelector, mergeExtensions
from CMGTools.RootTools.samples.autoAAAconfig import *

selectedComponents = [ DYJetsToLL_M50, WJetsToLNu ]

samples_1fake = [QCD_Mu15] + QCD_Mu5 + QCDPtEMEnriched + QCDPtbcToE + GJetsDR04HT
single_t = [TToLeptons_sch_amcatnlo,T_tch_powheg,TBar_tch_powheg,T_tWch_ext,TBar_tWch_ext] # single top + tW
tt_1l = [TTJets_SingleLeptonFromT,TTJets_SingleLeptonFromT_ext,TTJets_SingleLeptonFromTbar,TTJets_SingleLeptonFromTbar_ext] # TT 1l
w_jets = [WJetsToLNu_LO,WJetsToLNu_LO_ext,WJetsToLNu,WJetsToLNu_ext,WJetsToLNu_ext2v5] # W+jets
z_jets = [DYJetsToLL_M50_LO_ext,DYJetsToLL_M50_LO_ext2,DYJetsToLL_M50] # Z+jets
dibosons = [WW,WW_ext,WZ,WZ_ext,ZZ,ZZ_ext] # di-boson

samples_signal = w_jets
samples_1prompt = single_t + tt_1l + z_jets + dibosons

##configureSplittingFromTime(samples_1fake,30,6)
##configureSplittingFromTime(samples_1prompt,50,6)
##configureSplittingFromTime(samples_signal,100,6)

if   runOnSignal and signalZ==False:
    selectedComponents = [WJetsToLNu_ext2v5]
elif runOnSignal and signalZ:
    selectedComponents = [DYJetsToLL_M50, DYJetsToLL_M50_ext2]

for comp in selectedComponents:
    comp.splitFactor = len(comp.files)/2 #200

if scaleProdToLumi>0: # select only a subset of a sample, corresponding to a given luminosity (assuming ~30k events per MiniAOD file, which is ok for central production)
    target_lumi = scaleProdToLumi # in inverse picobarns
    for c in selectedComponents:
        if not c.isMC: continue
        nfiles = int(min(ceil(target_lumi * c.xSection / 30e3), len(c.files)))
        #if nfiles < 50: nfiles = min(4*nfiles, len(c.files))
        print "For component %s, will want %d/%d files; AAA %s" % (c.name, nfiles, len(c.files), "eoscms" not in c.files[0])
        c.files = c.files[:nfiles]
        c.splitFactor = len(c.files)
        c.fineSplitFactor = 1


if runData != False: # and not isTest: # For running on data

    is50ns = False
    dataChunks = []

    json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt' # 36.5/fb

    run_ranges = []; useAAA=False;
    eras='BCDEFGH'
    if runData in eras:
        eras=runData    
    for era in eras:
        if era=='B':
            processing = "Run2016B-07Aug17_ver2-v1"; short = "Run2016B"; dataChunks.append((json,processing,short,run_ranges,useAAA))
        else:
            processing = "Run2016%s-07Aug17-v1" % era; short = "Run2016%s" % era; dataChunks.append((json,processing,short,run_ranges,useAAA))
    
    DatasetsAndTriggers = []
    selectedComponents = [];
    exclusiveDatasets = False; # this will veto triggers from previous PDs in each PD, so that there are no duplicate events
 
    # DatasetsAndTriggers.append( ("DoubleMuon", triggers_mumu_iso + triggers_mumu_ss + triggers_mumu_ht + triggers_3mu + triggers_3mu_alt) )
    # DatasetsAndTriggers.append( ("DoubleEG",   triggers_ee + triggers_ee_ht + triggers_3e) )
    # DatasetsAndTriggers.append( ("MuonEG",     triggers_mue + triggers_mue_ht + triggers_2mu1e + triggers_2e1mu) )
    DatasetsAndTriggers.append( ("SingleMuon", triggers_1mu_iso ) ) # + triggers_1mu_noniso) )
    ##DatasetsAndTriggers.append( ("SingleElectron", triggers_1e) )

    if runDataQCD: # for fake rate measurements in data
        FRTrigs_mu = triggers_FR_1mu_noiso
        FRTrigs_el = triggers_FR_1e_noiso
        DatasetsAndTriggers = [
            ("DoubleMuon", FRTrigs_mu ),
            ("DoubleEG",   FRTrigs_el ),
            #("JetHT",   triggers_FR_jet )
        ]
        triggers_FR_muNoIso = [ 'HLT_Mu27_v*', 'HLT_Mu50_v*' ]
        triggerAna.myTriggerPrescales = { 'HLT_Mu50_v*':10 }
        exclusiveDatasets = False

    for json,processing,short,run_ranges,useAAA in dataChunks:
        if len(run_ranges)==0: run_ranges=[None]
        vetos = []
        for pd,triggers in DatasetsAndTriggers:
            for run_range in run_ranges:
                label = ""
                if run_range!=None:
                    label = "_runs_%d_%d" % run_range if run_range[0] != run_range[1] else "run_%d" % (run_range[0],)
                compname = pd+"_"+short+label
                comp = kreator.makeDataComponent(compname, 
                                                 "/"+pd+"/"+processing+"/MINIAOD",
                                                 "CMS", ".*root", 
                                                 json=json, 
                                                 run_range=(run_range if "PromptReco" not in processing else None), 
                                                 triggers=triggers[:], vetoTriggers = vetos[:],
                                                 useAAA=useAAA)
                if "PromptReco" in processing:
                    from CMGTools.Production.promptRecoRunRangeFilter import filterComponent
                    filterComponent(comp, verbose=1)
                print "Will process %s (%d files)" % (comp.name, len(comp.files))
                #comp.splitFactor = len(comp.files)/3
                comp.splitFactor = len(comp.files)/2
                comp.fineSplitFactor = 1
                selectedComponents.append( comp )
            if exclusiveDatasets: vetos += triggers
    if json is None:
        dmCoreSequence.remove(jsonAna)
    if runDataQCD: # for fake rate measurements in data
         configureSplittingFromTime(selectedComponents, 3.5, 2, maxFiles=15)

if True:
    from CMGTools.Production.promptRecoRunRangeFilter import filterComponent
    for c in selectedComponents:
        printnewsummary = False
        if "PromptReco" in c.name:
            printnewsummary = True
            filterComponent(c, 1)
            c.splitFactor = len(c.files)/6
    if printnewsummary: printSummary(selectedComponents)


if removeJetReCalibration:
    jetAna.recalibrateJets = False
    jetAnaScaleUp.recalibrateJets = False
    jetAnaScaleDown.recalibrateJets = False

if getHeppyOption("noLepSkim",False):
    ttHLepSkim.minLeptons=0

if forcedSplitFactor>0 or forcedFineSplitFactor>0:
    if forcedFineSplitFactor>0 and forcedSplitFactor!=1: raise RuntimeError, 'splitFactor must be 1 if setting fineSplitFactor'
    for c in selectedComponents:
        if forcedSplitFactor>0: c.splitFactor = forcedSplitFactor
        if forcedFineSplitFactor>0: c.fineSplitFactor = forcedFineSplitFactor

if selectedEvents!="":
    events=[ int(evt) for evt in selectedEvents.split(",") ]
    print "selecting only the following events : ", events
    eventSelector= cfg.Analyzer(
        EventSelector,'EventSelector',
        toSelect = events
        )
    dmCoreSequence.insert(dmCoreSequence.index(lheWeightAna), eventSelector)

#-------- SEQUENCE -----------

sequence = cfg.Sequence(dmCoreSequence+[
        puppiMetAna,
        ttHEventAna,        
        treeProducer,
    ])
preprocessor = None


if doRecoilVariables:
    print 'adding recoil variables'
    from CMGTools.WMass.analyzers.eventRecoilProducer import *
    eventRecoilProducer = cfg.Analyzer(EventRecoilProducer,
                                       name='eventRecoilProducer',
                                       pf='packedPFCandidates',
                                       pfType='std::vector<pat::PackedCandidate>',
                                       gen='packedGenParticles',
                                       genType='std::vector<pat::PackedGenParticle>')
    sequence.insert(sequence.index(ttHEventAna)+1,eventRecoilProducer)
    treeProducer.globalVariables += wmass_recoilVariables



#-------- HOW TO RUN -----------

print '====================================================================='
print '====================================================================='
print '====================================================================='
print '====================================================================='
print '====================================================================='
print runData
print 'THIS IS SELECTED COMPONENTS', selectedComponents

test = getHeppyOption('test')
if test in[ 'testw' , 'testz' , 'testdata' , 'testwnew' , 'testznew']:
    if test=='testdata':
        comp = selectedComponents[0]
        comp.files = ['/eos/cms/store/data/Run2016C/SingleElectron/MINIAOD/03Feb2017-v1/50000/AEA181FD-61EB-E611-B54D-1CC1DE18CFF6.root']
    if test=='testw':
        comp = WJetsToLNu_ext
        #comp.files = ['/eos/user/m/mdunser/w-helicity-13TeV/localFilesMINIAOD/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root']
        comp.files = comp.files[:1]
    elif test=='testw_nnpdf31':
        comp = WJetsToLNu
        # comp.files = ['/eos/cms/store/mc/RunIISummer16MiniAODv2/WplusJToENuJ_scalesUpTo8_NNPDF31_plus_CMSPDF_TuneCP5_13TeV_powheg2-minlo-pythia/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/30000/34F84C69-917E-E811-BB1F-FA163EB2ABE6.root'] # we+
        # comp.files = ['/eos/cms/store/mc/RunIISummer16MiniAODv2/WminusJToENuJ_scalesUpTo8_NNPDF31_plus_CMSPDF_TuneCP5_13TeV_powheg2-minlo-pythia/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/30000/2C26E10C-CE7F-E811-924B-FA163EE6807E.root'] # we-
        comp.files = ['/eos/cms/store/mc/RunIISummer16MiniAODv2/WplusJToMuNuJ_scalesUpTo8_NNPDF31_plus_CMSPDF_TuneCP5_13TeV_powheg2-minlo-pythia/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v4/10000/04783B3A-B57D-E811-B421-FA163EF059FC.root'] # wm+
    elif test=='testz_nnpdf31':
        comp = DYJetsToLL_M50
        comp.files = ['/eos/cms/store/mc/RunIISummer16MiniAODv2/ZJToEEJ_M-50_scalesUpTo8_NNPDF31_plus_CMSPDF_TuneCP5_13TeV_powheg2-minlo-pythia/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/70000/F22FBFAF-0179-E811-9038-FA163E54DF16.root']
    elif test=='testz':
        comp=DYJetsToLL_M50
        comp.files=comp.files[:1]#'/eos/cms/store/cmst3/user/psilva/Wmass/DYJetsMG_test/FCDD4D28-12C4-E611-8BFC-C4346BC8F6D0.root']
    elif test=='testznew':
        comp=DYJetsToLL_M50
        comp.files=['/eos/user/m/mdunser/w-helicity-13TeV/testfiles/ZJToMuMuJ_M-50_scalesUpTo8_NNPDF31_plus_CMSPDF_TuneCP5_13TeV_powheg2-minlo-pythia.root']
    elif test=='testwnew':
        comp=WJetsToLNu_LO
        comp.files=['/eos/cms/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/06C6936D-7EBC-E611-B990-0025905A60B4.root']
    else:
        #comp=MuonEG_Run2016B_18Apr2017
        #comp=SingleElectron_Run2016C_03Feb2017
        comp = selectedComponents[0]
        comp.files = comp.files[:1]
        #comp.files=['/eos/cms/store/data/Run2016G/SingleMuon/MINIAOD/18Apr2017-v1/120000/1C169863-7541-E711-81BD-1CC1DE1CEFE0.root']
        #comp.files=['/eos/user/m/mdunser/w-helicity-13TeV/testfiles/Run2016G-SingleMuon-MINIAOD-18Apr2017-v1-120000_1C169863-7541-E711-81BD-1CC1DE1CEFE0.root']
        #comp.files=['/eos/cms/store/data/Run2016D/MuonEG/MINIAOD/03Feb2017-v1/80000/92E2F2D7-BAEA-E611-A725-0090FAA583C4.root']
        #comp.files=['/eos/cms/store/data/Run2016C/SingleElectron/MINIAOD/03Feb2017-v1/50000/AEA181FD-61EB-E611-B54D-1CC1DE18CFF6.root']
    print comp.files
    comp.splitFactor = 1
    comp.fineSplitFactor = 1
    selectedComponents = [ comp ]
elif test == '2':
    sel = getHeppyOption('sel','.*')
    for comp in selectedComponents[:]:
        if sel and not any(re.search(p.strip(),comp.name) for p in sel.split(",")):
            selectedComponents.remove(comp)
            continue
        comp.files = comp.files[:1]
        comp.splitFactor = 1
        comp.fineSplitFactor = 1
elif test == '3':
    selectedComponents = selectedComponents[:1]
    for comp in selectedComponents:
        comp.files = comp.files[:1]
        comp.splitFactor = 1
        comp.fineSplitFactor = 4
elif test == '5':
    for comp in selectedComponents:
        comp.files = comp.files[:5]
        comp.splitFactor = 1
        comp.fineSplitFactor = 5
elif test == '21':
    for comp in selectedComponents:
        comp.files = comp.files[:7]
        comp.splitFactor = 1
        comp.fineSplitFactor = 3
elif test == '80X-MC':
    what = getHeppyOption("sample","DYLL")
    if what == "DYLL":
        DYJetsToLL_M50 = kreator.makeMCComponent("DYJetsToLL_M50", "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM", "CMS", ".*root", 2008.*3)
        selectedComponents = [ DYJetsToLL_M50 ]
        comp = selectedComponents[0]
        comp.triggers = []
        comp.files = [ '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/101D622A-85C4-E611-A7C2-C4346BC80410.root' ]
        tmpfil = os.path.expandvars("/tmp/$USER/101D622A-85C4-E611-A7C2-C4346BC80410.root")
        if not os.path.exists(tmpfil):
            os.system("xrdcp root://eoscms//eos/cms%s %s" % (comp.files[0],tmpfil))
        comp.files = [ tmpfil ]
        comp.splitFactor = 1
        if not getHeppyOption("single"): comp.fineSplitFactor = 4
    elif what == "WJets" or what == "WJetsNNPDF31":
        if what == "WJetsNNPDF31":
            WJetsToLNu = kreator.makeMCComponent("WJetsToLNu_NNPDF31", "/WplusJToENuJ_scalesUpTo8_NNPDF31_plus_CMSPDF_TuneCP5_13TeV_powheg2-minlo-pythia/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM", "CMS", ".*root", 3* 20508.9)            
            comp.files = [ '/store/mc/RunIISummer16MiniAODv2/WplusJToENuJ_scalesUpTo8_NNPDF31_plus_CMSPDF_TuneCP5_13TeV_powheg2-minlo-pythia/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/30000/B64A576F-7466-E811-8744-0025905A60B6.root' ]
            tmpfil = os.path.expandvars("/tmp/$USER/B64A576F-7466-E811-8744-0025905A60B6.root")
        else:
            WJetsToLNu = kreator.makeMCComponent("WJetsToLNu", "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", "CMS", ".*root", 3* 20508.9)
            comp.files = [ '/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/160DFF65-E8BF-E611-A4AD-0CC47AC08C1A.root' ]
            tmpfil = os.path.expandvars("/tmp/$USER/160DFF65-E8BF-E611-A4AD-0CC47AC08C1A.root")
        selectedComponents = [ WJetsToLNu ]
        comp = selectedComponents[0]
        comp.triggers = []
        if not os.path.exists(tmpfil):
            os.system("xrdcp root://eoscms//eos/cms%s %s" % (comp.files[0],tmpfil))
        comp.files = [ tmpfil ]
        comp.splitFactor = 1
        if not getHeppyOption("single"): comp.fineSplitFactor = 4
    else: raise RuntimeError, "Unknown MC sample: %s" % what
elif test == '80X-Data':
    DoubleMuon = kreator.makeDataComponent("DoubleMuon_Run2016B_run274315", "/DoubleMuon/Run2016B-18Apr2017_ver2-v1/MINIAOD", "CMS", ".*root", run_range = (274315,274315), triggers = triggers_mumu)
    DoubleEG = kreator.makeDataComponent("DoubleEG_Run2016B_run274315", "/DoubleEG/Run2016B-18Apr2017_ver2-v1/MINIAOD", "CMS", ".*root", run_range = (274315,274315), triggers = triggers_ee)
    DoubleMuon.files = [ 'root://xrootd-cms.infn.it//store/data/Run2016B/DoubleMuon/MINIAOD/18Apr2017_ver2-v1/100000/DA4A8A8A-7436-E711-BDB1-24BE05C6D731.root' ]
    DoubleEG.files = [ 'root://xrootd-cms.infn.it//store/data/Run2016B/DoubleEG/MINIAOD/18Apr2017_ver2-v1/00000/689FEB7E-DE3E-E711-8815-001E67A41EA0.root' ]
    selectedComponents = [ DoubleMuon, DoubleEG ]
    for comp in selectedComponents:
        comp.json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
        tmpfil = os.path.expandvars("/tmp/$USER/%s" % os.path.basename(comp.files[0]))
        if not os.path.exists(tmpfil): os.system("xrdcp %s %s" % (comp.files[0],tmpfil)) 
        comp.files = [tmpfil]
        comp.splitFactor = 1
        comp.fineSplitFactor = 1
elif test != None:
    raise RuntimeError, "Unknown test %r" % test

## FAST mode: pre-skim using reco leptons, don't do accounting of LHE weights (slow)"
## Useful for large background samples with low skim efficiency
if getHeppyOption("fast"):
    histoCounter.doLHE = False
    from CMGTools.TTHAnalysis.analyzers.ttHFastLepSkimmer import ttHFastLepSkimmer
    fastSkim = cfg.Analyzer(
        ttHFastLepSkimmer, name="ttHFastLepSkimmer1lep",
        muons = 'slimmedMuons', muCut = lambda mu : mu.pt() > 10 and mu.isLooseMuon(),
        electrons = 'slimmedElectrons', eleCut = lambda ele : ele.pt() > 15,
        minLeptons = 1, 
    )
    if jsonAna in sequence:
        sequence.insert(sequence.index(jsonAna)+1, fastSkim)
    else:
        sequence.insert(sequence.index(skimAnalyzer)+1, fastSkim)

if not keepLHEweights:
    if "LHE_weights" in treeProducer.collections: treeProducer.collections.pop("LHE_weights")
    if lheWeightAna in sequence: sequence.remove(lheWeightAna)
    histoCounter.doLHE = False

if runOnSignal and not signalZ:
    if ttHLepSkim in sequence: sequence.remove(ttHLepSkim)
    if triggerAna in sequence: sequence.remove(triggerAna)
    if genAna in sequence: genAna.saveAllInterestingGenParticles = True

## Auto-AAA
if not getHeppyOption("isCrab"):
    autoAAA(selectedComponents)

## output histogram
outputService=[]
from PhysicsTools.HeppyCore.framework.services.tfile import TFileService
output_service = cfg.Service(
    TFileService,
    'outputfile',
    name="outputfile",
    fname='treeProducerWMass/tree.root',
    option='recreate'
    )    
outputService.append(output_service)

# print summary of components to process
printSummary(selectedComponents)

# the following is declared in case this cfg is used in input to the heppy.py script
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
from CMGTools.TTHAnalysis.tools.EOSEventsWithDownload import EOSEventsWithDownload
event_class = EOSEventsWithDownload if not preprocessor else Events
EOSEventsWithDownload.aggressive = 2 # always fetch if running on Wigner
if getHeppyOption("nofetch") or getHeppyOption("isCrab"):
    event_class = Events
    if preprocessor: preprocessor.prefetch = False
config = cfg.Config( components = selectedComponents,
                     sequence = sequence,
                     services = outputService, 
                     preprocessor = preprocessor, 
                     events_class = event_class)
