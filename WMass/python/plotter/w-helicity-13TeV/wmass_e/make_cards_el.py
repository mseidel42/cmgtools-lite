import os
from datetime import datetime

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
parser.add_option("-d", "--dry-run", dest="dryRun",   action="store_true", default=False, help="Do not run the job, only print the command");
parser.add_option("-s", "--suffix", dest="suffix", type="string", default=None, help="Append a suffix to the default outputdir (helicity_<date>)");
parser.add_option("-q", "--queue", dest="queue", type="string", default="2nd", help="Select the queue to use");
parser.add_option(      "--syst", dest="addSyst", action="store_true", default=False, help="Add PDF and QCD scale systematics to the signal (need incl_sig directive in the MCA file)");
(options, args) = parser.parse_args()

BASECONFIG="w-helicity-13TeV/wmass_e"
PROG="w-helicity-13TeV/make_helicity_cards.py" 
MCA=BASECONFIG+'/mca-80X-wenu-helicity.txt'    
CUTFILE=BASECONFIG+'/wenu_80X.txt'
SYSTFILE=BASECONFIG+'/systsEnv.txt'
TREEPATH="/afs/cern.ch/work/e/emanuele/TREES/TREES_electrons_1l_V6_TINY"
QUEUE=str(options.queue)
VAR="\"ptElFull(LepGood1_calPt,LepGood1_eta):LepGood1_eta\""
BINNING="\"[-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.566,-1.4442,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.566,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5]*[30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45]\""
WEIGHTSTRING="\'puw2016_nTrueInt_36fb(nTrueInt)*lepSF(LepGood1_pdgId,LepGood1_pt,LepGood1_eta,LepGood1_SF1,LepGood1_SF2,LepGood1_SF3)\' "
LUMI=35.9

OUTDIR="helicity_%s" % datetime.now().strftime("%Y_%m_%d")
if options.suffix: OUTDIR += ("_%s" % options.suffix)


components=[" -s "," -b "]

for c in components:
    cmd="python " + " ".join([PROG,MCA,CUTFILE,VAR,BINNING,SYSTFILE,OUTDIR]) + \
        (" -W %s " % WEIGHTSTRING) + (" -P %s " % TREEPATH) + (" -q %s " % QUEUE) + c + \
        (" -l %f " % LUMI)
    if options.dryRun: cmd += '  --dry-run '
    if options.addSyst: cmd += '  --pdf-syst --qcd-syst '
    os.system(cmd)
