#!/usr/bin/env python
# may use as: cat wmass/wmass_e/zee_catlist.txt | xargs -i python wmass/wmass_plots.py plots/testZskim {} > runplots.sh
import sys
import re
import os

ODIR=sys.argv[1]

FASTTEST=''
#FASTTEST='--max-entries 1000 '

#dowhat = "plots" 
#dowhat = "dumps" 
dowhat = "yields" 

TREES = "--FMC Friends '{P}/friends/tree_Friend_{cname}.root' "
TREESONLYSKIMW = "-P /data1/emanuele/wmass/TREES_1LEP_53X_V3_WSKIM_V7/" # to be done 
TREESONLYSKIMZ = "-P /data1/emanuele/wmass/TREES_1LEP_53X_V3_ZEESKIM_V7/" # to be done
TREESONLYFULL = "-P /eos/cms/store/group/dpg_ecal/comm_ecal/localreco/TREES_1LEP_80X_V3/"

def base(selection,useSkim=True):

    if 'wenu' in selection: TREESONLYSKIM=TREESONLYSKIMW
    elif 'zee' in selection: TREESONLYSKIM=TREESONLYSKIMZ
    else:
        raise RuntimeError, 'Unknown selection'

    CORE=' '.join([TREES,TREESONLYSKIM])
    if 'pccmsrm29' in os.environ['HOSTNAME']: CORE = CORE.replace('/data1/emanuele/wmass','/u2/emanuele')

    CORE+=" -f -j 4 -l 36.26 --s2v --tree treeProducerWMass --obj tree "+FASTTEST
    if dowhat == "plots": CORE+=" --lspam '#bf{CMS} #it{Preliminary}' --legendWidth 0.20 --legendFontSize 0.035 --showRatio --maxRatioRange 0.75 1.25 --fixRatioRange "

    if selection=='wenu':
        GO="%s wmass/wmass_e/mca-80X-wenu.txt wmass/wmass_e/wenu.txt "%CORE
        GO="%s -W 'puWeight*LepGood_effSF[0]'"%GO
        if dowhat in ["plots","ntuple"]: GO+=" wmass/wmass_e/wenu_plots.txt "
    elif selection=='zee':
        GO="%s wmass/wmass_e/mca-80X-zee.txt wmass/wmass_e/zee.txt "%CORE
        GO="%s -W 'puWeight*LepGood_effSF[0]*LepGood_effSF[1]' --sp 'Z' "%GO
        if dowhat in ["plots","ntuple"]: GO+=" wmass/wmass_e/zee_plots.txt "
    else:
        raise RuntimeError, 'Unknown selection'

    return GO

def procs(GO,mylist):
    return GO+' '+" ".join([ '-p %s'%l for l in mylist ])
def sigprocs(GO,mylist):
    return procs(GO,mylist)+' --showIndivSigs --noStackSig'
def runIt(GO,name,plots=[],noplots=[]):
    if   dowhat == "plots":  print 'python mcPlots.py',"--pdir %s/%s"%(ODIR,name),GO,' '.join(['--sP \'%s\''%p for p in plots]),' '.join(['--xP \'%s\''%p for p in noplots]),' '.join(sys.argv[3:])
    elif dowhat == "yields": print 'echo %s; python mcAnalysis.py'%name,GO,' '.join(sys.argv[3:])
    elif dowhat == "dumps":  print 'echo %s; python mcDump.py'%name,GO,' '.join(sys.argv[3:])
    elif dowhat == "ntuple": print 'echo %s; python mcNtuple.py'%name,GO,' '.join(sys.argv[3:])

def add(GO,opt):
    return '%s %s'%(GO,opt)
def setwide(x):
    x2 = add(x,'--wide')
    x2 = x2.replace('--legendWidth 0.35','--legendWidth 0.20')
    return x2
def fulltrees(x,selection):
    if 'wenu' in selection: TREESONLYSKIM=TREESONLYSKIMW
    elif 'zee' in selection: TREESONLYSKIM=TREESONLYSKIMZ
    else:
        raise RuntimeError, 'Unknown selection'
    return x.replace(TREESONLYSKIM,TREESONLYFULL)

allow_unblinding = True

if __name__ == '__main__':

    torun = sys.argv[2]

    if (not allow_unblinding) and '_data' in torun and (not any([re.match(x.strip()+'$',torun) for x in ['.*_appl.*','cr_.*']])): raise RuntimeError, 'You are trying to unblind!'

    x=""
    if 'zee' in torun:
        x = base('zee')
        if '_ebeb' in torun: x = add(x,"-A alwaystrue ebeb 'max(abs(LepGood1_eta),abs(LepGood2_eta))<1.44'  --scaleSigToData --sP 'z_mll,mZ1' ")
        if '_notebeb' in torun: x = add(x,"-A alwaystrue notebeb 'max(abs(LepGood1_eta),abs(LepGood2_eta))>1.57' --scaleSigToData --sP 'z_mll,mZ1' ")
        if '_gg' in torun: x = add(x,"-A alwaystrue goldgold 'min(LepGood1_r9,LepGood2_r9)>0.94' --scaleSigToData --sP 'z_mll,mZ1' ")
        if '_notgg' in torun: x = add(x,"-A alwaystrue notgoldgold 'min(LepGood1_r9,LepGood2_r9)<0.94' --scaleSigToData --sP 'z_mll,mZ1' ")
        if '_w_reweight' in torun and dowhat=="plots": x = add(x,"--sP 'z_mll,pt1,pt2,ptZ,scaledptZ,costheta_cs,phi_cs,sumAiPi,y_vs_ctheta,y_vs_phi,y_vs_sumAiPi' ")
        if '_genpt' in torun: x = add(x,"--sP 'gen_ptv,gen_scaledptv' --xp 'data' -p 'Z' ")
    elif 'wenu' in torun:
        x = base('wenu')
        if '_w_reweight' in torun: x = x.replace("-W 'puWeight*SF_LepTight_1l'","-W 'puWeight*SF_LepTight_1l*zpt_w*aipi_w'")
        if '_genpt' in torun: x = add(x,"--sP 'gen_ptv,gen_scaledptv' --xp 'data' -p 'W' ")
        
    # skims not ready yet
    x = fulltrees(x,torun)

    
    plots = [] # if empty, to all the ones of the txt file
    if "gen" in torun: noplots = []
    else: noplots = ["^gen_.*"]
    runIt(x,'%s'%torun,plots,noplots)