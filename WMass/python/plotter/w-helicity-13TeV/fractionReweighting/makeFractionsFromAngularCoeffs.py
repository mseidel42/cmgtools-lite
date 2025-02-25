import ROOT, os, copy


## USAGE:
## python makeFractionsFromAngularCoeffs.py -i /afs/cern.ch/user/m/mdunser/www/private/w-helicity-13TeV/angularCoefficients/2018-03-23//moments_muEl/NLO/ -c fractions.root

def formatHisto(hist):
    hist.GetXaxis().SetTitleOffset(1.02)
    hist.GetXaxis().SetTitleSize(0.06)
    hist.GetXaxis().SetLabelSize(0.06)

    hist.GetYaxis().SetTitleOffset(1.02)
    hist.GetYaxis().SetTitleSize(0.06)
    hist.GetYaxis().SetLabelSize(0.06)

    hist.GetZaxis().SetLabelSize(0.06)


def makeFraction(a0, a4, pol, ch):
    rethist = a0.Clone('fraction{p}_{ch}'.format(p=pol,ch=ch))
    rethist.SetTitle('f_{{{p}}} {ch} - constructed '.format(p=pol,ch=ch))
    rethist.SetName ('fraction{p}_{ch}_constructed'.format(p=pol,ch=ch))

    chmult = 1. if ch == 'minus' else 1.  #### WWWWWAAAAHHHHHH

    for ix in range(1,rethist.GetNbinsX()+1):
        for iy in range(1,rethist.GetNbinsY()+1):
            a0val = a0.GetBinContent(ix, iy)
            a4val = a4.GetBinContent(ix, iy)
            newval = 0.
            if pol == 'L':
                newval = 1./4.*(2. - a0val - chmult*a4val)
            if pol == 'R':
                newval = 1./4.*(2. - a0val + chmult*a4val)
            if pol == '0':
                newval = 1./2.*a0val

            rethist.SetBinContent(ix, iy, newval)
    if pol == 'R':
        rethist.GetZaxis().SetRangeUser(0., 0.5)
    if pol == 'L':
        rethist.GetZaxis().SetRangeUser(0., 0.8)
    if pol == '0':
        rethist.GetZaxis().SetRangeUser(0., 0.4)
    
    formatHisto(rethist)

    return rethist


from optparse import OptionParser
parser = OptionParser(usage='%prog [options] cards/card*.txt')
parser.add_option('-v', dest='verbose'    , default=False, action='store_true', help='Save all the bin-by-bin distributions (it\'s a lot of them...')
parser.add_option('-i','--indir', dest='indir', default='', type='string', help='Specify the input directory. all the root files with the angular coeffs should be there.')
parser.add_option('-c','--compare', dest='compfile', default='', type='string', help='Compare the constructed fractions with the fractions from this file.')
(options, args) = parser.parse_args()

if options.compfile:
    compfile = ROOT.TFile(options.compfile, 'read')
    comphists = {}
    for ch in ['plus', 'minus']:
        for pol in ['L', 'R', '0']:
            comphists['fraction{p}_{ch}'.format(p=pol,ch=ch)] = copy.deepcopy(compfile.Get('fraction{p}_{ch}'.format(p=pol,ch=ch)))
    compfile.Close()


ROOT.gROOT.SetBatch()
#ROOT.gStyle.SetOptStat(0)

rootfiles = {}
fractions = {}
rootfile = ''
for inf in os.listdir(options.indir):
    if not '.root' in inf:
        continue
    #if not 'Coeff' in inf:
    #    continue
    #rootfiles[inf.split('Coeff')[0]] = ROOT.TFile(options.indir+'/'+inf, 'read')
    rootfile = ROOT.TFile(options.indir+'/'+inf, 'read')

canv = ROOT.TCanvas()
canv.GetPad(0).SetTopMargin(0.05)
canv.GetPad(0).SetBottomMargin(0.15)
canv.GetPad(0).SetLeftMargin(0.16)
canv.GetPad(0).SetRightMargin(0.15)

for ch in ['plus', 'minus']:
    for pol in ['L', 'R', '0']:
        
        tmp_a0 = copy.deepcopy(rootfile.Get('a0_{ch}'.format(ch=ch)))
        tmp_a4 = copy.deepcopy(rootfile.Get('a4_{ch}'.format(ch=ch)))

        tmp_hist = makeFraction(tmp_a0, tmp_a4, pol, ch)
        tmp_hist.GetXaxis().SetTitle('Y_{W}')
        tmp_hist.GetYaxis().SetTitle('p_{T}^{W}')
        fractions['fraction{p}_{ch}_constr'.format(p=pol,ch=ch)] = tmp_hist
        tmp_hist.Draw('colz')
        canv.SaveAs(options.indir+'/'+tmp_hist.GetName()+'.pdf')
        canv.SaveAs(options.indir+'/'+tmp_hist.GetName()+'.png')
        if options.compfile:
            tmp_compare = copy.deepcopy(comphists['fraction{p}_{ch}'.format(p=pol,ch=ch)])
            tmp_ratio = tmp_hist.Clone('RATIO_fraction{p}_{ch}'.format(p=pol,ch=ch))
            tmp_ratio.Reset()
            tmp_ratio.SetTitle('f_{{{p}}} {ch} ratio: angular / fit'.format(p=pol,ch=ch))
            for ix in range(1,tmp_compare.GetNbinsX()+1):
                for iy in range(1,tmp_compare.GetNbinsY()+1):
                    tmp_ratio.SetBinContent(ix, iy, tmp_hist.GetBinContent(ix, iy) / tmp_compare.GetBinContent(ix, iy) )
            tmp_ratio.GetZaxis().SetRangeUser(0.5, 1.5)
            tmp_ratio.Draw('colz')
            tmp_ratio.GetZaxis().SetRangeUser(0.5, 1.5)
            canv.SaveAs(options.indir+'/'+tmp_ratio.GetName()+'.pdf')
            canv.SaveAs(options.indir+'/'+tmp_ratio.GetName()+'.png')
            fractions['RATIO_fraction{p}_{ch}'.format(p=pol,ch=ch)] = tmp_ratio

    summedFractions = copy.deepcopy(fractions['fractionL_{ch}_constr'.format(ch=ch)]).Clone('fractionSum_{ch}'.format(ch=ch))
    summedFractions.Reset()
    summedFractions.SetName('fractionSum_{ch}'.format(ch=ch))
    summedFractions.SetTitle('#Sigma f_{{i}} {ch}'.format(ch=ch))
    summedFractions.Add(fractions['fractionL_{ch}_constr'.format(ch=ch)])
    summedFractions.Add(fractions['fractionR_{ch}_constr'.format(ch=ch)])
    summedFractions.Add(fractions['fraction0_{ch}_constr'.format(ch=ch)])
    #summedFractions.SetAxisRange(0.9, 1.1, 'Z')
    summedFractions.GetZaxis().SetRangeUser(0.9, 1.1)
    summedFractions.Draw('colz')
    summedFractions.GetZaxis().SetRangeUser(0.9, 1.1)
    canv.Update()
    canv.SaveAs(options.indir+'/'+summedFractions.GetName()+'.pdf')
    canv.SaveAs(options.indir+'/'+summedFractions.GetName()+'.png')
    fractions[summedFractions.GetName()] = copy.deepcopy(summedFractions)


outfile = ROOT.TFile(options.indir+'/'+'fractions_constructed.root', 'recreate')
print 'writing outfile', outfile
for h in fractions.values():
    h.Write()
outfile.Close()
    



