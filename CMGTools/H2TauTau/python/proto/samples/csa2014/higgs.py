import copy
import re 
import PhysicsTools.HeppyCore.framework.config as cfg
from CMGTools.RootTools.yellowreport.YRParser import yrparser13TeV
# from CMGTools.H2TauTau.proto.samples.sampleShift import sampleShift


HiggsVBF125 = cfg.MCComponent(
    name = 'HiggsVBF125',
    files = [],
    xSection = None, 
    nGenEvents = 0,
    triggers = [],
    effCorrFactor = 1 )



############# Gluon fusion ###############


HiggsGGH125 = cfg.MCComponent(
    name = 'HiggsGGH125',
    files = [],
    xSection = None, 
    nGenEvents = 0,
    triggers = [],
    effCorrFactor = 1 )




############# VH  ###############



HiggsVH125 = cfg.MCComponent(
    name = 'HiggsVH125',
    files = [],
    xSection = None, 
    nGenEvents = 0,
    triggers = [],
    effCorrFactor = 1 )







HiggsVBFtoWW125 = cfg.MCComponent(
    name = 'HiggsVBFtoWW125',
    files = [],
    xSection = None, 
    nGenEvents = 0,
    triggers = [],
    effCorrFactor = 1 )

HiggsGGHtoWW125 = cfg.MCComponent(
    name = 'HiggsGGHtoWW125',
    files = [],
    xSection = None, 
    nGenEvents = 0,
    triggers = [],
    effCorrFactor = 1 )

HiggsVHtoWW125 = cfg.MCComponent(
    name = 'HiggsVHtoWW125',
    files = [],
    xSection = None, 
    nGenEvents = 0,
    triggers = [],
    effCorrFactor = 1 )



#############

mc_higgs_vbf = [
    HiggsVBF125,
    HiggsVBFtoWW125,
    ]

mc_higgs_ggh = [
    HiggsGGH125,
    HiggsGGHtoWW125,
    ]

mc_higgs_vh = [
    HiggsVH125,
    HiggsVHtoWW125
    ]


mc_higgs = copy.copy( mc_higgs_vbf )
mc_higgs.extend( mc_higgs_ggh )
mc_higgs.extend( mc_higgs_vh )


pattern = re.compile('Higgs(\D+)(\d+)')
for h in mc_higgs:
    # import pdb ; pdb.set_trace()
    # import pdb; pdb.set_trace()
    m = pattern.match( h.name )
    process = m.group(1)
    if 'toWW' in process :
      process = process.replace('toWW','')
      isToWW = True
    else :
      isToWW = False     
    mass = float(m.group(2))
    xSection = 0.
    try:
        if process=='VH':
            xSection += yrparser13TeV.get(mass)['WH']['sigma']
            xSection += yrparser13TeV.get(mass)['ZH']['sigma']
        else:
            xSection += yrparser13TeV.get(mass)[process]['sigma']
    except KeyError:
        print 'Higgs mass', mass, 'not found in cross section tables. Interpolating linearly at +- 1 GeV...'
        if process=='VH':
            xSection += 0.5 * (yrparser13TeV.get(mass-1.)['WH']['sigma'] + xSection + yrparser13TeV.get(mass+1.)['WH']['sigma'])
            xSection += 0.5 * (yrparser13TeV.get(mass-1.)['ZH']['sigma'] + yrparser13TeV.get(mass+1.)['ZH']['sigma'])
        else:
            xSection += 0.5 * (yrparser13TeV.get(mass-1.)[process]['sigma'] + yrparser13TeV.get(mass+1.)[process]['sigma'])

    if not isToWW :
      br = yrparser13TeV.get(mass)['H2F']['tautau']
    else :
      br = yrparser13TeV.get(mass)['H2B']['WW']
    h.xSection = xSection*br
    h.branchingRatio = br
    print h.name, 'sigma*br =', h.xSection, 'sigma =', xSection, 'br =', h.branchingRatio


# print 'generating shifted samples'
# mc_higgs_up = sampleShift( locals(), mc_higgs, 'Up')
# mc_higgs_down = sampleShift( locals(), mc_higgs, 'Down')
