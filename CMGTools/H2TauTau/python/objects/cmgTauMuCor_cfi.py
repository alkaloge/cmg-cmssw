import FWCore.ParameterSet.Config as cms

cmgTauMuCor = cms.EDProducer(
    "TauMuUpdateProducer",
    # leg1 and leg2 are dummy collections here
    leg1Collection      = cms.InputTag(''),
    leg2Collection      = cms.InputTag(''),
    #metCollection       = cms.InputTag('recoilCorrectedMET'),
    diObjectCollection  = cms.InputTag('cmgTauMuSel'),
    # Read the pruned particles containing the leptons
    genCollection = cms.InputTag('prunedGenParticles'),
    nSigma              = cms.double(0),
    uncertainty         = cms.double(0.03),
    shift1ProngNoPi0    = cms.double(0.),
    # shift1Prong1Pi0     = cms.double(0.015),
    shift1Prong1Pi0     = cms.double(0.012), #Phil's update
    ptDependence1Pi0    = cms.double(0.000),
    # 0.012 according to the TWiki
    shift3Prong         = cms.double(0.012),
    ptDependence3Prong  = cms.double(0.000),
    shiftMet = cms.bool(True),
    shiftTaus = cms.bool(True)
)

# cmgTauMuCor = cms.EDProducer(
#     "TauMuUpdateProducer",
#     cfg  = cmgTauMuCorFactory.clone(),
#     cuts = cms.PSet()
#     )
