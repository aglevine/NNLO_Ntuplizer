import FWCore.ParameterSet.Config as cms
import os
import sys
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.parseArguments()

process = cms.Process("PAT")
###### Remove pattuple at output####
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import *

process.TotalEventCounter = cms.EDProducer("EventCountProducer")
process.TotalEventCounter = cms.EDProducer("EventCountProducer")
process.AfterPVFilterCounter = cms.EDProducer("EventCountProducer")
process.AfterNSFilterCounter = cms.EDProducer("EventCountProducer")
process.AfterPATCounter = cms.EDProducer("EventCountProducer")
process.AfterCandidatesCounter = cms.EDProducer("EventCountProducer")
process.AfterJetsCounter = cms.EDProducer("EventCountProducer")



process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
### process.load('Configuration.StandardSequences.Validation_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Generator_cff")
process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.load("DQMServices.Components.DQMFileSaver_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')




#options.inputFiles= 'file:evtfifo.hepmc2g'
#options.outputFile = 'ntuple.root'
#process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(options.maxEvents)
#)
print options.inputFiles
process.source = cms.Source(
    "MCFileSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
)
#process.TFileService = cms.Service(
#    "TFileService",
#    fileName = cms.string('NNLOTest.root')
#)
# Input source
#process.source = cms.Source("MCFileSource",
#    secondaryFileNames = cms.untracked.vstring(),
#    fileNames = cms.untracked.vstring('file:/hdfs/store/user/aglevine/DY_mumu_13TEV_PDF4LHC15nnlo100_Mll_50_150/hepmcDY_mumu_13TEV_PDF4LHC15nnlo100_Mll_50_150_Tune18_MPI_run5006.dat')
#    fileNames = cms.untracked.vstring('file:evtfifo.hepmc2g')
#    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/b/bbilin/BATCH_GEN_SHERPA/Sherpa_4july_z2jetNLO_4jLO/Z4JetNLO.hepmc2g')
#)

#process.TFileService = cms.Service("TFileService",
#                                    fileName = cms.string('ntuple.root' )
#                                    )

process.out = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring("keep *"),
    fileName = cms.untracked.string("test.root")
)



process.tupel = cms.EDAnalyzer("PatBasicAnalyzer",
  #photonSrc = cms.untracked.InputTag("selectedPatPhotonsPFlow"),
  electronSrc = cms.untracked.InputTag("selectedPatElectronsPFlow"),
  #muonSrc = cms.untracked.InputTag("selectedPatMuonsPFlow"),
  #tauSrc = cms.untracked.InputTag("selectedPatTausPFlow"),
  #jetSrc = cms.untracked.InputTag("selectedPatJetsPFlow"),
  #metSrc = cms.untracked.InputTag("selectedPatMetPFlow")
)

process.pgen_novtxsmear = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")+process.GeneInfo+process.genJetMET)

process.genParticles.src = cms.InputTag("source");
process.genParticlesForJets.excludeFromResonancePids = cms.vuint32(12, 14, 16)
process.genParticlesForJets.ignoreParticleIDs = cms.vuint32(1000022, 1000012, 1000014, 1000016, 2000012,
    2000014, 2000016, 1000039, 5100039, 4000012,
    4000014, 4000016, 9900012, 9900014, 9900016,
    39, 12, 14, 16)


process.p = cms.Path(
process.pgen_novtxsmear
*process.tupel
#*process.out
)

#process.e = cms.EndPath(process.out)

process.load("FWCore.MessageService.MessageLogger_cfi")
# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
# process all the events
process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )


#print $inputFileNames
