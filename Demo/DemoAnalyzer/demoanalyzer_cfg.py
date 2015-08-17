import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
#import PhysicsTools.PythonAnalysis.LumiList as LumiList

process = cms.Process("Demo")

# intialize MessageLogger and output report

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.FwkReport.reportEvery = 100


#process.load("FWCore.MessageService.MessageLogger_cfi")

#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
       # limit = cms.untracked.int32(-1)
       # )
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
                                        
# set the maximum number of events to be processed here
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

goodJSON = '/home/cms-opendata/CMSSW_4_2_8/src/Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_v2.txt'
#goodJSON = '/afs/desy.de/user/s/stefan/CMSPhysic/FromIrene/data_files/Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_v2.txt'

myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')


files2010data = FileUtils.loadListFromFile ('/home/cms-opendata/CMSSW_4_2_8/src/MinBsData0000.txt')
#files2010data = FileUtils.loadListFromFile ('/home/cms-opendata/CMSSW_4_2_8/src/JetData0000.txt')
#files2010data = FileUtils.loadListFromFile ('/home/cms-opendata/CMSSW_4_2_8/src/MuoniaData0000.txt')


#files2010data = FileUtils.loadListFromFile ('/afs/desy.de/user/s/stefan/CMSPhysic/FromIrene/data_files/CMS_Run2010B_Mu_AOD_Apr21ReReco-v1_0000_file_index.txt')

#files2010data = FileUtils.loadListFromFile ('/home/cms-opendata/CMSSW_4_2_8/src/CMS_Run2010B_Mu_AOD_Apr21ReReco-v1_0000_file_index.txt')

#files2010data = FileUtils.loadListFromFile ('/home/cms-opendata/CMSSW_4_2_8/src/CMS_Run2010B_MinimumBias_AOD_Apr21ReReco-v1_0000_file_index.txt')

#files2010data = FileUtils.loadListFromFile ('/home/cms-opendata/CMSSW_4_2_8/src/CMS_Run2010B_ZeroBias_AOD_Apr21ReReco-v1_0008_file_index.txt')

#files2010data = FileUtils.loadListFromFile ('/home/cms-opendata/CMSSW_4_2_8/src/CMS_Run2010B_MuOnia_AOD_Apr21ReReco-v1_0000_file_index.txt')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(*files2010data
    
    )
)
'''
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring('root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Mu/AOD/Apr21ReReco-v1/0000/00459D48-EB70-E011-AF09-90E6BA19A252.root'
    
    )
)
'''

    
#process.source.lumisToProcess=cms.untracked.VLuminosityBlockRange()
#process.source.lumisToProcess.extend(myLumis)

# this would be the "unofficial" way, which works and is an alternative
#please un-comment the import CfgTypes and PhysicsTools.PythonAnalysis.LumiList
#to use it.Also comment out the FWCore.PythonUtilities.Lumilist
#   (needs to be placed *after* the process.source input file definition!)


process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)





process.demo = cms.EDAnalyzer('DemoAnalyzer'
)

# change this output file name according to input file
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('OutputD02.root')
                                   )

process.p = cms.Path(process.demo)
