import AthenaPoolCnvSvc.ReadAthenaPool
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
from AthenaCommon.AppMgr import ServiceMgr
from AthenaCommon import CfgMgr

from glob import glob
filelist = ["../../../Run2METReco/WorkArea/run/AOD.pool.root"]
ServiceMgr.EventSelector.InputCollections = filelist

# Set up default configurations
from METUtilities.METUtilConfig import getMETUtilAlg,METUtilConfig

inputdict = {
    }
termdict = {
    }

config = METUtilConfig('RefFinal',termdict,inputdict)
jvftermdict = {
    'Jet':'RefJet_JVFCut',
    }
jvfconfig = METUtilConfig('RefFinal',jvftermdict,inputdict,'RefFinalJVF')
metAlg = getMETUtilAlg('METUtilAlg',[config,jvfconfig])

from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()
topSequence += metAlg

# The tools are accessible via the configurations in metFlags
topSequence.METUtilAlg.OutputLevel = DEBUG
from AthenaCommon.AppMgr import ToolSvc
ToolSvc.MET_Rebuilder_MyRefFinal.OutputLevel = VERBOSE
ToolSvc.MET_Rebuilder_RefFinalJVF.OutputLevel = VERBOSE

write_xAOD = True
if write_xAOD:
    from OutputStreamAthenaPool.MultipleStreamManager import MSMgr
    xaodStream = MSMgr.NewPoolStream( "StreamXAOD", "xAOD.pool.root" )
    xaodStream.AddItem('xAOD::MissingETContainer_v1#MET_RefFinal')
    xaodStream.AddItem('xAOD::MissingETAuxContainer_v1#MET_RefFinalAux.')
    xaodStream.AddItem('xAOD::MissingETContainer_v1#MET_MyRefFinal')
    xaodStream.AddItem('xAOD::MissingETAuxContainer_v1#MET_MyRefFinalAux.')
    xaodStream.AddItem('xAOD::MissingETContainer_v1#MET_RefFinalJVF')
    xaodStream.AddItem('xAOD::MissingETAuxContainer_v1#MET_RefFinalJVFAux.')
    xaodStream.Print()

    ServiceMgr.AthenaPoolCnvSvc.PoolAttributes += [
        "DEFAULT_SPLITLEVEL='99'" ]

    # Force POOL to just simply use the StoreGate keys as branch names:
    ServiceMgr.AthenaPoolCnvSvc.SubLevelBranchName = "<key>"

from Valkyrie.JobOptCfg import ValgrindSvc
svcMgr += ValgrindSvc( OutputLevel = DEBUG,
                       ProfiledAlgs = ["METUtilAlg"],
                       ProfiledIntervals = ["METUtilAlg.execute"])

from PerfMonComps.PerfMonFlags import jobproperties as pmon_properties
pmon_properties.PerfMonFlags.doSemiDetailedMonitoring=True

theApp.EvtMax = -1
ServiceMgr.EventSelector.SkipEvents = 0
