from driverConstants import *
from driverStandard import StandardAnalysis
import driverUtils, sys
options = {
    'ams':OFF,
    'analysisType':STANDARD,
    'applicationName':'analysis',
    'aqua':OFF,
    'ask_delete':OFF,
    'background':None,
    'beamSectGen':OFF,
    'biorid':OFF,
    'cavityTypes':[],
    'cavparallel':OFF,
    'complexFrequency':OFF,
    'contact':OFF,
    'cosimulation':OFF,
    'coupledProcedure':OFF,
    'cpus':1,
    'cse':OFF,
    'cyclicSymmetryModel':OFF,
    'directCyclic':OFF,
    'direct_port':'62997',
    'direct_solver':DMP,
    'dsa':OFF,
    'dynamic':OFF,
    'enrichment':OFF,
    'filPrt':[],
    'fils':[],
    'finitesliding':OFF,
    'foundation':ON,
    'geostatic':OFF,
    'heatTransfer':OFF,
    'importer':OFF,
    'importerParts':OFF,
    'includes':[],
    'initialConditionsFile':OFF,
    'input':'jetnoise_eigen',
    'inputFormat':INP,
    'job':'jetnoise_eigen',
    'keyword_licenses':[],
    'lanczos':ON,
    'libs':[],
    'listener_name':'luigi.cs.wright.edu',
    'listener_resource':'4400',
    'magnetostatic':OFF,
    'massDiffusion':OFF,
    'memory':'90%',
    'message':None,
    'messaging_mechanism':'DIRECT',
    'modifiedTet':OFF,
    'moldflowFiles':[],
    'moldflowMaterial':OFF,
    'mp_mode':THREADS,
    'mp_mode_requested':THREADS,
    'multiphysics':OFF,
    'noDmpDirect':[],
    'noMultiHost':[],
    'noMultiHostElemLoop':['frequency'],
    'no_domain_check':1,
    'outputKeywords':ON,
    'parameterized':OFF,
    'partsAndAssemblies':ON,
    'parval':OFF,
    'postOutput':OFF,
    'preDecomposition':OFF,
    'restart':OFF,
    'restartEndStep':OFF,
    'restartIncrement':0,
    'restartStep':0,
    'restartWrite':ON,
    'rezone':OFF,
    'runCalculator':OFF,
    'soils':OFF,
    'soliter':OFF,
    'solverTypes':['DIRECT'],
    'standard_parallel':ALL,
    'staticNonlinear':OFF,
    'steadyStateTransport':OFF,
    'step':ON,
    'subGen':OFF,
    'subGenLibs':[],
    'subGenTypes':[],
    'submodel':OFF,
    'substrLibDefs':OFF,
    'substructure':OFF,
    'symmetricModelGeneration':OFF,
    'thermal':OFF,
    'tmpdir':'C:\\Users\\WSUadm\\AppData\\Local\\Temp',
    'tracer':OFF,
    'visco':OFF,
}
analysis = StandardAnalysis(options)
status = analysis.run()
sys.exit(status)
