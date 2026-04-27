### MAIN FUNCTIONS CARRING OUT THE JOBS
import correlators_script as cs
import effective_masses_script as efs
import eigenvalues_script as evs
import fitting_script as fts
import add_rem_operators as rcs

### MAIN ANALYSIS FUNCTIONS
import bin_size_script as bs
import dispersion_relation_script as drs
import correlated_differences_script as cds

### LAYOUT AND EXTRA FUNCTIONS
import set_of_layout_functions as vfl
import set_of_analysis_functions as vfa

### EXTRA LIBRARIES
import ast
import numpy as np
import os
import h5py
import sys


import ensembles as ed

### ------ EXTRA VARIABLES -----
# runSorting = False
# runReduceDataSet = False

### TO CHECK 
# myNrIrreps = []
# myNrOps = []
# runPlotBinOnly = False


### ------- WHAT IT IS DONE --------
myArgs = vfl.parse_args()
myRuns = vfl.WhichRuns(myArgs, ed.ensembles)

runCorrs = myRuns.corrs
runEigenvals = myRuns.eigenvals
runRowsCols = myRuns.ops
runEffMass = myRuns.effmass
runFits = myRuns.fits
runDispRel = myRuns.disp
myTypeRs = myRuns.rs_type
myRebinOn = myRuns.rebin
myRb = myRuns.rb
myKbt = myRuns.kbt

### This is the amount of irreps to compute or when to start and when to finish the analysis
myNrIrreps = None # 2 # 1
myFirstIrrep = None # None # 1 # 2
myLastIrrep = None # None

### Other parameters
myKbtSamples = None #np.array(np.loadtxt('bootstrap_samples.txt')) #None
myEffMassDistance = myRuns.dist_eff_mass #None #2 #3


### ------ MAIN VARIABLES ---------
myEns = myRuns.ensemble
myWhichCorrelator = myRuns.correlator

### ------ CHECKING EVERYTHING ---------

vfl.VALIDATE_RUNS(myRuns)

### This is just naming schemes
reBin = f"_bin{myRb}" if myRuns.rebin else ""

### Information from the ensembles dictionary
myLocation = vfl.DIRECTORY_EXISTS(f'{ed.outputLocation}{myEns}/')


### ----  GETTING THE INFO FROM THE FOLLOWING FILES -------

myCnfgs = ed.ensembles[myEns]['ncfgs'] # None # 20 # 100
myWeight = vfa.REWEIGHTS(ed.ensembles[myEns]['weight_raw'], myCnfgs)

if not ed.ensembles[myEns]['allConfigs']:
    myTempCnfgs = ed.ensembles[myEns]['nfgsList']
    myCnfgs = len(myTempCnfgs)
    myWeight = np.asarray(vfa.RW_NORMALIZATION([myWeight[ii] for ii in myTempCnfgs], myCnfgs), dtype=np.float128)

myLatSize = ed.ensembles[myEns]['LatSize']


### -------- PRINTING INFO OF ENSEMBLE ---------

vfl.INFO_PRINTING(myWhichCorrelator, myEns)



### ------------ START --------------

##  Single Hadron correlators
if myWhichCorrelator =='s':
    myVersion =  f'_{myEns}_singles_test' 
    
    myArchivo = h5py.File(ed.ensembles[myEns]['fs'], 'r')
    myIrreps = list(myArchivo.keys())        
    
    ### Correlators analysis
    if runCorrs: 
        locationWorkedCorrelators = cs.SingleCorrelatorAnalysis(myArchivo, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps, own_kbt_list = myKbtSamples, first_irrep = myFirstIrrep , last_irrep = myLastIrrep)
    else:
        locationWorkedCorrelators = f'{myLocation}Single_correlators_{myTypeRs}{reBin}_v{myVersion}.h5'
    
    try:
        myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    except FileNotFoundError:
        print('Cannot find the correlator file for further analysis.')
        sys.exit()
        
    ### Effective Masses analysis
    if runEffMass: 
        efs.SingleCorrelatorEffectiveMass(myCorrelator, myTypeRs, dist_eff_mass=myEffMassDistance) 
        
    ### Fits analysis
    if runFits: 
        myOneTMin = myRuns.fit.one_tmin
        myTypeFit = myRuns.fit.type_fit
        myTypeCorrelation = myRuns.fit.type_correlation
        myTMaxList = ed.ensembles[myEns]['singleTMaxFits']
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myLocation}Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(f'{myFitsLocation}Single_correlators_{myTypeRs}{reBin}_fits_v{myVersion}.h5', 'a')
        
        fts.FitSingleCorrelators(myCorrelator, myFitCorrelator, myTypeRs, myTMaxList, myIrreps, one_tmin = myOneTMin, type_fit = myTypeFit, type_correlation = myTypeCorrelation)
        
        myFitCorrelator.close()
    
    ### Dispersion relation analysis
    if runDispRel:
        myDispMode = myRuns.disp_run.mode
        myTMinSPlot = ed.ensembles[myEns]['singleTMinResults']
        myDispLocation = vfl.DIRECTORY_EXISTS(f'{ed.location}/Plots/{myEns}/Dispersion_Relation/')
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myLocation}Fits_SingleHadrons/')
        myDispFile  = h5py.File(f'{myFitsLocation}Single_correlators_{myTypeRs}{reBin}_fits_v{myVersion}.h5', 'a')
        
        drs.DISPERSION_RELATION(myDispFile, myTypeRs, myTMinSPlot, myLatSize, myDispLocation, myDispMode, myEns)
        
        myDispFile.close()
        
        ### THIS NEEDS TO BE CHECKED
#     ### Binning analysis
#     # if any([runCorrs, runBinSizeAnalysis]):
#         # vfl.SINGLE_INFO_PRINTING(myIrreps, runBinSizeAnalysis, runCorrs)
#     if runBinSizeAnalysis:
#         myBinSizeIrrep = ast.literal_eval(input('Enter irrep for bin-size analysis: []'))
#         myBinSizeIrrepName = myIrreps[myBinSizeIrrep[0]]
#         myBinLocation = vfl.DIRECTORY_EXISTS(f'{myLocation}Bin_Size_Analysis_{myBinSizeIrrepName}_{myTypeRs}/')
#         
#         bs.BinSizeAnalysis(myArchivo, myBinLocation, myTypeRs, myBinSizeIrrep, myIrreps, myTMinFit, myMaxBinSize, myBinSizeFitRange, myChosenBinSize, myVersion, myWeight, one_tmin=myOneTMin, type_fit=myTypeFit, type_correlation=myTypeCorrelation, kbt=myKbt, number_cfgs=myCnfgs, own_kbt_list=myKbtSamples, isospin_label=ed.ensembles[myEns]['iso_label'], ensemble=myEns, plots_only=runPlotBinOnly)
#     
    myArchivo.close()
    myCorrelator.close()
    
        
## Multi-Hadron correlators
elif myWhichCorrelator=='m':  
    myIsospin = myRuns.isospin
    myNonInteractingLevels = ed.ensembles[myEns][myIsospin]['nonInteractingLevels']
    myChosenIsospin = ed.ensembles[myEns][myIsospin]['iso_tag']
    myVersion =  f'_{myEns}_{myChosenIsospin}_test' 
    myArchivo = h5py.File(ed.ensembles[myEns][myIsospin]['fm'], 'r')
    myIrreps = list(myArchivo.keys())
    
    ### Correlators analysis
    if runCorrs: 
        locationWorkedCorrelators = cs.MultiCorrelatorAnalysis(myArchivo, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps, own_kbt_list=myKbtSamples, first_irrep = myFirstIrrep , last_irrep = myLastIrrep)
    else:
        locationWorkedCorrelators = f'{myLocation}Matrix_correlators_{myTypeRs}{reBin}_v{myVersion}.h5'
    
    try:
        myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    except FileNotFoundError:
        print('Cannot find the correlator file for the further analysis steps.')
        sys.exit()
    
    ### GEVP calculation
    if runEigenvals or runRowsCols:
        mySorting = myRuns.gevp.sorting
        myTD = myRuns.gevp.td
        myRsSorting = myRuns.gevp.rs_sorting
        
        myT0Min = myRuns.t0min
        myT0Max = myRuns.t0max
    
    if runEigenvals:
        evs.EigenvaluesExtraction(myCorrelator, myTypeRs, myIrreps, myT0Min, myT0Max, sorting = mySorting, the_td = myTD, rs_sorting = myRsSorting)
    
    ### Operators Analysis
    if runRowsCols:
        myOperatorAnalysisMethod = myRuns.ops_run.method
        myChosenOpsList = ed.ensembles[myEns][myIsospin]['operatorsChoice']
        rcs.OperatorsAnalysis(myCorrelator, myTypeRs, myOperatorAnalysisMethod, myIrreps, myT0Min, myT0Max, ops_analysis_list = myChosenOpsList)
    
    ### Effective Masses analysis
    if runEffMass: 
        myGevpFlag = myRuns.fit.gevp_flag
        myOperatorsFlag = myRuns.ops
        myDiagonalCorrs = myRuns.diag_corr
        efs.MultiCorrelatorEffectiveMass(myCorrelator, myTypeRs, dist_eff_mass = myEffMassDistance, diag_corrs= myDiagonalCorrs, gevp=myGevpFlag, ops_analysis=myOperatorsFlag)

    ### Fits analysis
    if runFits:        
        myGevpFlag = myRuns.fit.gevp_flag
        myOperatorsFlag = myRuns.ops
        myT0 = myRuns.fit.t0
        myTypeFit = myRuns.fit.type_fit
        myTMaxList = ed.ensembles[myEns][myIsospin]['multiTMaxFits']
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myLocation}Fits_Matrices/')
        myFitCorrelator = h5py.File(f'{myFitsLocation}Matrix_correlators_{myTypeRs}{reBin}_fits_v{myVersion}.h5', 'a')
        
        fts.FitMultiCorrelators(myCorrelator, myFitCorrelator, myTypeRs, myTMaxList, myIrreps,  type_fit = myTypeFit, type_correlation = myTypeCorrelation, one_tmin = myOneTMin, one_t0 = myOneT0, chosen_t0 = myT0, gevp=myGevpFlag, operators_analysis = myOperatorsFlag, the_operator_analysis_method = myOperatorAnalysisMethod, first_irrep = myFirstIrrep , last_irrep = myLastIrrep)
        
        myFitCorrelator.close()
    
    ### Still under construction
    # if runCorrDiff:
        # myNonInteractingFileDict = ed.ensemble_dictionary[myEns]['singleMesons']
        # myFitsLocation = vfa.DIRECTORY_EXISTS(myAnalysisLocation + 'Fits_Matrices/')
        # myFitCorrelator = h5py.File(
            # myFitsLocation + f'Matrix_correlators_{myTypeRs}{reBin}_fits_v{myVersion}.h5', 'r')
        # myCorrDiffLocation = vfa.DIRECTORY_EXISTS(myAnalysisLocation + 'Correlated_Differences/')
        # cds.CORRELATED_DIFFERENCES(myNonInteractingLevels, myTMinPlot, myTMinMPlot, myNonInteractingFileDict, myFitCorrelator, myCorrDiffLocation, myTMinFit, myLatSize,myNrIrreps, myNrOps, myT0, myTD, myTypeRs, myCorrDiffPlot, myThresholdMasses, myEns, myIsoLabel, myTwoBody, myThreeBody, fmToMev, myRb)
    
    myArchivo.close()
    myCorrelator.close()


## Ratio Multi-Hadrons (THIS PART IS STILL UNDER CONSTRUCTION)
elif myWhichCorrelator=='mr':     
    myIsospin = myRuns.isospin
    myNonInteractingLevels = ed.ensembles[myEns][myIsospin]['nonInteractingLevels']
    myChosenIsospin = ed.ensembles[myEns][myIsospin]['iso_tag']
        
    ### Ratio of Correlators
    if runCorrs:
        
        ### Single-hadron corrs
        myVersion =  f'_{myEns}_singles_test' 
        myArchivo = h5py.File(ed.ensembles[myEns]['fs'], 'r')
        mySingleIrreps = list(myArchivo.keys())
        myArchivo.close()
        mySHCorrelator = h5py.File(f'{myLocation}Single_correlators_{myTypeRs}{reBin}_v{myVersion}.h5', 'r+')
        
        
        ### Multi-hadron corrs
        myVersion =  f'_{myEns}_{myChosenIsospin}_test' 
        myArchivo = h5py.File(ed.ensembles[myEns][myIsospin]['fm'], 'r')
        myMultiIrreps = list(myArchivo.keys())
        myArchivo.close()
        myMHCorrelator = h5py.File(f'{myLocation}Matrix_correlators_{myTypeRs}{reBin}_v{myVersion}.h5', 'r+')
        
        ### Taking the ratios
        locationWorkedCorrelators = cs.MultiCorrelatorAnalysisRatios(myMHCorrelator, mySHCorrelator, myLocation, myNonInteractingLevels, myVersion, myTypeRs, myMultiIrreps, mySingleIrreps, rebin_on = myRebinOn, rb = myRb, nr_irreps = myNrIrreps, first_irrep = myFirstIrrep , last_irrep = myLastIrrep)
        
        myMHCorrelator.close()
        mySHCorrelator.close()
    else:
        locationWorkedCorrelators = f'{myLocation}Ratio_matrix_correlators_{myTypeRs}{reBin}_v{myVersion}.h5'
        
    myRatioCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    
    ###(NOT WORKED OUT YET UNTIL HERE, ALSO NOT TESTED)
    ### Effective Masses analysis
    if runEffMass: 
        # efs.MultiCorrelatorEffectiveMass(myRatioCorrelator, myTypeRs)
        efs.RatioMultiCorrelatorEffectiveMass(myRatioCorrelator, myTypeRs)
    
    ### Fits analysis
    if runFits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_Ratios/')
        myFitRatioCorrelator = h5py.File(myFitsLocation + 'Ratio_matrix_correlators' + myChosenIsospin + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        fs.FitRatioMultiCorrelators(myRatioCorrelator, myFitRatioCorrelator, myTypeRs, listTMaxMultiHads, myIrreps, type_fit = myTypeFit, type_correlation = myTypeCorrelation, one_tmin = myOneTMin, one_t0 = myOneT0, chosen_t0 = myT0)
        
        myFitRatioCorrelator.close()
    
    myRatioCorrelator.close()
    
    
else: 
    print('Not proper choice.')
    sys.exit()


### --------- PRINTS WHERE IT IS SAVED --------        

print('-'*(len(locationWorkedCorrelators)+1))
print('Correlator analysis saved as: \n' + locationWorkedCorrelators )
print('_'*(len(locationWorkedCorrelators)+1))


### -------- CLOSES ALL OTHER FILES --------



