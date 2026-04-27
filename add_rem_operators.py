import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_functions as vf

def OperatorsAnalysis(the_matrix_correlator_data, the_type_rs, the_operator_analysis_method, the_irreps, **kwargs):
    
    print("                     OPERATORS ANALYSIS PROCESS \n")
    ### The list of total irreps
    the_m_irreps =  list(the_matrix_correlator_data.keys())
    
    ### Resampling scheme
    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
    
    ### Getting the t0 min and t0 max to do the GEVP
    if (the_t0_min is None or the_t0_max is None) and kwargs.get('the_td') is None:
        sys.exit('Error: T0 min or T0 max not valid. \nQuitting.')
    
    ### What type of sorting of the eigenstates
    the_sorting = kwargs.get('sorting')
    the_sorting_map = {
        None : (vfa.SORTING_EIGENVALUES_NEW, "Sorting states based on Eigenvalues."),
        'eigenvals' : (vfa. SORTING_EIGENVALUES_NEW, "Sorting states based on Eigenvalues."),
        'vecs_fix' : (vfa.SORTING_EIGENVECTORS, "Sorting states by Eigenvectors with a fixed reference time slice."),
        'vecs_fix_norm' : (vfa.SORTING_EIGENVECTORS_NORMALIZED, "Sorting states by normalized Eigenvectors with a fixed reference time slice."),
        'vecs_var' : (vfa.SORTING_EIGENVECTORS_CHANGING_TSLICE, "Sorting states by Eigenvectors with a varying reference time slice."),
        'vecs_var_norm' : (vfa.SORTING_EIGENVECTORS_NORMALIZED_CHANGING_TSLICE, "Sorting states by normalized Eigenvectors with a varying reference time slice." ), 
        'vecs_var_rs_mean' : (vfa.SORTING_EIGENVECTORS_RS_MEAN, "Sorting states by Eigenvectors with a varying reference time slice. The resamples are sorted in the same way as the mean."),}
    
    the_sorting_process, the_msg = the_sorting_map.get(the_sorting, the_sorting_map[None])
    print(the_msg)
    
    the_rs_sorting = kwargs.get('rs_sorting')
    if the_rs_sorting is None:
        the_rs_sorting_process = vfa.SORTING_EIGENVALUES_NEW
    else:
        the_rs_sorting_process = the_sorting_process
        
    ### If not all irreps are analysed, this number can be changed.
    
    the_operator_analysis_map = { 
        'adding': {'method': vfa.ADD_ROWS_COLS, 'start': 1, 'end': -1, 'name': 'Add_', 'msg': "Adding operators to the correlation matrix one by one starting with a 2x2 matrix."},
        'removing': {'method': vfa.REMOVE_ROWS_COLS, 'start': 0, 'end': 0, 'name': 'Remove_', 'msg': "Removing operators from the correlation matrix one at the time."},
        'from_list': {'method': vfa.CHOOSE_OPS, 'start': None, 'end': None, 'name': "Ops_chosen_", 'OpList': kwargs.get('ops_analysis_list'), 'msg':  "Using a custom list of operators"},
        }
    
    if the_operator_analysis_method not in the_operator_analysis_map:
        raise ValueError(f"Unknown operator analysis method: {the_operator_analysis_method}")
    
    da_method = the_operator_analysis_map[the_operator_analysis_method]
    
    the_op_method = da_method['method']
    the_start_op = da_method['start']
    the_last_op = da_method['end']
    the_name_method = da_method['name']
    print(da_method['msg'])
    
    the_chosen_op_list = da_method.get('OpList')        
        
    begin_time = time.time()
    for the_irrep in the_m_irreps:
        
        ### The data to analyse. 
        this_data = the_matrix_correlator_data[the_irrep]
        
        ### The operators list and the time slices
        the_op_list, the_nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        
        if 'Operators_Analysis' not in this_data.keys(): 
            the_group_rows_cols = this_data.create_group('Operators_Analysis')
        else:
            the_group_rows_cols = this_data.get('Operators_Analysis')
        
        this_data_corr = np.asarray(this_data['Correlators/Real/'])
        
         ### This is just reshaping the data to make it easier to analyse
        the_mod_data = vfa.RESHAPING_CORRELATORS(this_data_corr['Mean'])
        the_mod_data_rs = vfa.RESHAPING_CORRELATORS_RS_NT(this_data_corr['Resampled'])
        
        print('\n----------------------------------------------')
        print(f'     IRREP ({the_irreps.index(the_irrep)+1}/{len(the_irreps)}): {the_irrep}')
        
        ### If you want to choose from a list, then you have to provide the list.
        if the_operator_analysis_method == 'from_list':
            the_chosen_op_list_j = the_chosen_op_list[the_irreps.index(the_irrep)]
            
            if not the_chosen_op_list_j or len(the_chosen_op_list_j) == len(the_op_list):
                continue

            the_selections = [(the_chosen_op_list_j, the_name_ops_analysis + ''.join(map(str, the_chosen_op_list_j)))]
            the_ops_chosen_string = f"{the_name_ops_analysis}{'-'.join(map(str, sorted(the_chosen_op_list_j)))}"
        else:
            the_indices = np.arange(the_start_op, len(the_op_list) + the_last_op)            
            the_selections = [(ii, f"{the_name_ops_analysis}Op_{ii}")for ii in the_indices]
        
        for the_sel, the_group_name in the_selections:
            if the_group_name in the_group_rows_cols:
                del the_matrix_correlator_data[f"{the_irrep}/Operators_Analysis/{the_group_name}"]

            group_i = the_group_rows_cols.create_group(the_group_name)

            the_new_corrs = the_op_method(the_mod_data, the_mod_data_rs, the_sel)

            the_mean_corr = np.asarray(the_new_corrs[0], dtype=np.float64)
            the_rs_real = np.asarray(the_new_corrs[1], dtype=np.float64)

            the_mean_corr = vfa.RESHAPING_EIGENVALS_MEAN(the_mean_corr)
            the_rs_real = vfa.RESHAPING_EIGENVALS_RS(the_rs_real)

            print(f"Size of the Correlation matrix: {the_mean_corr.shape[-1]}x{the_mean_corr.shape[-1]}\nTime slices: {the_nt[0]} - {the_nt[-1]}\nResampling data: {the_resampling_scheme} {the_rs_real.shape[1]}\n----------------------------------------------")

            vfa.DOING_THE_GEVP([the_t0_min, the_t0_max], the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, the_rs_sorting_process, group_i)
        
    end_time = time.time()
            


### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


if __name__=="__main__":
    import ensembles as ed
    
    ### Ensemble you want to analyse
    myEns = str(sys.argv[1]).upper()
    
    ### Type of resampling 'bt' or 'jk'
    myTypeRs = str(sys.argv[2]).lower()
    
    ### Rebinning
    myRebinOn = str(sys.argv[3]).lower()
    myRb = 1
    
     ### How the eigenstates are sorted
    mySorting = 'eigenvals' #'eigenvals' # 'vecs_fix' # 'vecs_fix_norm' # 'vecs_var' # 'vecs_var_norm'
    myRsSorting = None # mySorting # None
    myTD = None
    
    myOperatorMethod = 'from_list' # 'adding' # 'removing' # 'from_list'
    
    ### Min and Max t0 to do the GEVP
    myT0Min = int(input('T0 min: '))
    myT0Max = int(input('T0 max: ')) 
    
    myLocation = vfl.DIRECTORY_EXISTS(f'{ed.outputLocation}{myEns}/')
    
    myNrIrreps=None #1 #2
    
    ### Rebin naming
    if myRebinOn=='rb': 
        reBin = f'_bin{myRb}'
    else: 
        reBin = ''
    
    ### Name of this version of analysis
    myVersion =  f'_{myEns}_{myChosenIsospin}_test' 
    
    myArchivo = h5py.File(ed.ensembles[myEns][myIsospin]['fm'], 'r')
    myIrreps = list(myArchivo.keys())
    myArchivo.close()
    myListOperators = ed.ensembles[myEns][myIsospin]['operatorsChoice']
    
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    myArchivo = f'{myLocation}Matrix_correlators_{myTypeRs}{reBin}_v{myVersion}.h5'
    myMatrixCorrelatorData = h5py.File(myArchivo,'r+')
    
    OperatorsAnalysis(myMatrixCorrelatorData, myTypeRs, myOperatorMethod, myIrreps, ops_analysis_list = myListOperators)
    
    myArchivo.close()
    
    
