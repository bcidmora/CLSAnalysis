import numpy as np
from scipy.linalg import eigh
from scipy.linalg import fractional_matrix_power
import h5py
import time
import sys
import os
import set_of_analysis_functions as vfa

def EigenvaluesExtraction(the_matrix_correlator_data, the_type_rs, the_irreps, the_t0_min, the_t0_max, **kwargs):   
    
    print("                     SOLVING GEVP \n")
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

    begin_time = time.time()       
    for the_irrep in the_m_irreps:
        
        ### The data to analise is extracted here
        this_data = the_matrix_correlator_data[the_irrep]
        
        ### The list of operators of the correlation matrix and the time slices.
        the_op_list, the_nt = list(this_data['Operators']), np.asarray(this_data['Time_slices'])
        
        ### The size of each correlation matrix
        the_size_matrix = len(the_op_list)
        
        ### The resampled correlators
        the_rs_real = np.asarray(this_data['Correlators/Real/Resampled'])
        
        ### The central values of the original correlators
        the_mean_corr_real = np.asarray(this_data['Correlators/Real/Mean'])    
        the_mean_corr = np.asarray(the_mean_corr_real, dtype = np.float64)
        
        print('\n----------------------------------------------')
        print(f'     IRREP ({the_irreps.index(the_irrep)+1}/{len(the_irreps)}): {the_irrep}')
        print(f'Size of the Correlation matrix: {the_size_matrix}x{the_size_matrix}\nTime slices: {the_nt[0]} - {the_nt[-1]}\nResampling data ({the_resampling_scheme}): {the_rs_real.shape[1]}\n----------------------------------------------')
        print('      OPERATORS LIST \n----------------------------------------------')
        
        for i in range(the_size_matrix):
            print(f'       {the_op_list[i].decode('utf-8')}')
        
        if 'GEVP' in this_data.keys(): 
            del the_matrix_correlator_data[f'{the_irrep}/GEVP']
        group_gevp = this_data.create_group('GEVP')
        
        if kwargs.get('the_td')==None:
            ### This is a loop over the t0s
            vfa.DOING_THE_GEVP([the_t0_min, the_t0_max], the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, the_rs_sorting_process, group_gevp)
        else:
            vfa.DOING_THE_GEVP_SINGLE_PIVOT([the_t0_min, the_t0_max], the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, the_rs_sorting_process, group_gevp,t_diag=int(kwargs.get('the_td')))
    end_time = time.time()
    print(f'TIME TAKEN: {(end_time-begin_time)/60} mins')


### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


if __name__=="__main__":
    
    ### Ensemble you want to analyse
    myEns = str(sys.argv[1]).upper()
    
    ### Type of resampling 'bt' or 'jk'
    myTypeRs = str(sys.argv[2]).lower()
    
    ### Rebinning
    myRebinOn = str(sys.argv[3]).lower()
    myRb = 10#sys.argv[4]
    
    ### How the eigenstates are sorted
    mySorting = 'eigenvals' #'eigenvals' # 'vecs_fix' # 'vecs_fix_norm' # 'vecs_var' # 'vecs_var_norm'
    
    ### Name of this version of analysis
    myVersion = 'test'
    
    ### Min and Max t0 to do the GEVP
    myT0Min = int(input('T0 min: '))
    myT0Max = int(input('T0 max: ')) 
    
    ### Root Location of your hdf5 file that contains the correlators already resampled and averaged.
    myLocation = os.path.expanduser('~')+'/Documents/Chris Files/CorrelatorData/%s/'%myEns
    
    ### Rebin naming
    if myRebinOn=='rb': 
        rb = int(myRb)
        reBin = '_bin'+str(rb)
    else:
        reBin = ''  
    
    if myEns == 'N451': from files_n451 import name
    elif myEns == 'N201': from files_n201 import name 
    elif myEns == 'D200': from files_d200 import name
    elif myEns == 'X451': from files_x451 import name
    
    import ensembles as ed
    
    myIrreps = name
    
    ### This is the file that cotains the averaged correlators and stuff
    # myArchivo = myLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion
    myArchivo = myLocation + 'Matrix_correlators_bt_bin10_v8_old_oficial_new_gevp.h5'
    
    myMatrixCorrelatorData = h5py.File(myArchivo,'r+')
    
    EigenvaluesExtraction(myMatrixCorrelatorData, myTypeRs, myIrreps, t0_min = myT0Min, t0_max = myT0Max, sorting=mySorting)
    # EigenvaluesExtraction(myCorrelator, myTypeRs, myIrreps, t0_min = myT0Min, t0_max = myT0Max, sorting=mySorting,the_td=myTD)
    myMatrixCorrelatorData.close()   
    
    print('-'*(len(myArchivo)+1))
    print('Saved as: \n' + myArchivo)
    print('_'*(len(myArchivo)+1))
