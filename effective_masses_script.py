import numpy as np
import h5py
import time
import sys
import os
import set_of_analysis_functions as vfa
import set_of_layout_functions as vfl
import warnings
warnings.filterwarnings('ignore')


def SingleCorrelatorEffectiveMass(the_single_correlator_data, the_type_rs,**kwargs):   
    
    print("                     EFFECTIVE MASSES COMPUTATION \n")
    
    ### Defining distance between time-slice elements of the correlator
    the_dist_eff_mass = int(kwargs.get('dist_eff_mass', 1))
    
    ### The irreps
    the_list_name_irreps = list(the_single_correlator_data.keys())
    
    begin_time = time.time()
       
    for j in range(len(the_list_name_irreps)):
        
        ### Extracting data from file
        this_data = the_single_correlator_data[the_list_name_irreps[j]]
        
        ### Mean values of the real part of the correlator to get the effective masses
        the_mean_corr_real = np.asarray(this_data['Correlators/Real/Mean'])
        
        ### The real part of the resampled data
        the_rs_real = np.transpose(np.asarray(this_data['Correlators/Real/Resampled']), (1,0))
        
        ### Effective Mass computation
        the_em_rs_f =  vfa.EFF_MASS(the_mean_corr_real, the_dist_eff_mass)
            
        ### Loop over the resampled data            
        the_em_rs = np.transpose(np.array([vfa.EFF_MASS(the_rs_real.real[l],the_dist_eff_mass) for l in range(len(the_rs_real))]), (1,0))
        
        ### Mean values of the resamples
        the_mrs_f_real_rs = np.mean(the_em_rs, axis=1)
        
        ### Sigma values for the resamples
        the_sigma_eff_mass = vfa.STD_DEV_MEAN(the_em_rs, the_mrs_f_real_rs, the_type_rs)
        
        ### If the Effective Masses were already computed, this part gets deleted and created a new branch.
        if 'Effective_masses' in this_data.keys(): 
            del the_single_correlator_data[f'{the_list_name_irreps[j]}/Effective_masses']
        
        group_em = this_data.create_group('Effective_masses')
        group_em.create_dataset('Mean', data = the_em_rs_f)
        group_em.create_dataset('Sigmas', data = the_sigma_eff_mass)
        
        print(f'Irrep nr.: {j+1} out of {len(the_list_name_irreps)}')
    end_time = time.time()
    print(f'TOTAL TIME TAKEN: {(end_time-begin_time)/60} mins')
    
    

def SingleCorrelatorEffectiveMassIB(the_single_correlator_data, the_type_rs,**kwargs):   
    
    print("                  EFFECTIVE MASSES COMPUTATION WITH IB\n")
    
    ### Defining distance between time-slice elements of the correlator
    the_dist_eff_mass = int(kwargs.get('dist_eff_mass', 1))
    
    ### The irreps
    the_list_name_irreps = list(the_single_correlator_data.keys())
    
    begin_time = time.time()     
    ### Extracting data from file
    this_data = the_single_correlator_data[the_list_name_irreps[0]]
    
    ### The real part of the resampled data
    the_rs_real = np.transpose(np.asarray(this_data['Correlators/Real/Resampled']), (1,0))
    
    ### Mean values of the real part of the correlator to get the effective masses
    the_mean_corr_real = np.asarray(this_data['Correlators/Real/Mean'])
    
    ### Effective Mass computation
    the_em_rs_f =  vfa.EFF_MASS(the_mean_corr_real, the_dist_eff_mass)
    
    ### Loop over the resampled data            
    the_em_rs = np.transpose(np.array([vfa.EFF_MASS(the_rs_real.real[l],the_dist_eff_mass) for l in range(len(the_rs_real))]), (1,0))
    
    the_mrs_f_real_rs = np.mean(the_em_rs, axis=1)
    
    ### Sigma values for the resamples
    the_sigma_eff_mass = vfa.STD_DEV_MEAN(the_em_rs, the_mrs_f_real_rs, the_type_rs)
    
    ### CORRECTIONS
    this_data_mass = the_single_correlator_data[the_list_name_irreps[1]]
    this_data_qed = the_single_correlator_data[the_list_name_irreps[2]]
    
    ### mean vals
    the_mean_corr_real_mass = np.asarray(this_data_mass['Correlators/Real/Mean'])
    the_mean_corr_real_qed = np.asarray(this_data_qed['Correlators/Real/Mean'])
            
    ### resamples
    the_rs_real_mass = np.transpose(np.asarray(this_data_mass['Correlators/Real/Resampled']), (1,0))
    the_rs_real_qed = np.transpose(np.asarray(this_data_qed['Correlators/Real/Resampled']), (1,0))
    
    ### Effective Mass computation
    the_em_rs_f_mass =  vfa.EFF_MASS_CORRECTIONS(the_mean_corr_real_mass, the_mean_corr_real)
    the_em_rs_f_qed =  vfa.EFF_MASS_CORRECTIONS(the_mean_corr_real_qed, the_mean_corr_real)

    ### Loop over the resampled data            
    the_em_rs_mass = np.transpose(np.asarray([vfa.EFF_MASS_CORRECTIONS(the_rs_real_mass.real[l],the_rs_real[l]) for l in range(len(the_rs_real))]), (1,0))
    
    the_em_rs_qed = np.transpose(np.asarray([vfa.EFF_MASS_CORRECTIONS(the_rs_real_qed.real[l],the_rs_real[l]) for l in range(len(the_rs_real))]), (1,0))
    
    the_mrs_f_real_rs_mass = np.mean(the_em_rs_mass, axis=1)
    the_mrs_f_real_rs_qed = np.mean(the_em_rs_qed, axis=1)
    
    the_sigma_eff_mass_ms = vfa.STD_DEV_MEAN(the_em_rs_mass, the_mrs_f_real_rs_mass, the_type_rs)
    the_sigma_eff_mass_qed = vfa.STD_DEV_MEAN(the_em_rs_qed, the_mrs_f_real_rs_qed, the_type_rs)
    
    the_effs_masses = [the_em_rs_f, the_em_rs_f_mass, the_em_rs_f_qed]
    the_sigmas = [the_sigma_eff_mass, the_sigma_eff_mass_ms, the_sigma_eff_mass_qed]
    
    for j in range(len(the_list_name_irreps)):
        this_data = the_single_correlator_data[the_list_name_irreps[j]]
        ### If the Effective Masses were already computed, this part gets deleted and created a new branch.
        if 'Effective_masses' in this_data.keys(): 
            del the_single_correlator_data[f'{the_list_name_irreps[j]}/Effective_masses']
        
        group_em = this_data.create_group('Effective_masses')
        
        group_em.create_dataset('Mean', data = the_effs_masses[j])
        group_em.create_dataset('Sigmas', data = the_sigmas[j])
        
        print(f'Irrep nr.: {j+1} out of {len(the_list_name_irreps)}')
    end_time = time.time()
    print(f'TOTAL TIME TAKEN: {(end_time-begin_time)/60} mins')
            
            
def MultiCorrelatorEffectiveMass(the_matrix_correlator_data, the_type_rs, **kwargs):
    
    print("                     EFFECTIVE MASSES COMPUTATION \n")
    
    ### Defining distance between time-slice elements of the correlator
    the_dist_eff_mass = int(kwargs.get('dist_eff_mass', 1))
    
    ### The irreps
    the_list_name_irreps = list(the_matrix_correlator_data.keys())    
    
    begin_time = time.time()
    for j in range(len(the_list_name_irreps)): 
        
        ### Extracting data from file
        this_data = the_matrix_correlator_data[the_list_name_irreps[j]]
        
        ### Getting the operators list and the time interval
        the_op_list, nt = list(this_data['Operators']), np.asarray(this_data['Time_slices'])
        
        ### This is the size of the correlator matrices
        the_size_matrix = len(the_op_list)
        
        ### Mean values of the real part of the correlator to get the effective masses
        the_mean_corr_real = np.asarray(this_data['Correlators/Real/Mean'])

        ### The real part of the resampled data
        the_rs_real = np.asarray(this_data['Correlators/Real/Resampled'])
        
        ### Reshaping data
        the_reshaped_mean_corr = vfa.RESHAPING_CORRELATORS(the_mean_corr_real)
        the_reshaped_rs_corr = vfa.RESHAPING_CORRELATORS_RS(the_rs_real)
        
        ### Loop over the nr. of operators = size of the correlator matrix            
        print("Effective Masses of correlator diagonal in process...")
        
        the_ncnfgs = len(the_rs_real[0])
        the_efm_mass = np.asarray([vfa.EFF_MASS(the_reshaped_mean_corr[ii, ii],the_dist_eff_mass) for ii in range(the_size_matrix)])
        
        the_sigma_efm = []
        for ii in range(the_size_matrix):
            
            ### Loop over the resamples                
            the_rs_eff = vfa.NT_TO_NCFGS(np.asarray([vfa.EFF_MASS(the_reshaped_rs_corr[ii, ii, xyz], the_dist_eff_mass) for xyz in range(the_ncnfgs)]))
            
            ### MEan value of the resampled data to compute the sigma vals.
            the_rs_mean_eff = vfa.MEAN(the_rs_eff)
            the_sigma_efm.append(vfa.STD_DEV_MEAN(the_rs_eff, the_rs_mean_eff, the_type_rs))
                
            ### If the branch Effective Masses exists in the file, then it gets deleted and created a new one with the new values.
        the_group_real = this_data['Correlators/Real']
        
        if 'Effective_masses' in the_group_real.keys(): 
            del the_matrix_correlator_data[f'{the_list_name_irreps[j]}/Correlators/Real/Effective_masses']
        the_group_real.create_dataset('Effective_masses', data = np.asarray(the_efm_mass))
        
        if 'Effective_masses_sigmas' in the_group_real.keys(): 
            del the_matrix_correlator_data[f'{the_list_name_irreps[j]}/Correlators/Real/Effective_masses_sigmas']
        the_group_real.create_dataset('Effective_masses_sigmas', data = np.asarray(the_sigma_efm))
            
        ### If the GEVP was performed, then it will identify this section.
        if 'GEVP' in this_data.keys():
            
            gevp_group = this_data.get('GEVP')
                
            print("Effective Masses of GEVP eigenvalues in process...")
            
            vfa.DOING_EFFECTIVE_MASSES_EIGENVALUES(gevp_group, the_dist_eff_mass, the_type_rs)
            
            print(f'Irrep nr.: {j+1} out of {len(the_list_name_irreps)}')
    end_time = time.time()
    print(f'TIME TAKEN: {(end_time-begin_time)/60} mins')
    


def RatioMultiCorrelatorEffectiveMass(the_matrix_correlator_data, the_type_rs, **kwargs):
    
    print("                     EFFECTIVE MASSES COMPUTATION \n")
    
    ### Defining distance between time-slice elements of the correlator
    the_dist_eff_mass = int(kwargs.get('dist_eff_mass', 1))
    
    ### The irreps
    the_list_name_irreps = list(the_matrix_correlator_data.keys())    
    
    begin_time = time.time()
    for j in range(len(the_list_name_irreps)): 
        
        ### Extracting data from file
        this_data = the_matrix_correlator_data[the_list_name_irreps[j]]
        
        ### Getting the operators list and the time interval
        the_op_list, nt = list(this_data['Operators']), np.asarray(this_data['Time_slices'])
        
        ### This is the size of the correlator matrices
        the_size_matrix = len(the_op_list)
                        
        gevp_group = this_data['GEVP']
        
        print("Effective Masses of GEVP eigenvalues in process...")
        
        for item in gevp_group.keys():
            if 'Effective_masses' in gevp_group.get(item).keys(): 
                del gevp_group[f'{item}/Effective_masses']

            group_em_t0 = gevp_group.get(item).create_group('Effective_masses')
            the_evalues_rs_f = np.asarray(gevp_group[f'{item}/Eigenvalues/Resampled'])
            the_evalues_mean_f = np.asarray(gevp_group[f'{item}/Eigenvalues/Mean'])
            
            ### Loop over the total number of eigenvalues
            the_eff_mass_mean , the_cov_eff_mass = [], []
            for ls in range(the_evalues_mean_f.shape[0]):
                
                the_eff_mass_mean_i = []
                the_cov_eff_mass_i  = []
                
                for nn in range(the_evalues_mean_f.shape[1]):
                    
                    the_average = np.array(vfa.EFF_MASS(the_evalues_mean_f[ls,nn], the_dist_eff_mass),dtype=np.float64)
                    the_eff_mass_mean_i.append(the_average)
                    
                    ### Loop over the resamples of the eigenvalues
                    the_eff_mass_rs = [np.asarray(vfa.EFF_MASS(the_evalues_rs_f[ls,nn,zz], the_dist_eff_mass), dtype=np.float64) for zz in range(the_evalues_rs_f.shape[2])]
                
                    ### Reshaping the data
                    the_eff_mass_rs = np.asarray(vfa.NCFGS_TO_NT(the_eff_mass_rs))
                
                    ### Here the statistical errors of the resampled data are computed
                    the_eff_rs_mean = vfa.MEAN(np.array(the_eff_mass_rs))
                    the_cov_eff_mass_i.append(vfa.STD_DEV_MEAN(the_eff_mass_rs, the_eff_rs_mean, the_type_rs))
                
                the_eff_mass_mean.append(np.asarray(the_eff_mass_mean_i))
                the_cov_eff_mass.append(np.asarray(the_cov_eff_mass_i)) 
            
            group_em_t0.create_dataset('Mean', data=np.asarray(the_eff_mass_mean))
            group_em_t0.create_dataset('Sigmas',data=np.asarray(the_cov_eff_mass))
        
        print('Irrep nr.: '+ str(j+1) + ' out of ' +str(len(the_list_name_irreps)))
    end_time = time.time()
    print('TIME TAKEN: ' + str((end_time-begin_time)/60) +' mins')


### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


if __name__=="__main__":
    import ensembles as ed
    
    myArgs = vfl.parse_args()
    myRuns = vfl.WhichRuns(myArgs, ed.ensembles)
    
    reBin = f"_bin{myRuns.rb}" if myRuns.rebin else ""
    
    myLocation = vfl.DIRECTORY_EXISTS(f'{ed.outputLocation}{myRuns.ensemble}/')
    
    vfl.INFO_PRINTING(myWhichCorrelator, myEns)
    
    if myRuns.correlator =='s':
        if not myRuns.ib_corr:
            myVersion =  f'{myRuns.ensemble}_singles_test' 
            myEffMassFunction = SingleCorrelatorEffectiveMass
        else:
            myVersion =  f'{myRuns.ensemble}_omega_test' 
            myEffMassFunction = SingleCorrelatorEffectiveMassIB
            
        mySingleCorrelatorData = h5py.File(f'{myLocation}Single_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5', 'r')
        myIrreps = list(mySingleCorrelatorData.keys())  
        
        myEffMassFunction(mySingleCorrelatorData,  myRuns.rs_type, dist_eff_mass = myRuns.dist_eff_mass) 
        mySingleCorrelatorData.close()
    
    elif myWhichCorrelator=='m':    
        myVersion =  f'_{myEns}_{myChosenIsospin}_test' 
        myNameArchivo = f'{myLocation}Matrix_correlators_{myRuns.rs_type}{reBin}_v{myVersion}.h5'
        
        myMatrixCorrelatorData = h5py.File(myNameArchivo, 'r+')
        MultiCorrelatorEffectiveMass(myMatrixCorrelatorData, myRuns.rs_type, dist_eff_mass = myRuns.dist_eff_mass, diag_corrs= myDiagonalCorrs, gevp=myGevpFlag, ops_analysis=myOperatorsFlag)
        myMatrixCorrelatorData.close()
    
    elif myWhichCorrelator=='mr':
        myNameArchivo = f'{myLocation}Matrix_correlators_ratios_{myRuns.rs_type}{reBin}_v{myVersion}.h5'
        myRatioMatrixCorrelatorData = h5py.File(myNameArchivo, 'r+')
        MultiCorrelatorEffectiveMass(myRatioMatrixCorrelatorData, myRuns.rs_type)
        myRatioMatrixCorrelatorData.close()
    
    print('-'*(len(myNameArchivo)+1))
    print('Saved as: \n' + myNameArchivo)
    print('_'*(len(myNameArchivo)+1))
