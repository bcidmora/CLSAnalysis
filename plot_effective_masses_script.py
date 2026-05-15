import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_layout_functions as vfl
import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp

def PlotSingleHadronsEffectiveMasses(the_single_correlator_data, the_rs_scheme, the_version, the_location, the_rebin, **kwargs):
    
    s_irreps = list(the_single_correlator_data.keys())
    
    # ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(s_irreps)
    
    s_irreps = s_irreps[the_first_irrep:the_last_irrep]
    
    ### Loop over the SH irreps
    for the_irrep in s_irreps:
        ### Central values of the correlators
        the_mean_corr = np.array(the_single_correlator_data[f'{the_irrep}/Effective_masses/Mean'])
        
        ### Statistical errors of the correlators
        the_sigmas_corr = np.array(the_single_correlator_data[f'{the_irrep}/Effective_masses/Sigmas'])
        
        ##3 Time extent
        the_nt_corr = np.array(the_single_correlator_data[f'{the_irrep}/Time_slices'])
        the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
        
        the_nt_ticks = np.arange(the_nt_corr[2], the_nt_corr[-1], 2)
        
        ### The SH operator that appears in the plot
        the_op = list(the_single_correlator_data[f'{the_irrep}/Operators'])[0]
        theOperatorNamePlot = vfp.OPERATORS_SH(the_op.decode('utf-8'))
        
        da_irrep = vfp.IrrepInfo(the_irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        nrMomentumIrrep = da_irrep.Momentum
        NameIrrepPlot = da_irrep.NamePlot 
    
        OperatorNamePlot = vfp.SH_OPERATORS_RELABEL(theOperatorNamePlot, NameIrrepPlot, nrMomentumIrrep )
        
        print('Effective Mass plot in progress...')
        # the_efm_fig = plt.figure(figsize=(4,3))
        the_efm_fig = plt.figure()
        vfp.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o',  OperatorNamePlot)
        plt.show()
        the_efm_fig.savefig(f'{the_location}EffectiveMass_{the_irrep}{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
        
def PlotMultiHadronsEffectiveMasses(the_matrix_correlator_data, the_quantum_number, the_rs_scheme, the_version, the_t0, the_location, the_rebin, **kwargs):
    ### Getting all the irreps in this ensemble
    m_irreps = list(the_matrix_correlator_data.keys())
    
    ### These variables are to plot the GEVP or the operators analysis eigenvalues
    the_diagonal_corrs_flag = kwargs.get('diag_corrs')
    all_corr_flag = kwargs.get('all_corr')
    
    ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(m_irreps)
    
    m_irreps = m_irreps[the_first_irrep:the_last_irrep]    
    
    ### Loop over the irreps
    for the_irrep in m_irreps:
        
        ### Operator list of that specific irrep
        the_op_list = list(the_matrix_correlator_data[f'{the_irrep}/Operators'])
        
        ### Time slices 
        the_nt = np.asarray(the_matrix_correlator_data[f'{the_irrep}/Time_slices'])
            
        ### The new time slices range for the effective masses
        the_nt_corr_efm = np.arange(the_nt[0]+0.5, the_nt[-1]+0.5, 1)
        
        ### These are the ticks in the x-axis
        the_nt_ticks = np.arange(the_nt[0]+1, the_nt[-1], int(len(the_nt)/5))
        
        ### Information about the irrep written for the plots
        da_irrep = vfp.IrrepInfo(the_irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot
        
        if the_diagonal_corrs_flag:

            ### Effective masses of correlators
            the_data_corr = np.asarray(the_matrix_correlator_data[f'{the_irrep}/Correlators/Real/Effective_masses'])
            
            ### Statistical error of the effective masses of the central values of the correlators
            the_data_sigmas_corr = np.asarray(the_matrix_correlator_data[f'{the_irrep}/Correlators/Real/Effective_masses_sigmas'])

            ### Loop over the size of the correlation matrix
            for bb in range(len(the_op_list)):
                ### The effective masses of the mean values of the diagonal of the correlators
                the_mean_efm = the_data_corr[bb]
                
                ### Their sigmas
                the_sigmas_efm = the_data_sigmas_corr[bb]
                
                ### The operator of this dataset
                the_op = the_op_list[bb]
                
                ### Convenient name for the plots
                theOperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
                OperatorNamePlot = vfp.MH_OPERATORS_RELABEL(theOperatorNamePlot)
                
                print('Effective mass diagonal correlators plots in process...')
                
                ### Checking for the data
                the_ymin = vfp.CHOOSING_YMIN_PLOT(the_mean_efm)
                
                efm_corr_fig = plt.figure(figsize=(6,4))
                # vfp.PLOT_CORRELATORS(the_nt_corr_efm[1:], the_mean_efm[1:], the_sigmas_efm[1:], the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o', NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to \;C_{%s}$'%(str(bb)+str(bb)) + ' = '+OperatorNamePlot, ymin=the_ymin)
                vfp.PLOT_CORRELATORS(the_nt_corr_efm[1:], the_mean_efm[1:], the_sigmas_efm[1:], the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o', f'{NameIrrepPlot} ({MomentumIrrep}) ' + rf'$\to$ {OperatorNamePlot}', ymin=the_ymin)
                plt.show()
                efm_corr_fig.savefig(f'{the_location}EffectiveMass_DiagonalCorrelators_{the_quantum_number}_{the_irrep}_{bb}{the_rebin}_{the_version}.pdf')
            
            if all_corr_flag:
                ### Here all the diagonal of the correlators are put together
                efm_corr_all_fig = plt.figure()
                print('Effective mass ALL diagonal correlators together plot in process...')
                
                ### Loop over the operators of the correlation matrix
                for bb in range(len(the_op_list)):  
                    ### The diagonal of the correlator
                    the_mean_efm = the_data_corr[bb]
                    the_sigmas_efm = the_data_sigmas_corr[bb]

                    the_op = the_op_list[bb]
                    theOperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
                    OperatorNamePlot = vfp.MH_OPERATORS_RELABEL(theOperatorNamePlot)
                    
                    plt.errorbar(the_nt_corr_efm, the_mean_efm, yerr = the_sigmas_efm, marker=the_markers_list[bb], ls='None', ms=8, markeredgewidth=2, lw=1.5, elinewidth=1.5, capsize=6, label = OperatorNamePlot, color=the_colors[bb])
                plt.xlabel(r'$t\,/\,a$', fontsize=36)
                plt.ylabel(r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=36)
                # plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to\;C_{ii}(t)$'+ ' [' + the_rs_scheme + ']', fontsize=20)
                # plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to\;C_{ii}(t)$', fontsize=20)
                plt.title( f'{NameIrrepPlot} ({MomentumIrrep}) ', fontsize=28)
                plt.xticks(the_nt_ticks,fontsize=18)
                plt.yticks(fontsize=18)
                plt.legend(fontsize=16, ncol=2, handletextpad=0.1, columnspacing=0.1)
                # if len(the_op_list)>6: the_n_cols = int(len(the_op_list)/3)
                # else: the_n_cols = int(len(the_op_list)/2)
                plt.ylim([0.3,1.5])
                # plt.legend(fontsize=14, ncol=the_n_cols, handletextpad=0.3)
                plt.tight_layout()
                plt.show()
                efm_corr_all_fig.savefig(f'{the_location}EffectiveMass_ALLDiagonalCorrelators_{the_quantum_number}_{the_irrep}{the_rebin}_{the_version}.pdf')
            
#             if the_ops_analysis_flag:#and the_chosen_op_list!=None and len(the_chosen_op_list)!=0:
#                 ### Here all the diagonal of the correlators are put together
#                 
#                 the_chosen_op_list = the_chosen_op_list[the_first_irrep:the_last_irrep]
#                 print(the_chosen_op_list)
#                 print('Effective mass CHOSEN diagonal correlators together plot in process...')
#                 
#                 the_temporary_op_list = the_chosen_op_list[m_irreps.index(the_irrep)]
#                 the_label_string= ''
#                 the_colors_index=0
#                 
#                 # efm_corr_all_chosen_fig = plt.figure(figsize=(6.5,5))
#                 efm_corr_all_chosen_fig = plt.figure()
#                 ### Loop over the operators of the correlation matrix
#                 for bb in the_temporary_op_list:  
#                     ### The diagonal of the correlator
#                     the_mean_efm = the_data_corr[bb]
#                     the_sigmas_efm = the_data_sigmas_corr[bb]
# 
#                     the_op = the_op_list[bb]
#                     theOperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
#                     OperatorNamePlot = vfp.MH_OPERATORS_RELABEL(theOperatorNamePlot)
#                     # print(the_op)
#                     # print(theOperatorNamePlot)
#                     # print(OperatorNamePlot)
#                     
#                     plt.errorbar(the_nt_corr_efm[1:], the_mean_efm[1:], yerr = the_sigmas_efm[1:], marker=the_markers_list[the_colors_index], ls='None', ms=8, markeredgewidth=2, lw=1.5, elinewidth=1.5, capsize=6, label = OperatorNamePlot, color=the_colors[the_colors_index])
#                     the_colors_index+=1
#                     the_label_string+=str(bb)
#                 plt.xlabel(r'$t\,/\,a$', fontsize=36)
#                 plt.ylabel(r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=36)
#                 # plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to\;C_{ii}(t)$'+ ' [' + the_rs_scheme + ']',fontsize=20)
#                 # plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to\;C_{ii}(t)$',fontsize=20)
#                 plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep,fontsize=28)
#                 plt.xticks(the_nt_ticks,fontsize=18)
#                 plt.yticks(fontsize=18)
#                 plt.legend(fontsize=16, ncol=2, handletextpad=0.1,columnspacing=0.01)
#                 # if len(the_op_list)>6: the_n_cols = int(len(the_op_list)/3)
#                 # else: the_n_cols = int(len(the_op_list)/2)
#                 plt.ylim([0.42,0.7])
#                 plt.xlim([2.8,19])
#                 # plt.legend(fontsize=14, ncol=the_n_cols, handletextpad=0.3)
#                 plt.tight_layout()
#                 plt.show()
#                 efm_corr_all_chosen_fig.savefig(the_location + 'EffectiveMass_Chosen%s_DiagonalCorrelators'%the_label_string + the_quantum_number + the_irrep + the_rebin + '_v%s.pdf'%the_version)

        ### If GEVP was performed, the eigenvalues are also going to be plotted.
        if 'GEVP' in list(the_matrix_correlator_data[the_irrep].keys()):
            
            the_data = np.asarray(the_matrix_correlator_data[f'{the_irrep}/GEVP/t0_{the_t0}/Effective_masses/Mean'])
            the_data_sigmas = np.asarray(the_matrix_correlator_data[f'{the_irrep}/GEVP/t0_{the_t0}/Effective_masses/Sigmas'])

            ### Loop over the eigenvalues
            for bb in range(len(the_data)):
                the_mean_corr = the_data[bb]
                the_sigmas_corr = the_data_sigmas[bb]
                
                print('Effective mass eigenvalues plots in process...')
                
                 ### Checking for the data
                the_ymin = vfp.CHOOSING_YMIN_PLOT(the_mean_corr)
                
                efm_fig = plt.figure()
                vfp.PLOT_CORRELATORS(the_nt_corr_efm, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o',  f'{NameIrrepPlot} ({MomentumIrrep}) ' + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %sa$)'%str(the_t0), ymin=the_ymin)
                plt.show()
                efm_fig.savefig(f'{the_location}EffectiveMass_Eigenvalues_{the_quantum_number}_{the_irrep}_{bb}_t0_{the_t0}{the_rebin}_{the_version}.pdf')
            
            if all_corr_flag:
                ### All Eigenvalues in only one plot
                efm_corr_all_fig = plt.figure()
                
                print('Effective mass ALL eigenvalues together plot in process...')

                the_min_position = np.where(the_data[0] == min(the_data[0][:-3]))
                the_max_position = np.where(the_data[0] == max(the_data[-1][:-3]))
        
                the_min_y = (the_data[0][the_min_position]-the_data_sigmas[0][the_min_position])*.95
                the_max_y= (the_data[-1][the_max_position]+the_data_sigmas[-1][the_max_position])*1.05
                
                ### Loop over the eigenvalues
                for bb in range(4):#len(the_data)):   
                    the_mean_efm = the_data[bb]
                    the_sigmas_efm = the_data_sigmas[bb]

                    the_op = the_op_list[bb]
                    theOperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
                    OperatorNamePlot = vfp.MH_OPERATORS_RELABEL(theOperatorNamePlot)
                    
                    plt.errorbar(the_nt_corr_efm[1:], the_mean_efm[1:], yerr = the_sigmas_efm[1:], marker=the_markers_list[bb], ls='None', ms=8, markeredgewidth=2, lw=1.5, elinewidth=1.5, capsize=6, label = r'$\lambda_{%s}$'%str(bb), color=the_colors[bb])
                plt.xlabel(r'$t\,/\,a$', fontsize=36)
                plt.ylabel(r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=36)
                plt.title( NameIrrepPlot+ '(%s)'%MomentumIrrep + r' $\to \lambda_{i}(t,t_{0} = %sa)$'%str(the_t0), fontsize=28)
                plt.xticks(the_nt_ticks, fontsize=18)
                plt.yticks(fontsize=18)
                # plt.ylim([the_max_y,the_min_y])
                if len(the_data)>6: the_n_cols = int(len(the_data)/3)
                else: the_n_cols = int(len(the_data)/2)
                plt.ylim([0.42,0.7])
                plt.xlim([2.8,19])
                plt.legend(fontsize=18, ncol=the_n_cols, handletextpad=0.05, columnspacing=0.3)
                plt.tight_layout()
                plt.show()
                efm_corr_all_fig.savefig(the_location + 'EffectiveMass_ALLEigenvalues' + the_quantum_number + the_irrep + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
        
            
            
def PlotRatioHadronsEffectiveMasses(the_ratio_correlator_data, the_quantum_number, the_rs_scheme, the_version, the_t0, the_location, the_rebin):
    mr_irreps = list(the_ratio_correlator_data.keys())
    for irrep in mr_irreps:
        the_op_list = list(the_ratio_correlator_data[irrep+'/Operators'])
        the_data = np.array(the_ratio_correlator_data[irrep + '/GEVP/t0_%s/Effective_masses/Mean'%the_t0])
        the_data_sigmas = np.array(the_ratio_correlator_data[irrep + '/GEVP/t0_%s/Effective_masses/Sigmas'%the_t0])
        for bb in range(len(the_op_list)):
            the_mean_corr = the_data[bb][the_t0-2:]
            the_sigmas_corr = the_data_sigmas[bb][the_t0-2:]
            the_nt_corr = np.array(the_ratio_correlator_data[irrep + '/Time_slices'])[the_t0-2:]
            the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
            the_nt_ticks = np.arange(5, the_nt_corr[-1], 5)

            the_op = the_op_list[bb]
            # OperatorNamePlot = vfp.OPERATORS_SH(the_op.decode('utf-8'))
            OperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
            da_irrep = vfp.IrrepInfo(irrep)
            MomentumIrrep = da_irrep.TotalMomPlot
            NameIrrepPlot = da_irrep.NamePlot
            NameIrrep = da_irrep.Name
            
            print('Effective Mass plot in progress...')
            efm_fig = plt.figure()
            plt.errorbar(the_nt, the_mean_corr, yerr = the_sigmas_corr, marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5, label = '%s'%the_rs_scheme)
            plt.xlabel(r'$t$')
            plt.ylabel(r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$')
            plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0)
            plt.xticks(the_nt_ticks)
            plt.tight_layout()
            #plt.show()
            efm_fig.savefig(the_location + 'EffectiveMass_Eigenvalues_ratios' + the_quantum_number + irrep + '_%s'%bb + the_rebin + '_v%s.pdf'%the_version)




### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


if __name__=="__main__":
    
    myEns = str(sys.argv[1]).upper()
    myWhichCorrelator = str(sys.argv[2]).lower()
    myTypeRs = str(sys.argv[3]).lower()
    myRebinOn = str(sys.argv[4]).lower()
    
    myRb = 2
    myVersion = 'test'
    myT0 = 4
    
    myDataLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/CorrelatorData/%s/'%myEns)
    
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    if myRebinOn=='rb': 
        reBin = '_bin'+str(myRb)
    else:
        reBin = '' 
        
    if myTypeRs=='jk':
        myResamplingScheme='Jackknife'
    elif myTypeRs=='bt':
        myResamplingScheme='Bootstrap'    

    if myWhichCorrelator=='s':
        mySingleCorrelatorData = h5py.File(myDataLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/Plots/%s/SingleHadrons/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotSingleHadronsEffectiveMasses(mySingleCorrelatorData, myResamplingScheme, myVersion, myPlotLocation, reBin)
        
        mySingleCorrelatorData.close()
    
    elif myWhichCorrelator=='m':
        myMatrixCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/Plots/%s/Matrices/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotMultiHadronsEffectiveMasses(myMatrixCorrelatorData, myResamplingScheme, myVersion, myT0, myPlotLocation, reBin)
        
        myMatrixCorrelatorData.close()
        
    elif myWhichCorrelator=='mr':
        myRatioCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/Plots/%s/Matrices/Ratios/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotRatioHadronsEffectiveMasses(myRatioCorrelatorData, myResamplingScheme, myVersion, myT0, myPlotLocation, reBin)
        
        myRatioCorrelatorData.close()
    
