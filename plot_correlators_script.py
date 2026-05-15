import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys

### LAYOUT AND EXTRA FUNCTIONS
import set_of_layout_functions as vfl
import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp

import warnings
warnings.filterwarnings('ignore')

def PlotSingleHadronCorrelators(the_single_correlator_data, the_type_rs, the_version, the_location, the_rebin, **kwargs):
    
    s_irreps = list(the_single_correlator_data.keys())
    
    ### How many irreps do you want to study        
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
    
    if the_type_rs=='jk':
        the_rs_scheme=r'Jackknife'
    elif the_type_rs=='bt':
        the_rs_scheme=r'Bootstrap'   
    
    for irrep in s_irreps:
        ### This is the central values of the correlator
        the_mean_corr = np.asarray(the_single_correlator_data[f'{irrep}/Correlators/Real/Mean'])
        
        ### These are the statistical errors of the correlators
        the_sigmas_corr = np.asarray(the_single_correlator_data[f'{irrep}/Correlators/Real/Sigmas'])
        
        ### This is the time extent
        the_nt = np.asarray(the_single_correlator_data[f'{irrep}/Time_slices'])
        the_nt_ticks = np.arange(the_nt[0]+1, the_nt[-1], int(len(the_nt)/5))

        ### Information about the operators
        the_op_list = list(the_single_correlator_data[f'{irrep}/Operators'])[0]
        theOperatorNamePlot = vfp.OPERATORS_SH(the_op_list.decode('utf-8'))        

        ### Information about the irrep
        da_irrep = vfp.IrrepInfo(irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        nrMomentumIrrep = da_irrep.Momentum
        NameIrrepPlot = da_irrep.NamePlot
        
        OperatorNamePlot = vfp.SH_OPERATORS_RELABEL(theOperatorNamePlot, NameIrrepPlot, nrMomentumIrrep )
        
        print('Correlator plot in progress...')
        the_corr_fig = plt.figure()
        vfp.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$\mathbb{Re}\;C(t)$', 'o', OperatorNamePlot)
        plt.show()
        the_corr_fig.savefig(f'{the_location}Correlator_{irrep[:4]}_{irrep[-1]}{the_rebin}_{the_version}.pdf', bbox_inches='tight')
    
        print('Correlator Log-plot in process...')
        the_log_corr_fig = plt.figure()
        vfp.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$\log \mathbb{Re}\;C(t)$', 'o', OperatorNamePlot, yscale='log')
        plt.show()
        the_log_corr_fig.savefig(f'{the_location}Correlator_{irrep[:4]}_{irrep[-1]}_log{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
        print('Correlator histogram in process...')
        tt = int(len(the_nt)/3)+1
        the_gauss_fig = plt.figure()
        the_nt_mean = the_mean_corr[tt]
        the_rs = np.asarray(the_single_correlator_data[f'{irrep}/Correlators/Real/Resampled'])[tt]
        the_nr_bins = int(len(the_rs)*.05)

        the_mean_rs = np.mean(the_rs)
        the_means_dif = np.abs(the_nt_mean - the_mean_rs)
        the_stat_error = the_sigmas_corr[tt]

        vfp.PLOT_HISTOGRAMS(the_rs, r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error), the_mean_rs, r'$ \bar{C}_{%s}(t) = $'%the_type_rs + r' $%s$'%'{:.10e}'.format(the_mean_rs), the_nt_mean, r'$ \bar{C}(t) = %s$'%'{:.10e}'.format(the_nt_mean), OperatorNamePlot + rf' $\to$ $t = {tt+the_nt[0]} a$', the_nr_bins,  r'Correlator')
        plt.show()
        the_gauss_fig.savefig(f'{the_location}Histogram_correlators_{irrep[:4]}_{irrep[-1]}{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        

def PlotMultiHadronCorrelators(the_matrix_correlator_data, the_quantum_number, the_type_rs, the_version, the_t0, the_location, the_rebin, **kwargs):
    
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
        
    ### Choosing the resampling scheme to put in the legend
    if the_type_rs=='jk':
        the_rs_scheme='Jackknife'
    elif the_type_rs=='bt':
        the_rs_scheme='Bootstrap'

    ### This is the nr. of bins to plot the histograms
    # the_nr_bins = 25
    
    ### Loop over the irreducible representations
    for irrep in m_irreps:
        ### The list of operators of this irrep
        the_op_list = list(the_matrix_correlator_data[f'{irrep}/Operators'])
        
        ### The correaltor dataset
        the_data_corr = vfa.RESHAPING_CORRELATORS(np.asarray(the_matrix_correlator_data[f'{irrep}/Correlators/Real/Mean']))
        
        ### The sigmas of this correlator dataset
        the_data_sigmas_corr = np.asarray(the_matrix_correlator_data[f'{irrep}/Correlators/Real/Sigmas'])
        
        ### This is the time interval
        the_nt = np.asarray(the_matrix_correlator_data[f'{irrep}/Time_slices'])

        ### These are going to be the ticks in the x-label of the plots
        the_nt_ticks = np.arange(the_nt[0]+1, the_nt[-1], int(len(the_nt)/5))
        
        ### Plotting eigenvalues in the full time range or only starting from t0
        if kwargs.get('full_range_nt')==None: the_start_nt = the_t0-the_nt[0]
        else: the_start_nt = the_nt[0]
        
        ### Information of the irrep to write it properly in the plots.
        da_irrep = vfp.IrrepInfo(irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot 
        NameIrrep = da_irrep.Name
        
        if the_diagonal_corrs_flag:
            
            ### Loop over the number of operators of this matrix
            for bb in range(len(the_op_list)):      
                
                ### The specific operator
                the_op = the_op_list[bb]
                theOperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
                OperatorNamePlot = vfp.MH_OPERATORS_RELABEL(theOperatorNamePlot)
                print('Correlator plots in process...')
                    
                ### Plotting the diagonal correlators
                corr_fig = plt.figure()
                vfp.PLOT_CORRELATORS(the_nt, the_data_corr[bb,bb], the_data_sigmas_corr[bb], the_rs_scheme, the_nt_ticks, r'$t$', r'$\mathbb{Re}\;C(t)$', 'o', f'{NameIrrepPlot} ({MomentumIrrep}) ' + r' $\to \;C_{%s}$'%(str(bb)+str(bb)) + f'= {OperatorNamePlot}')
                plt.show()
                corr_fig.savefig(f'{the_location}DiagonalCorrelator_{the_quantum_number}_{irrep}_{bb}{the_rebin}_{the_version}.pdf')
                
                ## Plotting the log of the diagonal correlators.
                print('Correlator Log-plots in progress...')
                corr_fig = plt.figure()
                vfp.PLOT_CORRELATORS(the_nt, the_data_corr[bb,bb], the_data_sigmas_corr[bb], the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$\log\mathbb{Re}\;C(t)$', 'o', f'{NameIrrepPlot} ({MomentumIrrep}) ' + r' $\to \;C_{%s}$'%(str(bb)+str(bb)) + f'= {OperatorNamePlot}', yscale='log')
                plt.show()
                corr_fig.savefig(f'{the_location}DiagonalCorrelator_{the_quantum_number}_{irrep}_{bb}_log{the_rebin}_{the_version}.pdf')
                
                ### Plotting the histogram at a certain time slice t
                print('Correlator histogram in progress...')
                tt = int(len(the_nt)/3)+1
                the_gauss_fig = plt.figure()
                the_nt_mean = the_data_corr[bb,bb,tt]
                the_rs = vfa.RESHAPING_CORRELATORS_RS_NT(np.asarray(the_matrix_correlator_data[f'{irrep}/Correlators/Real/Resampled']))[bb,bb,tt]
                the_nr_bins = int((np.asarray(the_matrix_correlator_data[f'{irrep}/Correlators/Real/Resampled']).shape[1])*.05)
                
                ### Here the mean value of the sampling data and the mean value are compared to check the quality of the resampled data
                the_mean_rs = np.mean(the_rs)
                the_means_dif = np.abs(the_nt_mean - the_mean_rs)
                the_stat_error = the_data_sigmas_corr[bb,tt]
                
                ### Plotting the histogram now
                vfp.PLOT_HISTOGRAMS(the_rs, r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error), the_mean_rs, r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%'{:.10e}'.format(the_mean_rs), the_nt_mean, r'$ \bar{C}(t) = $ %s'%'{:.10e}'.format(the_nt_mean), f'{NameIrrepPlot} ({MomentumIrrep}) ' +  r'$\to \;C_{%s}$'%(str(bb) + str(bb)) + r' $(t = %sa) = $'%(tt+the_nt[0]) + f' {OperatorNamePlot}', the_nr_bins, r'$C(t)$')
                plt.show()
                the_gauss_fig.savefig(f'{the_location}Histogram_DiagCorrelator_{the_quantum_number}_{irrep}_{bb}{the_rebin}_{the_version}.pdf')
            
            ### The Diagonal of the correlators are plotted all together with their errors to compare them directly. 
            if all_corr_flag:
                corr_fig = plt.figure()
                print('ALL Correlators Log-plot in progress...')
                ### Loop over each of the entries of the diagonal of the correlation matrix
                for bb in range(len(the_op_list)):
                    
                    ### Name of this operator
                    the_op = the_op_list[bb]
                    theOperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                    OperatorNamePlot = vf.MH_OPERATORS_RELABEL(theOperatorNamePlot)           
                    
                    # plt.errorbar(the_nt, the_data_corr[bb][bb], the_data_sigmas_corr[bb],  marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=3.5, label = r'$C_{%s}$ = '%(str(bb)+str(bb)) + OperatorNamePlot, color=the_colors[bb])
                    plt.errorbar(the_nt, the_data_corr[bb][bb], the_data_sigmas_corr[bb],  marker=the_markers_list[bb], ls='None', ms=6, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, capsize=5, label = OperatorNamePlot, color=the_colors[bb])
                plt.xlabel(r'$t\,/\,a$',fontsize=36)
                plt.ylabel(r'$\log\mathbb{Re}\;C(t)$', fontsize=36)
                # plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;C_{ii}(t)$',fontsize=20)
                plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep ,fontsize=26)
                # the_n_cols = int(len(the_op_list)/2)
                plt.legend(fontsize=16, ncol= 2, columnspacing=0.1, handletextpad=0.01)
                plt.yscale('log')
                plt.tight_layout()
                # if len(the_op_list)>7: the_n_cols = int(len(the_op_list)/3)
                # else: the_n_cols = int(len(the_op_list)/2)
                plt.xticks(the_nt_ticks,fontsize=18)
                plt.yticks(fontsize=18)
                plt.show()
                corr_fig.savefig(the_location + 'ALLDiagonalCorrelators' + the_quantum_number + irrep  + '_log' + the_rebin + '_v%s.pdf'%the_version)
            
#             if the_ops_analysis_flag:
#                 ### The Diagonal of the correlators are plotted all together with their errors to compare them directly. 
#                 corr_fig = plt.figure()
#                 print('Chosen Correlators Log-plot in progress...')
#                 
#                 the_colors_index = 0
#                 the_temporary_op_list = the_chosen_op_list[m_irreps.index(irrep)]
#                 the_label_string= ''
#                 
#                 ### Loop over each of the entries of the diagonal of the correlation matrix
#                 for bb in the_temporary_op_list:
#                     
#                     ### Name of this operator
#                     the_op = the_op_list[bb]
#                     theOperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
#                     OperatorNamePlot = vf.MH_OPERATORS_RELABEL(theOperatorNamePlot)
#                     the_label_string+=str(bb)
#                     
#                     plt.errorbar(the_nt, the_data_corr[bb][bb], the_data_sigmas_corr[bb],  marker=the_markers_list[bb], ls='None', ms=6, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=5, label = r'$C_{%s}$ = '%(str(bb)+str(bb)) + OperatorNamePlot, color=the_colors[the_colors_index])
#                     the_colors_index+=1
#                 plt.xlabel(r'$t\,/\,a$',fontsize=36)
#                 plt.ylabel(r'$\log\mathbb{Re}\;C(t)$', fontsize=36)
#                 # plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;C_{ii}$(t)',fontsize=20)
#                 plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep,fontsize=26)
#                 plt.yscale('log')
#                 plt.legend(fontsize=16, ncol= 2, handletextpad=0.01, columnspacing=0.1)
#                 plt.tight_layout()
#                 # if len(the_temporary_op_list)>7: the_n_cols = int(len(the_temporary_op_list)/3)
#                 # else: the_n_cols = int(len(the_temporary_op_list)/2)
#                 # plt.legend(fontsize=16, ncol= the_n_cols, handletextpad=0.01)
#                 plt.xticks(the_nt_ticks)
#                 # plt.show()
#                 corr_fig.savefig(the_location + 'Chosen%s_DiagonalCorrelators'%the_label_string + the_quantum_number + irrep  + '_log' + the_rebin + '_v%s.pdf'%the_version)
                
            
        ### Here the Eigenvalues are plotted all together too in a log-plot
        if 'GEVP' in list(the_matrix_correlator_data[irrep].keys()):
            
            the_data = np.asarray(the_matrix_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Eigenvalues/Mean'])
            the_data_sigmas = np.asarray(the_matrix_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Eigenvalues/Covariance_matrix'])
            
            for bb in range(len(the_data)):
                ### The mean value of the eigenvalue_{i}
                the_mean_corr = the_data[bb]
                
                ### The corresponding sigmas of this eigenvalue
                the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                
                print("..................................................")
                print(f'Eigenvalue = {bb} plot in process...')
                
                corr_fig = plt.figure()
                ### Plotting the eigenvalues one by one
                vfp.PLOT_CORRELATORS(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], the_sigmas_corr[the_start_nt:], the_rs_scheme, the_nt_ticks, r'$t$', r'$\lambda_{i}(t)$', 'o',  f'{NameIrrepPlot} ({MomentumIrrep}): ' + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %s$)'%str(the_t0))
                plt.show()
                corr_fig.savefig(f'{the_location}Eigenvalues_{the_quantum_number}_{irrep}_{bb}_t0_{the_t0}{the_rebin}_{the_version}.pdf')
                
                ### Plotting the eigenvalues log-plots one by one
                print('Eigenvalue = %s Log-plot in progress...'%str(bb))
                corr_fig = plt.figure()
                vfp.PLOT_CORRELATORS(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], the_sigmas_corr[the_start_nt:], the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$\log\,(\lambda_{i}(t))$', 'o',  f'{NameIrrepPlot} ({MomentumIrrep}) ' + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %sa$)'%str(the_t0), yscale='log')
                plt.show()
                corr_fig.savefig(f'{the_location}Eigenvalues_{the_quantum_number}_{irrep}_{bb}_log_t0_{the_t0}{the_rebin}_{the_version}.pdf')
                
                ### Plotting the histograms of these eigenvalues
                print(f'Eigenvalue = {bb} histogram in progress...')
                tt = int(len(the_nt)/3)+1
                the_gauss_fig = plt.figure()
                the_nt_mean = the_mean_corr[tt]
                the_rs = np.asarray(the_matrix_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Eigenvalues/Resampled'])[bb].transpose()[tt]
                the_nr_bins = int(np.asarray(the_matrix_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Eigenvalues/Resampled']).shape[1] * .05)

                the_mean_rs = np.mean(the_rs)
                the_means_dif = np.abs(the_nt_mean - the_mean_rs)
                the_stat_error = the_sigmas_corr[tt]
                
                vfp.PLOT_HISTOGRAMS(the_rs, r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error), the_mean_rs, r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%'{:.10e}'.format(the_mean_rs), the_nt_mean, r'$ \bar{C}(t) = $ %s'%'{:.10e}'.format(the_nt_mean), f'{NameIrrepPlot} ({MomentumIrrep}) ' +  r'$\to\; \lambda_{%s}$'%str(bb) + r'$(t_{0} = %sa)$'%str(the_t0) +  r' $[t = %sa]$'%(tt+the_nt[0]), the_nr_bins, r'Eigenvalue ($\lambda_{%s}$)'%str(bb))
                plt.show()
                the_gauss_fig.savefig(f'{the_location}Histogram_Eigenvalues_{the_quantum_number}_{irrep}_{bb}_t0_{the_t0}{the_rebin}_{the_version}.pdf')
            
            if all_corr_flag:
                print("..................................................\n")
                
                print('ALL Eigenvalues Log-plot in progress...')
                corr_fig = plt.figure()           
                for bb in range(len(the_data)):
                    the_mean_corr = the_data[bb]
                    the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                    
                    plt.errorbar(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], yerr = the_sigmas_corr[the_start_nt:], marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=3.5, label = r'$\lambda_{%s}$'%str(bb), color=the_colors[bb])
                plt.xlabel(r'$t\,/\,a$',fontsize=28)
                plt.ylabel(r'$\log\,(\lambda_{i}(t))$',fontsize=28)
                plt.title( f'{NameIrrepPlot} ({MomentumIrrep}) ' + r'$\to\;\lambda_{i} (t_{0} = %sa$)'%str(the_t0), fontsize=18)
                plt.yscale('log')
                plt.legend(fontsize=16, ncol= 2, handletextpad=0.01)
                plt.tight_layout()
                # if len(the_data)>6: the_n_cols = int(len(the_data)/3)
                # else: the_n_cols = int(len(the_data)/2)
                # plt.legend(fontsize=12, ncol=the_n_cols, handletextpad=0.3)
                plt.xticks(the_nt_ticks)
                plt.show()
                corr_fig.savefig(f'{the_location}ALLEigenvalues_{the_quantum_number}_{irrep}_log_t0_{the_t0}{the_rebin}_{the_version}.pdf')



def PlotRatioHadronCorrelators(the_ratio_correlator_data, the_quantum_number, the_type_rs, the_version, the_t0, the_location):
    mr_irreps = list(the_ratio_correlator_data.keys())
    the_nr_bins = 25
    for irrep in mr_irreps:
            the_op_list = list(the_ratio_correlator_data[irrep+'/Operators'])
            the_data = np.array(the_ratio_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Mean'%the_t0])
            the_data_sigmas = np.array(the_ratio_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Covariance_matrix'%the_t0])
            for bb in range(len(the_op_list)):
                the_mean_corr = the_data[bb]
                the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                the_nt = np.array(the_ratio_correlator_data[irrep + '/Time_slices'])
                the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 5)

                the_op = the_op_list[bb]
                OperatorNamePlot = vf.OPERATORS_SH(the_op.decode('utf-8'))
                da_irrep = vf.IrrepInfo(irrep)
                
                MomentumIrrep = da_irrep.TotalMomPlot
                NameIrrepPlot = da_irrep.NamePlot
                NameIrrep = da_irrep.Name
                
                print('Correlator plot in process...')
                the_corr_fig = plt.figure()
                plt.errorbar(the_nt[the_t0-the_nt[0]:], the_mean_corr[the_t0-the_nt[0]:], yerr = the_sigmas_corr[the_t0-the_nt[0]:], marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
                plt.xlabel(r'$t$')
                plt.ylabel(r'$\mathbb{Re}\;C(t)$')
                plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0, fontsize=18)
                plt.tight_layout()
                plt.xticks(the_nt_ticks)
                #plt.show()
                the_corr_fig.savefig(the_location + 'Ratios/Eigenvalues_ratios' + the_quantum_number + irrep  + '_%s_'%bb + the_rebin + 'v%s.pdf'%the_version)
                
                print('Correlator Log-plot in process...')
                the_log_corr_fig = plt.figure()
                plt.errorbar(the_nt[the_t0-the_nt[0]:], the_mean_corr[the_t0-the_nt[0]:], yerr = the_sigmas_corr[the_t0-the_nt[0]:], marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
                plt.xlabel(r'$t$')
                plt.ylabel(r'$\log\mathbb{Re}\;C(t)$')
                plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0, fontsize=18)
                plt.yscale('log')
                plt.tight_layout()
                plt.xticks(the_nt_ticks)
                #plt.show()
                the_log_corr_fig.savefig(the_location + 'Ratios/Eigenvalues_ratios' + the_quantum_number + irrep + '_%s_log'%bb + the_rebin + '_v%s.pdf'%the_version)
                
                print('Correlator histogram in progress...')
                tt = int(len(the_nt)/3)+1
                the_gauss_fig = plt.figure()
                the_mean = the_mean_corr[tt]
                the_rs = np.array(mr[irrep + '/GEVP/t0_%s/Eigenvalues/Resampled'%the_t0])[bb].transpose()[tt]

                the_mean_rs = np.mean(the_rs)
                the_means_dif = np.abs(the_mean - the_mean_rs)
                the_stat_error = the_sigmas_corr[tt]
                the_nr_samples = np.array(mr[irrep + '/GEVP/t0_%s/Eigenvalues/Resampled'%the_t0])[bb].shape[0]
                
                plt.hist(the_rs, bins=the_nr_bins, label =  r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error))
                plt.vlines(the_mean_rs, 0,  the_nr_samples*(the_nr_bins/100)*.65, colors= 'red', label = r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%the_mean_rs)
                plt.vlines(the_mean, 0,  the_nr_samples*(the_nr_bins/100)*.65, colors='black', label = r'$ \bar{C}(t) = $ %s'%the_mean)
                plt.title( NameIrrep + ' (%s): '%MomentumIrrep +  r'$\lambda_{%s}$'%bb +  r' $t =$ %s'%(tt+the_nt[0]) + r' ($t_{0} = %s$)'%the_t0, fontsize=18)
                plt.ylabel('Frequency')
                plt.xlabel(r'Eigenvalue ($\lambda_{%s}$)'%bb)
                plt.legend()
                plt.tight_layout()
                plt.ylim([0,the_nr_samples*(the_nr_bins/100)*.6])
                # plt.show()
                the_gauss_fig.savefig(the_location + 'Histogram_Eigenvalues_ratios' + the_quantum_number + irrep + '_%s_'%bb + the_rebin + 'v%s.pdf'%the_version)



### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------




if __name__=="__main__":
    import ensembles as ed
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
        rb = int(myRb)
        reBin = '_bin'+str(rb)
    else:
        reBin = '' 
        
    if myTypeRs=='jk':
        myResamplingScheme='Jackknife'
    elif myTypeRs=='bt':
        myResamplingScheme='Bootstrap'    

    if myWhichCorrelator=='s':
        mySingleCorrelatorData = h5py.File(myDataLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/Plots/%s/SingleHadrons/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotSingleHadronCorrelators(mySingleCorrelatorData, myTypeRs, myVersion, myPlotLocation, reBin)
        
        mySingleCorrelatorData.close()
    
    elif myWhichCorrelator=='m':
        myMatrixCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/Plots/%s/Matrices/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotMultiHadronCorrelators(myMatrixCorrelatorData, myTypeRs, myVersion, myT0, myPlotLocation, reBin)
        
        myMatrixCorrelatorData.close()
        
    elif myWhichCorrelator=='mr':
        myRatioCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/Plots/%s/Matrices_Ratios/'%myEns +  '%s/'%myResamplingScheme)
        
        myRatioCorrelatorData.close()
    
