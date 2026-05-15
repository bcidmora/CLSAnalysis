###     PLOTTING SCRIPTS
import plot_correlators_script as pcorr
import plot_effective_masses_script as peff
import plot_fits_script as pfit
# import set_of_functions as vf
import plot_fitted_eff_masses as pfem

### LAYOUT AND EXTRA FUNCTIONS
import set_of_layout_functions as vfl
import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp

### EXTRA LIBRARIES
import sys
import os
import h5py
from PyPDF2 import PdfMerger

import ensembles as ed


### ------- WHAT IT IS DONE --------

myArgs = vfp.parse_args()
myRuns = vfp.WhichRuns(myArgs, ed.ensembles)


### ------ MAIN VARIABLES ---------
reBin = f"_bin{myRuns.rb}" if myRuns.rebin else ""

### This is the information about the resampling
if myRuns.rs_type=='jk':
    myResamplingScheme='Jackknife'
elif myRuns.rs_type=='bt':
    myResamplingScheme='Bootstrap' 
    

### Information from the ensembles dictionary
myDataLocation = vfl.DIRECTORY_EXISTS(f'{ed.outputLocation}{myRuns.ensemble}/')
myPlotLocation = vfl.DIRECTORY_EXISTS(f'{ed.location}/Plots/{myRuns.ensemble}/')

### -------- PRINTING INFO OF ENSEMBLE ---------

vfl.INFO_PRINTING(myRuns.correlator, myRuns.ensemble)

### ------------ START ----------------

if myRuns.correlator =='s':
    myVersion =  f'{myRuns.ensemble}_singles_test' 
    mySingleTmins = ed.ensembles[myRuns.ensemble]['singleTMinResults']
    
    ### Original list of irreps
    myArchivo = h5py.File(ed.ensembles[myRuns.ensemble]['fs'], 'r')
    myIrreps = list(myArchivo.keys())
    myArchivo.close()

    ### Correlators data
    mySingleCorrelatorData = h5py.File(f'{myDataLocation}/Single_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5','r')
    
    ### Directory where the plots will be saved
    myPlotLocation += f'SingleHadrons/{myResamplingScheme}/'
    
    if myRuns.corrs:    
        pcorr.PlotSingleHadronCorrelators(mySingleCorrelatorData, myRuns.rs_type, myVersion, myPlotLocation, reBin, nr_irreps = myRuns.the_irreps.nr_irreps, first_irrep = myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep)
    
    if myRuns.effmass: 
        peff.PlotSingleHadronsEffectiveMasses(mySingleCorrelatorData, myResamplingScheme, myVersion, myPlotLocation, reBin, nr_irreps=myRuns.the_irreps.nr_irreps, first_irrep=myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep)
    
    if myRuns.fits: 
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myDataLocation}/Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(f'{myFitsLocation}Single_correlators_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
        
        pfit.PlotSingleHadronsFits(myFitCorrelator, myRuns.fit_param.type_correlation, myRuns.fit_param.type_fit, mySingleTmins, myVersion, myPlotLocation, reBin, myIrreps, first_irrep = myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep, zoom_fit = myRuns.zoom_fit, chi_plots = myRuns.plot_chi, total_chi = myRuns.total_chi_plot, delta_chi = myRuns.delta_chi_plot, all_fits_comparison = myRuns.all_fits, nr_irreps=myRuns.the_irreps.nr_irreps)
        
        myFitCorrelator.close()
    
    if myRuns.fitmass:        
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myDataLocation}Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(f'{myFitsLocation}Single_correlators_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
        
        pfem.PlotSingleHadronsEffectiveMassesFits(myFitCorrelator, mySingleCorrelatorData, myResamplingScheme, myRuns.fit_param.type_correlation, myRuns.fit_type, mySingleTmins, myVersion, myPlotLocation, reBin, myIrreps, first_irrep=myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep)
    
    ### Puts all the plots in one PDF file. It checks if the file exists first
    if myRuns.join:
        ### Loop over all the irreps in this ensemble
        irreps = list(mySingleCorrelatorData.keys())
        all_ef_mass_x = []
        for aa in irreps:
            ### Corrs
            the_corr_plot = f'{myPlotLocation}Correlator_{aa[:4]}_{aa[-1]}{reBin}_{myVersion}.pdf'
            
            ### Corrs log-plots
            the_corr_log_plot = f'{myPlotLocation}Correlator_{aa[:4]}_{aa[-1]}_log{reBin}_{myVersion}.pdf'
            
            ### Histogram Corrs
            the_hist_plot = f'{myPlotLocation}Histogram_correlators_{aa[:4]}_{aa[-1]}{reBin}_{myVersion}.pdf'
            
            ### Effective Masses Corrs
            the_eff_mass_plot = f'{myPlotLocation}EffectiveMass_{aa[:4]}_{aa[-1]}{reBin}_{myVersion}.pdf'
            
            ### Fits Corrs
            the_tmin_plot = f'{myPlotLocation}Tmin_Fits_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
            
            ### Zoom Fits Corrs
            the_tmin_zoom_plot = f'{myPlotLocation}Tmin_Fits_Zoom_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
            
            ### Chi^{2} Fits Corrs
            the_chi_plot = f'{myPlotLocation}Tmin_Chisqr_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
            
            ### Total Chi^{2} Fits Corrs
            the_total_chi_plot = f'{myPlotLocation}Tmin_TotalChisqr_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
            
            ### Delta Chi^{2} Fits Corrs
            the_delta_chi_plot = f'{myPlotLocation}Tmin_DeltaChisqr_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
            
            ### All fits in one plot
            the_diff_fits_plot = f'{myPlotLocation}Different_Fits_{aa[:4]}_{aa[-1]}_{reBin}_{myVersion}.pdf'
            
            ### Fitted effective masses
            the_fitted_mass_plot = f'{myPlotLocation}Fitted_Effective_Masses_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp_{reBin}_{myVersion}.pdf'
            
            the_plots = [the_corr_plot, the_corr_log_plot, the_hist_plot, the_eff_mass_plot, the_tmin_plot, the_tmin_zoom_plot, the_chi_plot, the_total_chi_plot, the_delta_chi_plot, the_diff_fits_plot, the_fitted_mass_plot]
        
            x=[]
            for item in the_plots:
                if os.path.isfile(item): x.append(item)
                if 'Fitted_Effective_Masses' in item: all_ef_mass_x.append(item)
                
            merger = PdfMerger()
            for pdf in x:
                merger.append(open(pdf, 'rb'))
            
            with open(f'{myPlotLocation}{myRuns.ensemble}_{aa}{reBin}_{myRuns.rs_type}_{myVersion}.pdf', "wb") as fout:
                merger.write(fout)
        
        merger_eff = PdfMerger()
        for pdf_eff in all_ef_mass_x:
            merger_eff.append(open(pdf_eff, 'rb'))
        with open(f'{myPlotLocation}{myRuns.ensemble}_All_Fitted_Masses{reBin}_{myRuns.rs_type}_{myVersion}.pdf', "wb") as fout:
            merger_eff.write(fout)
        print('Now all the plots are in one file for each irrep')
    
    mySingleCorrelatorData.close()
    
elif myRuns.correlator=='m':
    myIsospin = myRuns.isospin
    myChosenIsospin = ed.ensembles[myRuns.ensemble][myIsospin]['iso_tag']
    myArchivo = h5py.File(ed.ensembles[myRuns.ensemble][myIsospin]['fm'], 'r')
    myIrreps = list(myArchivo.keys())
    myArchivo.close()
    
    if myRuns.ops_flag:
        myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_reduced_test' 
    else:
        myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_test' 
    
    myMatrixCorrelatorData = h5py.File(f'{myDataLocation}Matrix_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5','r')
    
    ### Directory where the plots will be saved
    myPlotLocation += f'Matrices/{myResamplingScheme}/'
    
     ### Plots the correlators. Look at the booleans
    if myRuns.corrs: 
        pcorr.PlotMultiHadronCorrelators(myMatrixCorrelatorData, myChosenIsospin, myRuns.rs_type, myVersion, myRuns.fit_param.t0, myPlotLocation, reBin, nr_irreps = myRuns.the_irreps.nr_irreps, first_irrep = myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep, diag_corrs = myRuns.diag_flag)
    
    ### Plots the effective masses of the eigenvalues from the GEVP and/or the operators analysis
    if myRuns.effmass: 
        peff.PlotMultiHadronsEffectiveMasses(myMatrixCorrelatorData, myChosenIsospin, myResamplingScheme, myVersion, myRuns.fit_param.t0, myPlotLocation, reBin, nr_irreps=myRuns.the_irreps.nr_irreps, first_irrep=myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep, diag_corrs= myRuns.diag_flag)

        
    if myRuns.fits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_Matrices/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators' + myChosenIsospin + myRuns.rs_type + reBin + '_fits_v%s.h5'%myVersion, 'a')
                             
        pfit.PlotMultiHadronsFits(myFitCorrelator, myChosenIsospin, myRuns.fit_param.type_correlation, myRuns.fit_param.type_fit, myRuns.rs_type, multiTMinsFitPlots, myRuns.fit_param.t0, myVersion, myPlotLocation, reBin, myIrreps, gevp=myRuns.gevp_flag, zoom_fit=myRuns.zoom_fit, chi_plots=myRuns.plot_chi, first_irrep=myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep, total_chi=myRuns.total_chi_plot, delta_chi=myRuns.delta_chi_plot, ops_analysis=myRuns.ops_flag, ops_analysis_method=myOperatorsMethod)
        
        myFitCorrelator.close()
    
    ### Plot the fitted effective masses of the eigenvalues
    if myRuns.fitmass:
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_Matrices/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators' + myChosenIsospin + myRuns.rs_type + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        pfem.PlotMultiHadronsEffectiveMassesFits(myFitCorrelator, myMatrixCorrelatorData, myChosenIsospin, myResamplingScheme, myRuns.fit_param.type_correlation, myRuns.fit_param.type_fit, multiTMinsFitPlots, myRuns.fit_param.t0, myVersion, myPlotLocation, reBin, myIrreps) 
        myFitCorrelator.close()

    ### This part puts all the plots together in one pdf file. 
    if myRuns.join:
        ### Loop over all the irreps in this ensemble
        irreps = list(myMatrixCorrelatorData.keys())
        for aa in irreps:
            ##3 Loop over all the operators in this ensemble
            ops = list(myMatrixCorrelatorData[aa+'/Operators'])
            x=[]
            for bb in range(len(ops)):
                
                ### Diagonal Corrs
                # if os.path.isfile(myPlotLocation + 'DiagonalCorrelator' + myChosenIsospin + aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion):
                    # x.append(myPlotLocation + 'DiagonalCorrelator' + myChosenIsospin + aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion)
                
                ### Diagonal corrs log-plot
                # if os.path.isfile(myPlotLocation + 'DiagonalCorrelator' + myChosenIsospin + aa + '_%s_log'%str(bb) + reBin + '_v%s.pdf'%myVersion):
                #     x.append(myPlotLocation + 'DiagonalCorrelator' + myChosenIsospin + aa + '_%s_log'%str(bb) + reBin + '_v%s.pdf'%myVersion)
                
                # ### Histogram diagonal corrs
                # if os.path.isfile(myPlotLocation + 'Histogram_DiagCorrelator' + myChosenIsospin + aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion):
                #     x.append(myPlotLocation + 'Histogram_DiagCorrelator' + myChosenIsospin + aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion)
               
                ### Effective Mass diagonal corrs
                if os.path.isfile(myPlotLocation + 'EffectiveMass_DiagonalCorrelators' + myChosenIsospin + aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'EffectiveMass_DiagonalCorrelators' + myChosenIsospin + aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion)
                
                # ### Eigenvalues 
                # if os.path.isfile(myPlotLocation + 'Eigenvalues' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion):
                #     x.append(myPlotLocation + 'Eigenvalues' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
                
                ### Eigenvalues log-plots
#                 if os.path.isfile(myPlotLocation + 'Eigenvalues' + myChosenIsospin + aa +  '_%s_log'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion):
#                     x.append(myPlotLocation + 'Eigenvalues' + myChosenIsospin + aa +  '_%s_log'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
#                 
                # ### Histogram Eigenvlaues
                # if os.path.isfile(myPlotLocation +'Histogram_Eigenvalues' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin +  '_v%s.pdf'%myVersion):
                #     x.append(myPlotLocation +'Histogram_Eigenvalues' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin +  '_v%s.pdf'%myVersion)
                
                # ### Effective Mass Eigenvlaues
                # if os.path.isfile(myPlotLocation + 'EffectiveMass_Eigenvalues' + myChosenIsospin+ aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + reBin + '_v%s.pdf'%myVersion):
                #     x.append(myPlotLocation + 'EffectiveMass_Eigenvalues' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
                
#                 ### Fits Eigenvalues
#                 if os.path.isfile(myPlotLocation + 'Tmin_Fits' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + '_%sexp'%myRuns.fit_param.type_fit + reBin + '_v%s.pdf'%myVersion):
#                     x.append(myPlotLocation + 'Tmin_Fits' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + '_%sexp'%myRuns.fit_param.type_fit + reBin + '_v%s.pdf'%myVersion)
#                 
#                 ### Chi^{2} Fits Eigenvalues 
#                 if os.path.isfile(myPlotLocation + 'Tmin_Chisqr' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + '_%sexp'%myRuns.fit_param.type_fit + reBin + '_v%s.pdf'%myVersion):
#                     x.append(myPlotLocation + 'Tmin_Chisqr' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + '_%sexp'%myRuns.fit_param.type_fit + reBin + '_v%s.pdf'%myVersion)
#                 
#                 ### Zoom Chi^{2} Fits Eigenvalues 
#                 if os.path.isfile(myPlotLocation + 'Tmin_Chisqr_Zoom' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + '_%sexp'%myRuns.fit_param.type_fit + reBin + '_v%s.pdf'%myVersion):
#                     x.append(myPlotLocation + 'Tmin_Chisqr_Zoom' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + '_%sexp'%myRuns.fit_param.type_fit + reBin + '_v%s.pdf'%myVersion)
#                     
#                  ### Total Chi^{2} Fits Eigenvalues 
#                 if os.path.isfile(myPlotLocation + 'Tmin_TotalChisqr' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + '_%sexp'%myRuns.fit_param.type_fit + reBin + '_v%s.pdf'%myVersion):
#                     x.append(myPlotLocation + 'Tmin_TotalChisqr' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + '_%sexp'%myRuns.fit_param.type_fit + reBin + '_v%s.pdf'%myVersion)
#                     
#                 ### Delta Chi^{2} Fits
#                 if os.path.isfile(myPlotLocation + 'Tmin_DeltaChisqr' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + '_%sexp'%myRuns.fit_param.type_fit + reBin + '_v%s.pdf'%myVersion):
#                     x.append(myPlotLocation + 'Tmin_DeltaChisqr' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + '_%sexp'%myRuns.fit_param.type_fit + reBin + '_v%s.pdf'%myVersion)
#                     
#                 ### Fits Effetive Masses of Eigenvalues
#                 if os.path.isfile(myPlotLocation + 'Fitted_Effective_Masses' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + reBin + '_v%s.pdf'%myVersion):
#                     x.append(myPlotLocation + 'Fitted_Effective_Masses' + myChosenIsospin + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0)  + reBin + '_v%s.pdf'%myVersion)
#                 
            ### All Correlators together log-plot
            if os.path.isfile(myPlotLocation + 'ALLDiagonalCorrelators' + myChosenIsospin + aa + '_log'+ reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'ALLDiagonalCorrelators' + myChosenIsospin + aa + '_log'+ reBin + '_v%s.pdf'%myVersion)

            ### Effective Masses all diagonal correlators together
            if os.path.isfile(myPlotLocation + 'EffectiveMass_ALLDiagonalCorrelators' + myChosenIsospin + aa + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'EffectiveMass_ALLDiagonalCorrelators'+ myChosenIsospin + aa + reBin + '_v%s.pdf'%myVersion)
            
            # ### All eigenvalues together log-plots
            # if os.path.isfile(myPlotLocation + 'ALLEigenvalues' + myChosenIsospin + aa + '_log'+ '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion):
            #     x.append(myPlotLocation + 'ALLEigenvalues'+ myChosenIsospin + aa + '_log'+ '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
        
            ### Effective Masses all eigenvalues together
            if os.path.isfile(myPlotLocation + 'EffectiveMass_ALLEigenvalues' + myChosenIsospin + aa + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'EffectiveMass_ALLEigenvalues'+ myChosenIsospin + aa + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
            
#             ### If the operators analysis was performed, then it will look for the plots
#             if 'Operators_Analysis' in list(myMatrixCorrelatorData[aa].keys()):
#                 
#                 ### What type of operators analysis is in
#                 if any('Ops_chosen' in the_keys for the_keys in list(myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys())):
#                     
#                     the_list_of_chosen_ops = list(filter(lambda x: 'Ops_chosen' in x, myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys()))
#                     
#                     ### Loop over all the operators analysis
#                     for the_op_item in the_list_of_chosen_ops:
#                         the_len_data = len(myMatrixCorrelatorData[aa+'/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Mean'%myRuns.fit_param.t0])
#                         ### Loop over the eigenvalues of this analysis
#                         for bb in range(the_len_data):
#                             if os.path.isfile(myPlotLocation + 'EffectiveMass_Eigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion):
#                                 x.append(myPlotLocation + 'EffectiveMass_Eigenvalues'+ myChosenIsospin + the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
#                     
#                     ### All eigenvalues together in one plot
#                     if os.path.isfile(myPlotLocation + 'EffectiveMass_ALLEigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion):
#                             x.append(myPlotLocation + 'EffectiveMass_ALLEigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
#                 
#                 if any('Add_Op' in the_keys for the_keys in list(myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys())):
#                     
#                     
#                     the_list_of_chosen_ops = list(filter(lambda x: 'Add_Op' in x, myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys()))
#                     
#                     ### Loop over all the operators analysis
#                     for the_op_item in the_list_of_chosen_ops:
#                         the_len_data = len(myMatrixCorrelatorData[aa+'/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Mean'%myRuns.fit_param.t0])
#                         
#                         ### Loop over the eignvalues
#                         for bb in range(the_len_data):
#                             if os.path.isfile(myPlotLocation + 'EffectiveMass_Eigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion):
#                                 x.append(myPlotLocation + 'EffectiveMass_Eigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
#                     
#                     ### All the eigenvalues together in one plot
#                     if os.path.isfile(myPlotLocation + 'EffectiveMass_ALLEigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion):
#                             x.append(myPlotLocation + 'EffectiveMass_ALLEigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
#                         
#                 if any('Remove_Op' in the_keys for the_keys in list(myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys())):
#                     the_list_of_chosen_ops = list(filter(lambda x: 'Remove_Op' in x, myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys()))
#                     
#                     ### Loop over all the operators of the analysis
#                     for the_op_item in the_list_of_chosen_ops:
#                         the_len_data = len(myMatrixCorrelatorData[aa+'/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Mean'%myRuns.fit_param.t0])
#                         
#                         ### Loop over the eigenvalues of this selection
#                         for bb in range(the_len_data):
#                             if os.path.isfile(myPlotLocation + 'EffectiveMass_Eigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion):
#                                 x.append(myPlotLocation + 'EffectiveMass_Eigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
#                     
#                     ### All eigenvalues together in one plot
#                     if os.path.isfile(myPlotLocation + 'EffectiveMass_ALLEigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion):
#                         x.append(myPlotLocation + 'EffectiveMass_ALLEigenvalues' + myChosenIsospin + the_op_item +'_' + aa + '_t0_%s'%str(myRuns.fit_param.t0) + reBin + '_v%s.pdf'%myVersion)
                        
            merger = PdfMerger()
            for pdf in x:
                merger.append(open(pdf, 'rb'))
            
            with open(myPlotLocation + myRuns.ensemble + myChosenIsospin + "_%s"%aa + "_t0%s"%str(myRuns.fit_param.t0) + reBin + "_%s"%myRuns.rs_type + "_v%s.pdf"%myVersion, "wb") as fout:
                merger.write(fout)


        print('Now all the plots are in one file for each irrep')
        
    myMatrixCorrelatorData.close()
    
    
elif myRuns.correlator=='mr':
    myNonInteractingLevels = ed.ensembles[myRuns.ensemble][myRuns.isospin]['nonInteractingLevels']
    myRatioCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_ratios_' + myRuns.rs_type + reBin +'_v%s.h5'%myVersion,'r')
    myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/Plots/%s/Matrices_Ratios/'%myRuns.ensemble +  '%s/'%myResamplingScheme)

    if myRuns.correlator: 
        pcorr.PlotMultiHadronCorrelators(myRatioCorrelatorData, myRuns.rs_type, myVersion, myRuns.fit_param.t0, myPlotLocation, reBin)

    if myRuns.effmass: 
        peff.PlotMultiHadronsEffectiveMasses(myRatioCorrelatorData, myResamplingScheme, myVersion, myRuns.fit_param.t0, myPlotLocation,reBin)
        
    if myRuns.fit:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_Ratios/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators_ratios_' + myRuns.rs_type + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        pfit.PlotMultiHadronsFits(myFitCorrelator, myRuns.fit_param.type_correlation, myRuns.fit_param.type_fit,  multiTMinsFitPlots, myRuns.fit_param.t0, myVersion, myPlotLocation, reBin)
        
        myFitCorrelator.close()

    if myRuns.join:
        print('Now all the plots are in one file')
        
    myRatioCorrelatorData.close()

### --------- PRINTS WHERE IT IS SAVED ---------
        
print('-'*(len(myPlotLocation)+1))
print('Correlator Plots saved : \n' + myPlotLocation)
print('_'*(len(myPlotLocation)+1))

