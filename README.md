The set of scripts are separated by the main actions one can do: resampling, effective masses, GEVP, fits. Some notes before starting anything:

1. The ensembles.py file must be edited using the name of the correct files for the corresponding ensemble. For each ensemble, there are subdivisions that include different quantum numbers and separates the single hadron from the multi hadron correlators. This file contains the following:
  - Information about the ensemble (beta, extent, lattice spacing, nr. gauge configs, reweighting factors, etc.)
  - Information about the single hadrons (name of the file, TMins and TMaxs for the fits)
  - Information about the multi-hadrons (isospin, name of the file, TMins and TMaxs, non-interacting levels, operator choices, etc.)

2. The main_analysis_script.py is the file that contains the information for all the possible tasks to be done depending on the choice of correlators (single "s", multi "m" or ratio of correlators "mr"):
  - Resampling
  - Effective Masses
  - GEVP
  - Fits

3. The set_of_analysis_functions.py contains all the routines related to all the calculations.
4. The set_of_layout_functions.py includes all the routines that are related to showing information, looking for directories and related to the inputs to run the scripts. In here is defined what are the inputs it takes:
  - Ensemble (mandatory): --ensemble or -e 
  - Correlator (mandatory): --correlator or -c
  - Resampling scheme (mandatory): --rs-type or -rs
  - Isospin (only for multihadron choice): --isospin or -i
  - Binning: --rebin and -rb (both needed if binning is wanted)

  - Correlator analysis: --corrs
  - Effective Masses: --effmass (for "m" and "mr" it also needs --gepv-flag for the full set or --ops-flag for the reduced op set)
  - Fitting: --fits (it needs --fit-type, if not provided it does the 1-exp fit, --fit-correlation is needed, otherwise it performs the "Correlated" fir by default)
  - GEVP (full op set): --eigenvals (needs --t0min and --t0max)
  - GEVP (reduced op set): --ops (needs --t0min and --t0max)
  - Dispersion relation: --disp-rel
  - Correlated differences: --corr-dif (still under construction)
  - Binning analysis: --binning (Zoe's implementation)

  - Irreps: --start-irrep or -fi; --last-irrep or -li; --nr-irreps or -ir.

To run the script, the directories must be modified, after that one can do the following:
 
 python3 main_analysis_script.py -e x451 -c m -i s -rs bt -kbt 500 --rebin -rb 10 --corrs --effmass --eigenvals --t0min 3 --t0max 4 --fits --fit-type 1 -fi 1 -li 2

 Which means the following:
   * Ensemble: X451
   * Correlator: Multi-hadrons
   * Isospin: Isosinglet
   * Resampling scheme: Bootstrap with 500 samples
   * Rebin: Nbin size = 10
   * Does the Resampling (--corrs), GEVP of the full set (--eigenvals) with t0min and t0max, it computes the effective masses of the diagonal correlators and of the eigenvalues, does the fitting using 1-exp fit forms, and it starts on the first irrep 1 and it finishes on the second.
    


Here a description of the ensembles file:

- 'aLat' : lattice spacing
- 'betaLat' : coupling value
- 'LatSize' : lattice extent
- 'ncfgs' : Number of TOTAL configs of that ensemble
- 'allConfigs': If not all configs of ensemble, then False
- 'nfgsList': list of configs to process (notice they start form 0, so they are shifted by 1 unit)
- 'weight_raw' : reweighting factors file(s) as a list
- 'fs' : string with the location of the single hadrons hdf5 file
- 'singleTMaxFits' : Choices of tmax for the fits (list)
- 'singleTMinResults' :  Choices of Tmin for the fit-plots (list)
- 's' : { 'iso_name' : singlet, doublet, etc.
- 'iso_tag' : isosinglet, ...
- 'iso_label' : value of I
- 'fm' : string with the directory of the isosinglet hdf5 file
- 'multiTMaxFits' : [[]], List of Tmaxs for the fits for the full operator basis
- 'multiTMinResults' : [[]], List of Tmins for the fit-plots for the full operator basis
- 'nonInteractingLevels' : [[]], Non-interacting levels per irrep, look at the other ensembles with examples
- 'multiTMaxFitsRatios' : [[]], List of TMax for the fits to the ratio of correlators
- 'multiTMinResultsRatios' : [[]], List of tmins for the fit-plots of the ratio of correlators
- 'operatorsChoice': [[]], List of operators to keep for each irrep. If empty, then all are kept.
- 'multiTMaxFitsOpChoices' : [[]], List of Tmax for the fits of the selected operators (eigenvals)
- 'multiTMinResultsOpChoices' : [[]], List of tmins for the fit-plots of the eigenvalues derived form a selected fit }
