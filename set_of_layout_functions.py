import os
from pathlib import Path
import set_of_plot_functions as vfp

## Comments:
# This function checks if a directory exists, if not then it creats it.
def DIRECTORY_EXISTS(a_dir):
    if not os.path.isdir(a_dir):
        Path(a_dir).mkdir(parents=True, exist_ok=True)
        new_dir = a_dir
    else:
        new_dir=a_dir
    return new_dir
        
        
### Comments:
# This function prints the information of the ensemble and the type of correlators one have.
def INFO_PRINTING(the_corr_type, the_ensemble):
    if the_corr_type=='s':
        print('.............................................................................')
        print('                         SINGLE HADRON CORRELATORS')
        print('                               ENSEMBLE '+ the_ensemble)
        print('.............................................................................')
    elif the_corr_type=='m':
        print('.............................................................................')
        print('                         MULTIHADRON CORRELATORS')
        print('                               ENSEMBLE '+ the_ensemble)
        print('.............................................................................')
    elif the_corr_type=='mr':
        print('.............................................................................')
        print('                       MULTIHADRON RATIO CORRELATORS')
        print('                               ENSEMBLE '+ the_ensemble)
        print('.............................................................................')
        
    

### Comments:
# This function is intended to provide a meny with the choices of hadrons to do the dispersion relation plots
def HADRON_MENU(the_irreps):
    the_rest_irreps = [the_irrep for the_irrep in the_irreps if the_irrep.startswith("PSQ0_")]
    the_hadrons = sorted({vfp.IrrepInfo(the_irrep).Hadron for the_irrep in the_rest_irreps})
    
    if len(the_hadrons)==1:
        print("Available hadrons: ", the_hadrons)
        return the_hadrons
    else:
        print("Available hadrons:")
        for ii, hh in enumerate(the_hadrons):
            print(f"[{ii+1}] {hh}")
        print(f"[{len(the_hadrons)+1}] All")
            
        the_choice = input("Select hadrons (separated by commas): ")
        the_indices = [int(the_x.strip()) for the_x in the_choice.split(",")]
        
        if (len(the_hadrons) + 1) in the_indices:
            return the_hadrons
        else:
            return [the_hadrons[ii-1] for ii in the_indices]


### Comments:
def SINGLE_INFO_PRINTING(the_irreps, bin_analysis, corr_analysis):
    info_list = []

    if bin_analysis:
        print('Choose one irrep for the bin-size analysis with PSQ = 0.')
    if corr_analysis:
        print('Choose the irreps you want to include in the correlator analysis.')
    else: print('Choose the irreps you want to include in the effective mass- or fit analysis.')
    for i in range(len(the_irreps)):
        info_list.append(the_irreps[i] + ' = ' + str(i))
    print(*info_list, sep='\n')


### Comments:
# prints the name of the irreps and the index of how they appear in the file as well as the operator combinations and the corresponding indices
def IRREP_OP_INFO_PRINTING(the_new_archivo):

    the_irreps = list(the_new_archivo.keys())
    info_list = []
    print('Irreps and the operator combinations you can choose:')
    for i in range(len(the_irreps)):
        irrep_info = []
        op_list = list(the_new_archivo[the_irreps[i]].keys())
        irrep_info.append(the_irreps[i] + ' = ' + str(i))
        for ii in range(len(op_list)):
            irrep_info.append('      ' + op_list[ii] + ' = ' + str(ii))
        info_list.append(irrep_info)
    columns_per_row = 4
    column_width = 60  # Adjust for spacing, here taken 60 because then long operator lists can be shown correctly
    for i in range(0, len(info_list), columns_per_row):
        chunk = info_list[i:i + columns_per_row]
        # Print headlines
        for sublist in chunk:
            print(f"{sublist[0]:<{column_width}}", end='')
        print()
        # Find the tallest sublist
        max_height = max(len(sublist) for sublist in chunk) - 1
        # Print remaining elements
        for row in range(max_height):
            for sublist in chunk:
                if row + 1 < len(sublist):
                    print(f"{sublist[row + 1]:<{column_width}}", end='')
                else:
                    print(" " * column_width, end='')
            print()
        # Print separation between blocks
        print()  # Blank line for spacing
        

### Comments:
# Prints the index and the name of the irrep in the file and reminds one to choose a rest frame irrep
def IRREP_BIN_SIZE_INFO_PRINTING(the_irreps):

    info_list = []

    print('Choose one irrep for the bin-size analysis with PSQ = 0:')
    for i in range(len(the_irreps)):
        info_list.append(the_irreps[i] + ' = ' + str(i))
    print(*info_list, sep='\n')

if __name__=="__main__":
   print('Nothing to run here, unless you want to change something.')
