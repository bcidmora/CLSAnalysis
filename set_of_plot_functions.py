import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.linalg import fractional_matrix_power
from iminuit import Minuit
import re
import math

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


### ------------------  PLOTTING STUFF --------------------------------------------------

    
def PLOT_SINGLE_HADRON_NAMES(the_hadron_name):
    the_plot_hadron_name = ''
    if the_hadron_name[0]=='P':
        the_plot_hadron_name=r'$\pi$'
    elif the_hadron_name[0]=='k':
        the_plot_hadron_name=r'$K$'
    elif the_hadron_name[0]=='N':
        the_plot_hadron_name=r'$N$'
    elif the_hadron_name[0]=='L':
        the_plot_hadron_name=r'$\Lambda$'
    elif the_hadron_name[0]=='S':
        the_plot_hadron_name=r'$\Sigma$'
    elif the_hadron_name[0]=='X':
        the_plot_hadron_name=r'$\Xi$'
    elif the_hadron_name[0]=='D':
        the_plot_hadron_name=r'$D$'
    elif the_hadron_name[0]=='O':
        the_plot_hadron_name=r'$\Omega$'
    elif 'dm' in the_hadron_name:
        the_plot_hadron_name=r'$\Omega_{\Delta m_{s}}$'
    elif 'qed' in the_hadron_name:
        the_plot_hadron_name=r'$\Omega_{e^{2}}$'
    return the_plot_hadron_name


def SQUARED_MOM(the_mom_str):
    the_mod_str = list(the_mom_str.split(','))
    the_sqrd_mom = (int(the_mod_str[0][the_mod_str[0].index('(')+1:])**2) + (int(the_mod_str[1])**2) + (int(the_mod_str[2][:the_mod_str[2].index(')')])**2)
    return the_sqrd_mom


### Comments: This function searches the min value to plot the y-axis and not be too shifted.
def CHOOSING_YMIN_PLOT(the_mean_efm):
    the_ymin=0. 
    if not np.isinf(min(the_mean_efm)) and not np.isnan(min(the_mean_efm)):
        if min(the_mean_efm)<0.:
            the_ymin=the_mean_efm[0]/3
        elif the_mean_efm[int(2* (len(the_mean_efm)/3))]>(the_mean_efm[0]*1.5) or the_mean_efm[int(2* (len(the_mean_efm)/3))]<(the_mean_efm[0]*.5):
            the_ymin=(the_mean_efm[0]/2)*.65 
    else: 
        # the_ymin=0.
        the_ymin= the_mean_efm[int(len(the_mean_efm)/2)]*.68
    return the_ymin



### Comments: This function searches the max value to plot the y-axis and not be too shifted.
def CHOOSING_YMAX_PLOT(the_mean_corr):
    da_t =0
    while np.isinf(the_mean_corr[da_t]) or np.isnan(the_mean_corr[da_t]): 
        da_t+=1
    return the_mean_corr[da_t]*1.07


### Comments:
# This function receives an operator name and returns a string in a nice way to put in the plots. 
def OPERATORS_SH(operator_name):
    new_op = list(operator_name.split(' '))
    OperatorPlot = ''
    if 'GI{' not in operator_name:
        if new_op[0].lower()=='pion':
            OperatorPlot = 'P[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='kaon':
            OperatorPlot = 'k[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='nucleon':
            OperatorPlot = 'N[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='lambda':
            OperatorPlot = 'L[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='sigma':
            OperatorPlot = 'S[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='xi':
            OperatorPlot = 'X[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='dmeson':
            OperatorPlot = 'D[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='omega':
            OperatorPlot = 'O[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='mass_shift':
            OperatorPlot = 'dm[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='qed':
            OperatorPlot = 'qed[%s]'%new_op[-1].replace('_','')
    elif 'GI{' in operator_name:
        OperatorPlot = str(new_op[-2])
    return str(OperatorPlot)


### Comments:
# This function receives an operator name and returns a string in a nice way to put in the plots. 
def OPERATORS_SH_ISOSPIN(the_hadron):
    the_hadron_isolabel = ''
    if the_hadron.upper()=='P':
        the_hadron_isolabel = '1'
    elif the_hadron.upper()=='K':
        the_hadron_isolabel = r'$\frac{1}{2}$'
    elif the_hadron.upper()=='N':
        the_hadron_isolabel = r'$\frac{1}{2}$'
    elif the_hadron.upper()=='L':
        the_hadron_isolabel = '0'
    elif the_hadron.upper()=='S':
        the_hadron_isolabel = '1'
    elif the_hadron.upper()=='X':
        the_hadron_isolabel = r'$\frac{1}{2}$'
    return the_hadron_isolabel

def PLOT_HADRON_LABELINGS(the_irrep_name):
    the_irrep_name = the_irrep_name.replace(" ","")
    the_irrep_name_plot = ''
    if the_irrep_name=="G1g": 
        the_irrep_name_plot = r'$G_{1g}$'
    elif the_irrep_name=="G1u": 
        the_irrep_name_plot = r'$G_{1u}$'
    elif the_irrep_name=="G2g": 
        the_irrep_name_plot = r'$G_{2g}$'
    elif the_irrep_name=="G2u": 
        the_irrep_name_plot = r'$G_{2u}$'
    elif the_irrep_name=="Hg": 
        the_irrep_name_plot = r'$H_{g}$'
    elif the_irrep_name=="Hu": 
        the_irrep_name_plot = r'$H_{u}$'
    elif the_irrep_name=="H": 
        the_irrep_name_plot = r'$H$'
    elif the_irrep_name=="G1": 
        the_irrep_name_plot = r'$G_{1}$'
    elif the_irrep_name=="G2": 
        the_irrep_name_plot = r'$G_{2}$'
    elif the_irrep_name=="G": 
        the_irrep_name_plot = r'$G$'
    elif the_irrep_name=="F1": 
        the_irrep_name_plot = r'$F_{1}$'
    elif the_irrep_name=="F2": 
        the_irrep_name_plot = r'$F_{2}$'
    elif the_irrep_name=="A1um": 
        the_irrep_name_plot = r'$A_{1u}^{-}$'
    elif the_irrep_name=="A1u": 
        the_irrep_name_plot = r'$A_{1u}$'
    elif the_irrep_name=="A1g": 
        the_irrep_name_plot = r'$A_{1g}$'
    elif the_irrep_name=="A2m": 
        the_irrep_name_plot = r'$A_{2}^{-}$'
    elif the_irrep_name=="A2g": 
        the_irrep_name_plot = r'$A_{2g}$'
    elif the_irrep_name=="A1u": 
        the_irrep_name_plot = r'$A_{2u}$'
    elif the_irrep_name=="Eg": 
        the_irrep_name_plot = r'$E_{g}$'
    elif the_irrep_name=="Eu": 
        the_irrep_name_plot = r'$E_{u}$'
    elif the_irrep_name=="T1g": 
        the_irrep_name_plot = r'$T_{1g}$'
    elif the_irrep_name=="T1u": 
        the_irrep_name_plot = r'$T_{1u}$'
    elif the_irrep_name=="T2g": 
        the_irrep_name_plot = r'$T_{2g}$'
    elif the_irrep_name=="T2u": 
        the_irrep_name_plot = r'$T_{2u}$'
    elif the_irrep_name=="A1": 
        the_irrep_name_plot = r'$A_{1}$'
    elif the_irrep_name=="A2": 
        the_irrep_name_plot = r'$A_{2}$'
    elif the_irrep_name=="B1": 
        the_irrep_name_plot = r'$B_{1}$'
    elif the_irrep_name=="E": 
        the_irrep_name_plot = r'$E$'
    elif the_irrep_name=="B2": 
        the_irrep_name_plot = r'$B_{2}$'
    return the_irrep_name_plot





### Comments:
# This function receives an operator name and returns a string in a nice way to put in the plots. 
# def OPERATORS_MH(the_operator_name):
#     new_op = list(the_operator_name.split(' '))
#     if "CG" in the_operator_name: the_shift = 1
#     else: the_shift=0
#     OperatorPlot = ''
#     if 'GI{' not in the_operator_name:
#         if '_' in new_op[0]:
#             the_hads = list(new_op[0].split('_'))
#             for ii in range(1,len(the_hads)):
#                 the_mom = str(SQUARED_MOM(new_op[2+the_shift+((ii-1)*3)]))
#                 the_irrep = str(new_op[3+the_shift+((ii-1)*3)])
#                 the_site = str(new_op[4+the_shift+((ii-1)*3)][:-1]).replace('_','')
#                 if the_hads[ii].lower()=='pion':
#                     OperatorPlot += 'P[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#                 elif the_hads[ii].lower()=='kaon' or the_hads[ii].lower()=='kbar':
#                     OperatorPlot += 'k['+ the_mom + '_' + the_irrep + '_' + the_site + ']'
#                 elif the_hads[ii].lower()=='nucleon':
#                     OperatorPlot += 'N[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#                 elif the_hads[ii].lower()=='lambda':
#                     OperatorPlot += 'L[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#                 elif the_hads[ii].lower()=='sigma':
#                     OperatorPlot += 'S[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#                 elif the_hads[ii].lower()=='xi':
#                     OperatorPlot += 'X[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#                 elif the_hads[ii].lower()=='dmeson':
#                     OperatorPlot += 'D[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#             
#         elif '_' not in new_op[0]:
#             the_mom = str(SQUARED_MOM(new_op[1]))
#             if '_' not in new_op[2]:           
#                     the_irrep = str(new_op[2]) 
#             else:
#                 the_irrep = str(new_op[2][:new_op[2].index('_')]) 
#             the_site = str(new_op[-1]).replace('_','')
#             if new_op[0].lower()=='pion':
#                 OperatorPlot = 'P[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#             elif new_op[0].lower()=='kaon' or new_op[0].lower()=='kbar':
#                 OperatorPlot = 'k['+ the_mom + '_' + the_irrep + '_' + the_site + ']'
#             elif new_op[0].lower()=='nucleon':
#                 OperatorPlot = 'N[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#             elif new_op[0].lower()=='lambda':
#                 OperatorPlot = 'L[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#             elif new_op[0].lower()=='sigma':
#                 OperatorPlot = 'S[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#             elif new_op[0].lower()=='xi':
#                 OperatorPlot = 'X[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#             elif new_op[0].lower()=='dmeson':
#                 OperatorPlot = 'D[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
#     elif 'GI{' in the_operator_name:
#         if '_' in new_op[4]:
#             OperatorPlot = new_op[4]
#         else:
#             if 'P=(' in new_op[2]:
#                 the_mom = str(SQUARED_MOM(new_op[2]))
#             else: 
#                 the_mom = str(new_op[2][new_op[2].index('='):])
#             the_irrep = new_op[3]
#             if '}' in new_op[4]: new_op[4] = new_op[4][:-1]
#             OperatorPlot = new_op[4][:new_op[4].index('[')+1] + the_mom + '_' + the_irrep + '_' + new_op[4][new_op[4].index('[')+1:-1] + ']'
#     return str(OperatorPlot)




def OPERATORS_MH(the_operator_name):
    new_op = list(the_operator_name.split(' '))
    if "CG" in the_operator_name: 
        the_shift = 1
    else: 
        the_shift = 0
    OperatorPlot = ''
    
    if 'GI{' not in the_operator_name:
        if '_' in new_op[0]:
            the_hads = list(new_op[0].split('_'))
            for ii in range(1,len(the_hads)):
                the_mom = str(SQUARED_MOM(new_op[2+the_shift+((ii-1)*3)]))
                the_irrep = str(new_op[3+the_shift+((ii-1)*3)])
                
                # Only remove trailing ']' if present
                the_site = str(new_op[4+the_shift+((ii-1)*3)])
                if the_site.endswith(']'):
                    the_site = the_site[:-1]
                the_site = the_site.replace('_','')
                
                if the_hads[ii].lower()=='pion':
                    OperatorPlot += 'P[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='kaon' or the_hads[ii].lower()=='kbar':
                    OperatorPlot += 'k['+ the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='nucleon':
                    OperatorPlot += 'N[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='lambda':
                    OperatorPlot += 'L[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='sigma':
                    OperatorPlot += 'S[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='xi':
                    OperatorPlot += 'X[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='dmeson':
                    OperatorPlot += 'D[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
        
        elif '_' not in new_op[0]:
            the_mom = str(SQUARED_MOM(new_op[1]))
            if '_' not in new_op[2]:           
                the_irrep = str(new_op[2]) 
            else:
                the_irrep = str(new_op[2][:new_op[2].index('_')]) 
            
            # Only remove trailing ']' if present
            the_site = str(new_op[-1])
            if the_site.endswith(']'):
                the_site = the_site[:-1]
            the_site = the_site.replace('_','')
            
            if new_op[0].lower()=='pion':
                OperatorPlot = 'P[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='kaon' or new_op[0].lower()=='kbar':
                OperatorPlot = 'k['+ the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='nucleon':
                OperatorPlot = 'N[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='lambda':
                OperatorPlot = 'L[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='sigma':
                OperatorPlot = 'S[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='xi':
                OperatorPlot = 'X[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='dmeson':
                OperatorPlot = 'D[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
    
    elif 'GI{' in the_operator_name:
        if '_' in new_op[4]:
            OperatorPlot = new_op[4]
        else:
            if 'P=(' in new_op[2]:
                the_mom = str(SQUARED_MOM(new_op[2]))
            else: 
                the_mom = str(new_op[2][new_op[2].index('='):])
            the_irrep = new_op[3]
            
            # Only remove '}' or ']' if present
            if new_op[4].endswith('}'):
                new_op[4] = new_op[4][:-1]
            if new_op[4].endswith(']'):
                new_op[4] = new_op[4][:-1]
            
            OperatorPlot = new_op[4][:new_op[4].index('[')+1] + the_mom + '_' + the_irrep + '_' + new_op[4][new_op[4].index('[')+1:] + ']'
    
    return str(OperatorPlot)



def SH_OPERATORS_RELABEL(a_string_ops, an_irrep,the_mom):
    # print(a_string_ops)
    hadron_raw = a_string_ops.split('[')[0]
    ss_tag = a_string_ops.split('[')[1][:-1]     
    ss_number = ss_tag.replace("SS", "")  
    hadron_tex = PLOT_SINGLE_HADRON_NAMES(hadron_raw)
    irrep_clean = an_irrep.strip('$')
    hadron_clean = hadron_tex.strip('$')
    if r'\Omega' not in hadron_clean and r'\Delta m_{s}' not in hadron_clean and r'e^{2}' not in hadron_clean:
        return rf'$ {hadron_clean}[{irrep_clean}({the_mom})]_{{\mathrm{{{ss_number}}}}} $'
    else:
        return rf'$ {hadron_clean}[{irrep_clean}({the_mom})] $'
    
    
    

# def MH_OPERATORS_RELABEL(a_string_ops):
#     def repl(m):
#         the_hadron, the_inside = m.group(1), m.group(2)
#         the_parts = the_inside.split('_')
#         if len(the_parts) < 3:
#             return m.group(0)
#         a_num, a_irrep, a_suffix = the_parts[0], the_parts[1], '_'.join(the_parts[2:])
#         a_suffix = a_suffix.replace("SS", "")
#         the_hadron = PLOT_SINGLE_HADRON_NAMES(the_hadron).strip('$')
#         a_irrep = PLOT_HADRON_LABELINGS(a_irrep).strip('$')
#         return fr"${the_hadron}[{a_irrep}({a_num})]_{{\mathrm{{{a_suffix}}}}}$"
#     return re.sub(r'([A-Za-z])\[(.*?)\]', repl, a_string_ops)



def MH_OPERATORS_RELABEL(a_string_ops):
    def try_fix_inside(inside):
        if inside.count('_') >= 2:
            return inside
        m = re.search(r'(SS\d+)$', inside)
        if not m:
            return inside
        return inside[:m.start(1)] + "_" + m.group(1)

    def repl(m):
        the_hadron, the_inside = m.group(1), m.group(2)
        the_inside = the_inside.lstrip('=') 
        the_inside = try_fix_inside(the_inside)
        the_parts = the_inside.split('_')
        if len(the_parts) < 3:
            return m.group(0)
        a_num, a_irrep, a_suffix = the_parts[0], the_parts[1], '_'.join(the_parts[2:])
        a_suffix = a_suffix.replace("SS", "")
        a_suffix = a_suffix.replace("SD", "")
        the_hadron = PLOT_SINGLE_HADRON_NAMES(the_hadron).strip('$')
        a_irrep = PLOT_HADRON_LABELINGS(a_irrep).strip('$')
        return fr"${the_hadron}[{a_irrep}({a_num})]_{{\mathrm{{{a_suffix}}}}}$" ## FOR THE NOISES
        # return fr"${the_hadron}({a_num})_{{\mathrm{{{a_suffix}}}}}$"
        # return fr"${the_hadron}({a_num})$"

    return re.sub(r'([A-Za-z])\[(.*?)\]', repl, a_string_ops)


### Comments:
# This class rewrites the names of the irreps in TeX type of text such that they can be put in the plots in a nice way.
class IrrepInfo:
     def __init__(self,nombre):
        self.name = nombre.split('_')            
        self.Name = self.name[1]
        self.NamePlot = PLOT_HADRON_LABELINGS(self.Name)
        self.Momentum = self.name[0][-1]
        self.Hadron = self.name[2]
        self.HadronIsospin = OPERATORS_SH_ISOSPIN(self.Hadron)
        # self.TotalMomPlot = self.name[0][0]+ r'$^{2}=%s$'%self.name[0][-1]
        # self.TotalMomPlot = r'$\vec{\mathbf{d}}^{2}=%s$'%self.name[0]
        self.TotalMomPlot = r'$\vec{\mathbf{d}}^{2}=%s$'%self.Momentum
        
        
        
### Comments:
# This class gets the info of the non-interacting levels. It only works when there is a threshold of 2 states nearby. It needs modification for the 3particle threshold.
class NonInteractingLevels:
    def __init__(self,nombre,all_hads):
        self.FirstState = str(nombre[:4])
        self.SecondState =  str(nombre[4:])
        FirstRaw = 'PSQ'+ str(self.FirstState[2])+'__'+ str(self.FirstState[0])
        SecondRaw = 'PSQ'+ str(self.SecondState[2])+'__'+ str(self.SecondState[0])
        for item in all_hads:
            if FirstRaw[:4] in item and FirstRaw[-1]==item[-1]: first_name = item;break
            else: continue
        for item in all_hads:
            if SecondRaw[:4] in item and SecondRaw[-1]==item[-1]: second_name = item;break
            else: continue
        self.First = first_name
        self.Second = second_name


### Comments:
# This function writes the errors in the plot in the way of parenthesis according to the precision given.
def WRITTING_ERRORS_PLOTS(an_error, the_precision):
    out_precision = True
    the_new_precision = the_precision
    the_error_string = "("
    if an_error>=1:
        print("The error is of order 1")
        out_precision=False
        an_error = str(f'{np.round(an_error, 1):.{1}f}')
        the_error_string+=an_error
        the_new_precision = 1
    else:
        an_error = str(f'{np.round(an_error, the_precision):.{the_precision}f}')
        if len(an_error)>(the_precision+2): 
            an_error=an_error[:the_precision+2]
        for ii in range(len(an_error)):
            if len(the_error_string)>=2:
                the_error_string+=an_error[ii]
            else:
                if an_error[ii]!="0" and an_error[ii]!=".": 
                    the_error_string+=an_error[ii]
                else: continue
        if len(the_error_string)>3:
            while len(the_error_string)>3:
                the_error_string= 10 * round(int(the_error_string[1:])/10)
                the_error_string='('+str(the_error_string)[:-1]
                out_precision=False
                the_new_precision=the_new_precision-1
    return [the_error_string+")", out_precision, the_new_precision]


def PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, the_x_axis_label, the_y_axis_label, the_marker, the_title_info, **kwargs):
    the_min_position = np.where(the_mean_corr == min(the_mean_corr[:-3]))
    the_max_position = np.where(the_mean_corr == max(the_mean_corr[:-3]))
    
    the_min_y = (the_mean_corr[the_min_position]-the_sigmas_corr[the_min_position])*.95
    the_max_y= (the_mean_corr[the_max_position]+the_sigmas_corr[the_max_position])*1.05 
    
    plt.errorbar(the_nt, the_mean_corr, yerr = the_sigmas_corr, marker=the_marker, ls='None', ms=4, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=2.85, label = the_rs_scheme, color='#5d83d5')
    plt.xlabel(the_x_axis_label,fontsize=28)
    plt.ylabel(the_y_axis_label,fontsize=28)
    plt.title(the_title_info,fontsize=20)
    plt.xticks(the_nt_ticks,fontsize=18)
    plt.yticks(fontsize=18)
    if kwargs.get('yscale')!=None: plt.yscale(str(kwargs.get('yscale')))
    else:
        plt.ylim([the_min_y, the_max_y])
        # if kwargs.get('ymin')!=None:
            # plt.ylim(ymin=kwargs.get('ymin'), ymax=CHOOSING_YMAX_PLOT(the_mean_corr)*1.05)
    plt.legend(fontsize=16, handletextpad=0.01)
    plt.tight_layout()
    # plt.show()
    
    
def PLOT_HISTOGRAMS(the_rs, the_label , the_mean_rs, the_label_mean_rs, the_nt_mean, the_label_mean_nt, the_title_info, the_bins,  the_x_axis_label):
    counts, bins, patches = plt.hist(the_rs, bins=the_bins, label =  the_label, color='#5d83d5')
    padding = counts.max() * 0.1  # 10% padding on top
    plt.vlines(the_mean_rs, 0, 200, colors= '#b90f22', label = the_label_mean_rs)
    plt.vlines(the_nt_mean, 0, 200, colors='black', label = the_label_mean_nt)
    plt.title( the_title_info,fontsize=20)
    plt.ylabel('Frequency',fontsize=24)
    plt.xlabel(the_x_axis_label, fontsize=24)
    plt.legend(fontsize=14, handletextpad=0.01)
    plt.tight_layout()
    plt.ylim(0, counts.max() + padding)
    # plt.show()



def PLOT_FITS(the_nt, the_plot_data, the_sigmas_data, the_chosen_tmin, the_label, the_xlabel, the_ylabel, the_title, the_nt_ticks, **kwargs):
    if kwargs.get('zoom'):
        the_ll = int(kwargs.get('the_ll'))
        the_ul = int(kwargs.get('the_ul'))
        plt.errorbar(the_nt[the_chosen_tmin-the_ll:the_chosen_tmin+the_ul], the_plot_data[the_chosen_tmin-the_ll:the_chosen_tmin+the_ul], yerr = the_sigmas_data[the_chosen_tmin-the_ll:the_chosen_tmin+the_ul], marker='o', ls='None', ms=6, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=6, color='#5d83d5')#, label='Bootstrap')
    else:
        plt.errorbar(the_nt, the_plot_data, yerr = the_sigmas_data, marker='o', ls='None', ms=6, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=6,color='#5d83d5')#, label='Bootstrap')
        plt.xticks(the_nt_ticks,fontsize=18)    
    plt.errorbar([the_nt[the_chosen_tmin]], [the_plot_data[the_chosen_tmin]], yerr = [the_sigmas_data[the_chosen_tmin]], marker='o', ls='None', ms=7, markeredgewidth=1.75, lw=1.75, color = '#b90f22', elinewidth=1.75, zorder=3, markerfacecolor = 'white', capsize=6, label = the_label)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=16, ncols=2, columnspacing=0.1,handletextpad=0.01)
    plt.xlabel(the_xlabel,fontsize=36)
    plt.ylabel(the_ylabel,fontsize=36)
    plt.title(the_title,fontsize=26)
    # plt.xticks(the_nt_ticks)
    plt.tight_layout()




def PLOT_CHI_FITS(the_nt, the_plot_data, the_chosen_tmin, the_label, the_xlabel, the_ylabel, the_title, the_nt_ticks, **kwargs):
    plt.plot(the_nt, the_plot_data, marker='o', ls='None', ms=4, markeredgewidth=1.75, lw=1.75, zorder=3, color='#5d83d5')
    plt.plot([the_nt[the_chosen_tmin]], [the_plot_data[the_chosen_tmin]], marker='o', ls='None', ms=4, markeredgewidth=1.75, lw=1.75, color = '#b90f22', zorder=3, markerfacecolor = 'white', label = the_label)
    plt.xticks(the_nt_ticks)
    plt.legend(fontsize=14, handletextpad=0.01)
    plt.xlabel(the_xlabel,fontsize=16)
    plt.ylabel(the_ylabel,fontsize=16)
    plt.title(the_title,fontsize=16)
    # plt.xticks(the_nt_ticks)
    plt.tight_layout()
    
    
import matplotlib.ticker as mtick


def PLOT_FITTED_EFF_MASSES(the_nt, the_mean_corr, the_sigmas_corr, the_fit_data, the_fit_sigmas, the_chosen_tmin,
                           the_rs_scheme, the_label, the_title, the_nt_ticks, the_color_eff_mass, the_color_fit):
    the_min_position = np.where(the_mean_corr == min(the_mean_corr[:-3]))
    the_max_position = np.where(the_mean_corr == max(the_mean_corr[:-3]))

    the_min_y = (the_mean_corr[the_min_position] - the_sigmas_corr[the_min_position]) * .95
    the_max_y = (the_mean_corr[the_max_position] + the_sigmas_corr[the_max_position]) * 1.05

    plt.errorbar(the_nt, the_mean_corr, yerr=the_sigmas_corr, marker='o', ls='None', ms=8, markeredgewidth=1.75,
                 lw=1.75, elinewidth=1.75, capsize=5, label=the_rs_scheme, color=the_color_eff_mass)

    plt.axhline(y=the_fit_data[the_chosen_tmin], color=the_color_fit, ls='-', lw=1.75, label=the_label)

    plt.fill_between(x=[the_nt[0] - 1, the_nt[-1] + 3],
                     y1=the_fit_data[the_chosen_tmin] - the_fit_sigmas[the_chosen_tmin],
                     y2=the_fit_data[the_chosen_tmin] + the_fit_sigmas[the_chosen_tmin],
                     color=the_color_fit,
                     alpha=0.2)
    plt.xlabel(r'$t\,/\, a$', fontsize=36)
    plt.ylabel(r'$a_{t}\,m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=36)
    plt.title(the_title, fontsize=26)
    plt.xticks(the_nt_ticks, fontsize=18)
    plt.ylim([the_min_y, the_max_y])
    plt.yticks(fontsize=18)
    # plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    # # plt.xlim([the_nt[1]*.75,len(the_nt)])
    plt.xlim([the_nt[0] - 1, the_nt_ticks[-1] + 1])
    plt.legend(fontsize=18, handletextpad=0.3)
    plt.tight_layout()
    # plt.show()
    
    
    
    
## ----------------------------------------------- PLOT DESIGN STUFF FOR THE BIN SIZE AND THE DISPERSION RELATION --------------------------------

### All the colors one wants to use
pink = (238/255, 18/255, 137/255)
violet = (224 / 255, 102 / 255, 255 / 255)
dark_violet = (157 / 255, 50 / 255, 168 / 255)
green = (50 / 255, 205 / 255, 50 / 255)
orange_light = (255 / 255, 140 / 255, 0)
blue = (100 / 255, 149 / 255, 237 / 255)
dark_blue = (0, 0, 205 / 255)
brown = (139 / 255, 69 / 255, 19 / 255)
orange_dark = (255 / 255, 69 / 255, 0)
green_dark = (0, 100 / 255, 0)
türkis = (0, 250 / 255, 154 / 255)
rosa = (238 / 255, 180 / 255, 180 / 255)
dark_grey = (77 / 255, 77 / 255, 77 / 255)
mid_grey = (130 / 255, 130 / 255, 130 / 255)
light_grey = (191 / 255, 191 / 255, 191 / 255)

### Colors and markers need to appear in those lists
# colors = [pink, violet, dark_violet, green, dark_blue, orange_light, blue]
colors  = [ "#5d83d5", "#b90f22", "#ffa500", "#008000", "#c44601", "#f57600", "#5ba300","#e6308a", "#8a2be2", "#00ced1", "#ffd700", "#ff69b4", "#7cfc00", "#dc143c", "#4682b4", "#ff8c00", "#00fa9a", "#9370db", "#1e90ff", "#ff1493", "#9acd32"]
all_colors = [pink, violet, dark_violet, green, dark_blue, orange_light, blue, brown, orange_dark, green_dark, türkis, rosa, mid_grey, dark_grey]
special_colors =[dark_violet, orange_light]
greys = [dark_grey, mid_grey, light_grey]
markers = ['o', 's', '^', 'v', '_', 'd', '>', '<', 'p', 'h', 'x', '+', '*', 'P', 'X']
line_style = ['-', '--', '-.', ':']
special_markers=['^', 'd']

# sors the latex labels for the particles
particel_dictionary = {'K': '$K$',
                       'P': r'$\pi$'}

### Comments:
# Functions to round the numbers appearing in plots to the given digit
def ROUND_UP_3(x):
    factor = 10**3
    return math.ceil(x * factor) / factor

def ROUND_UP_2(x):
    factor = 10**2
    return math.ceil(x * factor) / factor

def ROUND_UP_4(x):
    factor = 10**4
    return math.ceil(x * factor) / factor

# Functions unseful for Plot titles
def GET_IRREP_LOGO(text):
    momentum_j = text[3]
    momentum_j_label = r'$P^2$ = ' + str(momentum_j)
    return momentum_j_label

### Comments:
def GET_RESAMPLING_BINNING(file):
    rs = file.filename.split('/')[-1].split('.')[0].split('_')[2]
    bin_raw = file.filename.split('/')[-1].split('.')[0].split('_')[3][3:]
    bin = str()
    if rs =='jk': rs='Jackknife'
    elif rs =='bt': rs='Bootstrap'
    if bin_raw == '1' : bin='no rebinning'
    elif bin_raw == 'E250': bin= r'$N_{\mathrm{bin}}$ = 1.5'
    else: bin= r'$N_{\mathrm{bin}}$ = ' + str(bin_raw)
    return rs, bin, bin_raw

### Comments:
# Converts A1 to A_1, T1u to T_{1u} etc.
def IRREP_TO_INDEX(name: str) -> str:

    if len(name) == 1:
        return name

    base = name[0]
    index = name[1:]

    # Use braces when index has more than one character
    if len(index) > 1:
        return f"${base}_{{{index}}}$"
    else:
        return f"${base}_{index}$"


def CONVERT_LABEL_TO_LATEX(level_str, label_dict):
    parts = re.findall(r'([A-Za-z]+)\((\d+)\)', level_str)
    latex_parts = [f"{label_dict.get(p, p)}({n})" for p, n in parts]
    return ''.join(latex_parts)


if __name__=="__main__":
   print('Nothing to run here, unless you want to change something.')

