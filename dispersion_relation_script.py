import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import h5py
import statistics
import numpy as np
import math
import os
from datetime import datetime
import glob

# from pandas.core.interchange.from_dataframe import primitive_column_to_ndarray

import set_of_plot_functions as vfp
import matplotlib.gridspec as gridspec
import os
import matplotlib.ticker as ticker
import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp

plt.rcParams["font.family"] = "sans"
plt.rcParams["mathtext.fontset"] = "dejavusans"

def DISPERSION_RELATION(fit_file, the_rs, the_tmin_choice, the_lat_size, the_disp_location, the_plot_type, the_iso_label, the_ensemble):

    print("\n                     DISPERSION RELATION \n")

    now = datetime.now()
    version = now.strftime("%d_%m_%H%M")

    irreps = list(fit_file.keys())
    norm = (2 * np.pi / the_lat_size) ** 2

    energies, energies_sampled, momenta = [], [], []

    momenta = sorted([int(irrep[3]) for irrep in irreps])
    mom_expected = set(range(min(momenta), max(momenta)+1))
    mom_expected.discard(7)
    mom_obtained = set(momenta)
    missing = mom_expected - mom_obtained
    duplicates = [m for m in set(momenta) if momenta.count(m) > 1]

    assert not missing, f'The momentum p^2 = {missing} irrep is missing in your fit file.'
    assert len(momenta) == len(the_tmin_choice), f'Length of chosen TMin list (len = {len(the_tmin_choice)}) does not match the number of different momenta (len = {len(momenta)}).'
    assert not duplicates, f'There are too many irreps with the same momentum squared p^2 = {duplicates}.'

    tmin_dict = dict(zip(momenta, the_tmin_choice))
    x_ticks_dict = {i: m for i, m in enumerate(momenta)}
    x_pos = list(x_ticks_dict.keys())


    for m in momenta:

        try:
            match = next(irrep for irrep in irreps if f'PSQ{m}' in irrep)
        except StopIteration:
            raise ValueError(f'Digit {m} occurs in more than one irrep. Routine only works for momenta less equal 9.')

        t_min_ll = int(list(fit_file[match].get('1exp/Tmin/Correlated/Mean'))[0][0])

        energy_i = list(fit_file[match].get('1exp/Tmin/Correlated/Mean'))[2][tmin_dict[m] - t_min_ll]
        energy_i_sample = list(fit_file[match].get('1exp/Tmin/Correlated/Resampled'))[tmin_dict[m] - t_min_ll]

        energies.append(energy_i)
        energies_sampled.append(energy_i_sample)

    energies = np.array(energies)
    energies_sampled = np.array(energies_sampled)
    momenta = np.array(momenta)

    if the_plot_type == 'abs' or the_plot_type == 'both':

        energy_dispersion = []
        error_dispersion = []
        energy_calculated = []
        energies_sampled_normalized = []

        for k in range(momenta.size):
            energy_dispersion.append(vfa.CONTINUUM_DISP_REL(momenta[k], energies[0], norm))

            energy_calculated.append(energies[k] / energies[0])

            useful = []
            useful_2 = []
            for kk in range(energies_sampled[k].size):
                useful.append(energies_sampled[k][kk] / energies_sampled[0][kk])
                useful_2.append(vfa.CONTINUUM_DISP_REL(momenta[k], energies_sampled[0][kk], norm))
            energies_sampled_normalized.append(useful)
            error_dispersion.append(useful_2)

        sigmas_dispersion = []
        sigmas_calculated = []
        for k in range(momenta.size):
            sigmas_dispersion.append(vfa.STD_DEV(error_dispersion[k], energy_dispersion[k], the_rs))
            sigmas_calculated.append(vfa.STD_DEV(energies_sampled_normalized[k], energy_calculated[k], the_rs))

        data_min = []
        data_max = []
        range_size = []
        center = []

        for i in range(momenta.size):
            mini = min(energy_calculated[i] + sigmas_calculated[i], energy_calculated[i] - sigmas_calculated[i],
                       energy_dispersion[i] + sigmas_dispersion[i], energy_dispersion[i] - sigmas_dispersion[i])
            maxi = max(energy_calculated[i] + sigmas_calculated[i], energy_calculated[i] - sigmas_calculated[i],
                       energy_dispersion[i] + sigmas_dispersion[i], energy_dispersion[i] - sigmas_dispersion[i])
            center_i = np.mean([mini, maxi])
            center.append(center_i)
            data_min.append(mini)
            data_max.append(maxi)
        data_min = np.array(data_min)
        data_max = np.array(data_max)
        range_size = max(list(data_max - data_min))

        if len(momenta) < 6:
            fig = plt.figure(figsize=(10, 7))
        else: fig = plt.figure(figsize=(21, 7))
        cols = len(momenta)
        gs = gridspec.GridSpec(2, cols, height_ratios=[1, 1])  # 2 Reihen, 4 Spalten
        ax_top = fig.add_subplot(gs[0, :])
        axes_bottom = [fig.add_subplot(gs[1, i]) for i in range(cols)]  # Untere Reihe: 4 einzelne Plots
        for i, ax in enumerate(axes_bottom):
            ax.errorbar(x_pos[i] , energy_dispersion[i], yerr=sigmas_dispersion[i], fmt='^', markersize=8,
                        capsize=4,
                        color=vfp.colors[4], ecolor=vfp.greys[1])
            ax.errorbar(x_pos[i], energy_calculated[i], yerr=sigmas_calculated[i], fmt='o', markersize=8,
                        capsize=4, color=vfp.colors[0], ecolor=vfp.greys[0])

            ax.set_ylim(center[i] - range_size / 2 - 0.005, center[i] + range_size / 2 + 0.005)

            ax.tick_params(axis='both', which='major',
                           length=10, direction='inout', labelsize=16, top=True, right=True)
            ax.set_xticks([x_pos[i]])
            ax.set_xticklabels([momenta[i]])

        ax_top.errorbar(x_pos, energy_dispersion, yerr=sigmas_dispersion, fmt='^', markersize=8,
                        capsize=4,
                        color=vfp.colors[4], ecolor=vfp.greys[1], label='$D(am, p^2)$')
        ax_top.errorbar(x_pos, energy_calculated, yerr=sigmas_calculated, fmt='o', markersize=8,
                        capsize=4,
                        color=vfp.colors[0], ecolor=vfp.greys[0], label='Fit results')

        ax_top.set_xlabel(r'$p^2$', fontsize=20)

        ax_top.set_ylabel(r'$aE$ / $am$', fontsize=20)
        axes_bottom[0].set_ylabel(r'$aE$ / $am$', fontsize=20)

        ax_top.set_title('Dispersion Relation', fontsize=20)
        ax_top.legend(fontsize=14)
        ax_top.tick_params(axis='both', which='major',
                           length=10, direction='inout', labelsize=16, top=True, right=True)
        ax_top.set_xticks(x_pos, [x_ticks_dict[i] for i in x_pos])
        fig.suptitle(rf'$I$ = {the_iso_label}, ' + vfp.GET_RESAMPLING_BINNING(fit_file)[0] + ', ' +
                     vfp.GET_RESAMPLING_BINNING(fit_file)[1] + ', ' + the_ensemble, fontsize=20)

        fig.tight_layout()
        fig.savefig(the_disp_location + f'v_{version}_absolute_dispersion_relation.pdf', format='pdf')
        fig.show()


    if the_plot_type == 'rel' or the_plot_type == 'both':

        energies_relative = []
        energies_sampled_relative = []
        sigmas_relative = []

        for i in range(momenta.size):
            energies_relative.append(
                vfa.RELATIVE_DISP_REL(momenta[i], energies[0], energies[i], norm))
            useful = []
            for k in range(energies_sampled[i].size):
                useful.append(vfa.RELATIVE_DISP_REL(momenta[i], energies_sampled[0][k],
                                                           energies_sampled[i][k], norm))
            energies_sampled_relative.append(useful)
            sigmas_relative.append(
                vfa.STD_DEV(energies_sampled_relative[i], energies_relative[i], the_rs))


        fig, ax = plt.subplots(figsize=(8, 4))
        ax.errorbar(x_pos, energies_relative, yerr=sigmas_relative, fmt='o', markersize=6, capsize=3,
                    color=vfp.colors[0], ecolor=vfp.greys[0])
        ax.tick_params(axis='both', which='major',
                       length=10, direction='inout', labelsize=12, top=True, right=True)
        ax.set_xlabel(r'$p^2$', fontsize=16)
        ax.axhline(1, 0, 10, linestyle=':', color=vfp.greys[2])
        ax.set_ylabel(r'$aE$ / $D(am, p^2)$', fontsize=16)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.set_xticks(x_pos, [x_ticks_dict[i] for i in x_pos])
        ax.set_title('Dispersion Relation', fontsize=16)
        fig.suptitle(rf'$I$ = {the_iso_label}, ' + vfp.GET_RESAMPLING_BINNING(fit_file)[0] + ', ' +
                     vfp.GET_RESAMPLING_BINNING(fit_file)[1] + ', ' + the_ensemble, fontsize=16)
        fig.tight_layout()
        fig.savefig(the_disp_location + f'v_{version}_relative_dispersion_relation.pdf', format='pdf')
