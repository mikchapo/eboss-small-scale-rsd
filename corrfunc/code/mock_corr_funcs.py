import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import sys
from seaborn import color_palette
from astropy.table import Table, vstack
from astropy.io import fits


def calc_step_size(binning):
    step_sizes = []
    for dim in binning:
        if dim[3] == 'log':
            step_size = (np.log10(dim[1]) - np.log10(dim[0])) / dim[2]
            
        else:
            step_size = (dim[1] - dim[0]) / dim[2]
        
        step_sizes.append(step_size)
        
    return step_sizes


def generate_empty_bins(binning):
    binning = np.array(binning, dtype=[('min', 'f4'), ('max', 'f4'), ('N_bin', 'i4'), ('scaling', 'U5')])
    # print(binning)
    column_bases = []
    for dim in binning:
        if dim['scaling'] == 'log':
            column_bases.append(np.logspace(np.log10(dim['min']), np.log10(dim['max']), dim['N_bin'], endpoint=False))
            
        elif dim['scaling'] == 'lin':
            column_bases.append(np.linspace(dim['min'], dim['max'], dim['N_bin'], endpoint=False))
                       
        else:
            print("Invalid scaling type")
            return False
    # print(column_bases)

    N_bins_tot = np.prod(binning['N_bin'])
    
    columns = []
    for i in range(binning.shape[0]):     
        column = np.tile(column_bases[i], (np.prod(binning['N_bin'][i+1:]),1))
        column = column.flatten(order='F')
        column = np.tile(column, np.prod(binning['N_bin'][:i]))
            
        columns.append(column)

    columns.append(np.zeros(N_bins_tot))   
    bin_array = np.column_stack(columns)
    # print(bin_array.shape)

    return bin_array


def rebin_pair_counts(original, compression_factor, binning):
    binning = np.array(binning, dtype=[('min', 'f4'), ('max', 'f4'), ('N_bin', 'i4'), ('scaling', 'U5')])
    new_binning = np.copy(binning)
    new_binning["N_bin"][0] = binning["N_bin"][0] / compression_factor
    rebinned = generate_empty_bins(new_binning)
    print("Shape of original:", original.shape)
    print("Shape of rebinned:", rebinned.shape)
    for i in range(binning["N_bin"][0]):
        i_new = int(i / compression_factor)
        for j in range(binning["N_bin"][1]):
            rebinned[i_new*new_binning["N_bin"][1]+j, 2] += (
                        original[i*binning["N_bin"][1]+j, 2])
    return rebinned


def landy_szalay(DD, DR, RR):
    corr_func = np.copy(DD)
    corr_func[:,2] = (np.divide(DD[:,2], RR[:,2], out=np.zeros_like(DD[:,2]), where=RR[:,2]!=0) - 2. *
                      np.divide(DR[:,2], RR[:,2], out=np.zeros_like(DR[:,2]), where=RR[:,2]!=0) + 1.)
    return corr_func


def L_0(x):
    return 1.


def L_2(x):
    return 0.5 * (3. * x * x - 1.)


def L_4(x):
    return (35.* x * x * x * x - 30. * x * x + 3.) / 8.


L_polynomials = [L_0, L_2, L_4]

def corr_func_multipole_calc(l, corr_func, mu_bins=100):
    L = L_polynomials[np.floor(l / 2).astype(int)]
    
    corr_func_multipole = np.zeros((np.unique(corr_func[:,0]).size, 2))
    corr_func_multipole[:,0] = np.unique(corr_func[:,0])
    for multi in corr_func_multipole:
        for corr in corr_func:
            if corr[0] == multi[0]:
    
                multi[1] += (1./ float(mu_bins)) * corr[2] * L(corr[1] + 0.5 / float(mu_bins))
    
    corr_func_multipole[:,1] = corr_func_multipole[:,1] * (2. * l + 1.)
    
    return corr_func_multipole


def corr_func_projected_calc(corr_func, dr_pi=1., r_pi_max=80.):
    corr_func_projected = np.zeros((np.unique(corr_func[:,0]).size, 2))
    corr_func_projected[:,0] = np.unique(corr_func[:,0])
#     corr_func_cutoff = np.zeros((np.unique(corr_func[:,0]).size, 2))
#     corr_func_cutoff[:,0] = np.unique(corr_func[:,0])
    
    for proj in corr_func_projected:
        for corr in corr_func:
            if corr[0] == proj[0] and corr[1] < r_pi_max:
                proj[1] += dr_pi * corr[2]
    
    corr_func_projected[:,1] = corr_func_projected[:,1] * 2.
    
    return corr_func_projected


base_output_path = "../output/mock_data_vectors/mock_"
base_data_path = "../data/"
mps_dir = "counts_mps_EZmocks_aemulus/"
wp_dir = "counts_wp_EZmocks_aemulus/"

original_mps_binning = [(0.1, 60., 180, 'log'),
                        (0., 1., 100, 'lin')]

original_wp_binning = [(0.1, 60., 180, 'log'),
                        (0., 120., 120, 'lin')]

mps_binning = [(0.1, 60., 9, 'log'),
                (0., 1., 100, 'lin')]

wp_binning = [(0.1, 60., 9, 'log'),
                (0., 120., 120, 'lin')]

mps_binning = np.array(mps_binning, dtype=[('min', 'f4'), ('max', 'f4'), ('N_bin', 'i4'), ('scaling', 'U5')])
wp_binning = np.array(wp_binning, dtype=[('min', 'f4'), ('max', 'f4'), ('N_bin', 'i4'), ('scaling', 'U5')])

wp_step_sizes = calc_step_size(wp_binning)
start_time = dt.datetime.now()


for i in range(1):
    print("Starting mock %i" % i)

    if i==0:
        print("Started loading pair counts", dt.datetime.now() - start_time)

    if i<9:
        with open(base_data_path + mps_dir + "DD_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-000%i_mike.dat" % (i+1), 'r') as f:
            DD_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DD_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-000%i_mike.dat" % (i+1), 'r') as f:
            DD_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DR_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-000%i_mike.dat" % (i+1), 'r') as f:
            DR_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DR_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-000%i_mike.dat" % (i+1), 'r') as f:
            DR_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "RR_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-000%i_mike.dat" % (i+1), 'r') as f:
            RR_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "RR_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-000%i_mike.dat" % (i+1), 'r') as f:
            RR_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')

        with open(base_data_path + wp_dir + "DD_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-000%i_aemulus.dat" % (i+1), 'r') as f:
            DD_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DD_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-000%i_aemulus.dat" % (i+1), 'r') as f:
            DD_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-000%i_aemulus.dat" % (i+1), 'r') as f:
            DR_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-000%i_aemulus.dat" % (i+1), 'r') as f:
            DR_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "RR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-000%i_aemulus.dat" % (i+1), 'r') as f:
            RR_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "RR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-000%i_aemulus.dat" % (i+1), 'r') as f:
            RR_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')

    elif i<99:
        with open(base_data_path + mps_dir + "DD_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-00%i_mike.dat" % (i+1), 'r') as f:
            DD_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DD_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-00%i_mike.dat" % (i+1), 'r') as f:
            DD_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DR_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-00%i_mike.dat" % (i+1), 'r') as f:
            DR_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DR_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-00%i_mike.dat" % (i+1), 'r') as f:
            DR_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "RR_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-00%i_mike.dat" % (i+1), 'r') as f:
            RR_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "RR_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-00%i_mike.dat" % (i+1), 'r') as f:
            RR_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')

        with open(base_data_path + wp_dir + "DD_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-00%i_aemulus.dat" % (i+1), 'r') as f:
            DD_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DD_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-00%i_aemulus.dat" % (i+1), 'r') as f:
            DD_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-00%i_aemulus.dat" % (i+1), 'r') as f:
            DR_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-00%i_aemulus.dat" % (i+1), 'r') as f:
            DR_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "RR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-00%i_aemulus.dat" % (i+1), 'r') as f:
            RR_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "RR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-00%i_aemulus.dat" % (i+1), 'r') as f:
            RR_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')

    elif i!=999:
        with open(base_data_path + mps_dir + "DD_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-0%i_mike.dat" % (i+1), 'r') as f:
            DD_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DD_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-0%i_mike.dat" % (i+1), 'r') as f:
            DD_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DR_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-0%i_mike.dat" % (i+1), 'r') as f:
            DR_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DR_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-0%i_mike.dat" % (i+1), 'r') as f:
            DR_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "RR_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-0%i_mike.dat" % (i+1), 'r') as f:
            RR_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "RR_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-0%i_mike.dat" % (i+1), 'r') as f:
            RR_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')

        with open(base_data_path + wp_dir + "DD_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-0%i_aemulus.dat" % (i+1), 'r') as f:
            DD_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DD_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-0%i_aemulus.dat" % (i+1), 'r') as f:
            DD_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-0%i_aemulus.dat" % (i+1), 'r') as f:
            DR_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-0%i_aemulus.dat" % (i+1), 'r') as f:
            DR_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "RR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-0%i_aemulus.dat" % (i+1), 'r') as f:
            RR_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "RR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-0%i_aemulus.dat" % (i+1), 'r') as f:
            RR_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
    else:
        with open(base_data_path + mps_dir + "DD_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-1000_mike.dat", 'r') as f:
            DD_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DD_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-1000_mike.dat", 'r') as f:
            DD_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DR_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-1000_mike.dat", 'r') as f:
            DR_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "DR_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-1000_mike.dat", 'r') as f:
            DR_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "RR_s-mu_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-1000_mike.dat", 'r') as f:
            RR_mps_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + mps_dir + "RR_s-mu_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-1000_mike.dat", 'r') as f:
            RR_mps_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')

        with open(base_data_path + wp_dir + "DD_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-1000_aemulus.dat", 'r') as f:
            DD_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DD_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-1000_aemulus.dat", 'r') as f:
            DD_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-1000_aemulus.dat", 'r') as f:
            DR_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "DR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-1000_aemulus.dat", 'r') as f:
            DR_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "RR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_NGC_mock-1000_aemulus.dat", 'r') as f:
            RR_wp_ngc = np.loadtxt(f, dtype=np.float, delimiter='\t')
        with open(base_data_path + wp_dir + "RR_rp-pi_log__eBOSS_ezmock_no-sel_LRG_SGC_mock-1000_aemulus.dat", 'r') as f:
            RR_wp_sgc = np.loadtxt(f, dtype=np.float, delimiter='\t')

    if i==0:
        print("Finished loading pair counts", dt.datetime.now() - start_time)

    if i==0:
        print("Started rebinning", dt.datetime.now() - start_time)

    DD_mps_ngc = rebin_pair_counts(DD_mps_ngc, 20, original_mps_binning)
    DD_mps_sgc = rebin_pair_counts(DD_mps_sgc, 20, original_mps_binning)
    DR_mps_ngc = rebin_pair_counts(DR_mps_ngc, 20, original_mps_binning)
    DR_mps_sgc = rebin_pair_counts(DR_mps_sgc, 20, original_mps_binning)
    RR_mps_ngc = rebin_pair_counts(RR_mps_ngc, 20, original_mps_binning)
    RR_mps_sgc = rebin_pair_counts(RR_mps_sgc, 20, original_mps_binning)

    DD_wp_ngc = rebin_pair_counts(DD_wp_ngc, 20, original_wp_binning)
    DD_wp_sgc = rebin_pair_counts(DD_wp_sgc, 20, original_wp_binning)
    DR_wp_ngc = rebin_pair_counts(DR_wp_ngc, 20, original_wp_binning)
    DR_wp_sgc = rebin_pair_counts(DR_wp_sgc, 20, original_wp_binning)
    RR_wp_ngc = rebin_pair_counts(RR_wp_ngc, 20, original_wp_binning)
    RR_wp_sgc = rebin_pair_counts(RR_wp_sgc, 20, original_wp_binning)

    if i==0:
        print("Finished rebinning", dt.datetime.now() - start_time)


    if i==0:
        print("Starting NGC/SGC combination", dt.datetime.now() - start_time)

    pair_counts_reduced = np.zeros((9,6))
    for i in range(9):
        pair_counts_reduced[i,0] += np.sum(DD_wp_ngc[i*150: i*150+80, -1])
        pair_counts_reduced[i,1] += np.sum(DR_wp_ngc[i*150: i*150+80, -1])
        pair_counts_reduced[i,2] += np.sum(RR_wp_ngc[i*150: i*150+80, -1])
        pair_counts_reduced[i,3] += np.sum(DD_wp_sgc[i*150: i*150+80, -1])
        pair_counts_reduced[i,4] += np.sum(DR_wp_sgc[i*150: i*150+80, -1])
        pair_counts_reduced[i,5] += np.sum(RR_wp_sgc[i*150: i*150+80, -1])

    print(pair_counts_reduced)

    DD_mps_ngc[:,2] += DD_mps_sgc[:,2]
    DR_mps_ngc[:,2] += DR_mps_sgc[:,2]
    RR_mps_ngc[:,2] += RR_mps_sgc[:,2]

    DD_wp_ngc[:,2] += DD_wp_sgc[:,2]
    DR_wp_ngc[:,2] += DR_wp_sgc[:,2]
    RR_wp_ngc[:,2] += RR_wp_sgc[:,2]

    if i==0:
        print("Finished NGC/SGC combination", dt.datetime.now() - start_time)


    if i==0:
        print("Starting correlation functions", dt.datetime.now() - start_time)

    mps_corr_func = landy_szalay(DD_mps_ngc, DR_mps_ngc, RR_mps_ngc)
    wp_corr_func = landy_szalay(DD_wp_ngc, DR_wp_ngc, RR_wp_ngc)

    mono_corr_func = corr_func_multipole_calc(0, mps_corr_func, mu_bins=mps_binning['N_bin'][1])
    quad_corr_func = corr_func_multipole_calc(2, mps_corr_func, mu_bins=mps_binning['N_bin'][1])
    proj_corr_func = corr_func_projected_calc(wp_corr_func, dr_pi=wp_step_sizes[1], r_pi_max=wp_binning['max'][1])

    if i==0:
        print("Finished correlation functions", dt.datetime.now() - start_time)

    data_vector = np.concatenate((mono_corr_func, quad_corr_func, proj_corr_func))
    np.savetxt(base_output_path + "%i_test.dat" % i, data_vector)

    print("Finished mock %i" % i, dt.datetime.now() - start_time)
