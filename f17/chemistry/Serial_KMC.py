import random
import math

species_out_flname = ''
traj_deriv_out_flname = ''
N = 0
t_trunc = 0.0
t = 0.0
t_prev = 0.0
dt = 0.0
rxn_to_fire_ind = 0
ind_rec = 0
props = 0.0

W = 0.0
prop_ders = 0.0
prop_ders_sum = 0.0
spec_profile = 0.0
traj_deriv_profile = 0.0

class Traj_stats(object):
    def __init__(self):
        self.species_avg_out_flname = ''
        self.SA_out_flname = ''
        self.spec_proiles_averages = 0.0
        self.traj_deriv_avgs = 0.0
        self.sensitivities = 0.0
        self.microscale_sensitivity_contributions = 0.0

countiii = 0
temp = ''
readIn = ''
N_record = 0
N_traj = 0
write_traj_files = false
rand_seed = 1
two_time_scale = false
n_fast_pairs = 0
n_fast_rxns = 0
n_slow_rxns = 0
n_micro_steps = 1000
n_rxns = 0
n_specs = 0
n_params = 0
t_final = 0
eps = 0.0
spec_names = ''
param_names = ''
N_0 = 0
stoich_mat = 0
rate_const = 0.0
t_rec = 0.0
fast_pairs = 0
fast_rxns = 0
slow_rxns = 0
fwd_rxn_ind = 0
rev_rxn_ind = 0
slow_ind = 0
fileName = "input.csv"
sim = Traj_stats()


# READ IN AND SET DATA - NOT CLEAR ON WHAT GOES WITH WHAT
# https://github.com/vipud/hpc/blob/master/f17/chemistry/C/KMC-acc.c
# 106-208
sim.species_avg_out_flname = "species_avgs_out.txt"
sim.SA_out_flname = "sensitivities_out.txt"
for(i in range(0, N_record)):
    for(j in range(0, n_specs)):
        sim.spec_profile_averages[i][j] = 0
for(i in range(0, N_record)):
    for(j in range(0, n_params)):
        sim.traj_deriv_avgs[i][j] = 0
for(i in range(0, N_record)):
    for(j in range(0, n_specs)):
        for(k in range(0, n_params)):
            sim.sensitivities[i][j][k] = 0
randomNumbers = [N_traj][2][2000]
for(i in range(0, N_traj)):
    random.seed(rand_seed + i)
    for(j in range(0,2)):
        for(k in range(0,500)):
            randomNumbers[i][j][k] = ((double)random.random()/32767)
species_out_flname = "species_out.txt"
traj_deriv_out_flname = "traj_deriv_out.txt"

for(i in range(i, N_traj)):
    t[i] = 0
    t_prev[i] = 0
    ind_rec[i] = 0
    for(j in range(0, n_specs)):
        N[i][j] = N_0[j]

    for(j in range(0, n_params)):
        W[i][j] = 0
    for(j in range(0, N_record)):
        for(k in range(0, n_spec)):
            spec_profile[i][j][k] = 0
    for(j in range(0, N_record)):
        for(k in range(0, n_params)):
            traj_deriv_profile[i][j][k] = 0


#Logic
for(x in range (0, N_traj)):
    r = 0
    r_rxn_choose = 0.0
    r_timestep = 0.0
    asum = 0.0
    while(t[x] < t_final):
        r_rxn_choose = randomNumbers[x][0][r]
        r_timestep = randomNumbers[x][1][r]
        r+=1
        for(i in range(0, n_rxns)):
            props[x][i] = rate_const[i]
            for(j in range(0, n_specs)):
                if(stoich_mat[i][j] < 0):
                    props[x][i] *= math.pow(N[x][j], -stoich_mat[i][j])
            for(k in range(0, n_params)):
                if(i == k):
                    prop_ders[x][i][k] = 1.0
                    for(l in range(0, n_specs)):
                        if(stpoch_matrix[i][l] < 0):
                            prop_ders[x][i][k] *= math.pow(N[x][l],
                             -stoich_mat[i][l])
                else:
                    prop_ders[x][i][k] = 0
            for(i in range(0, n_params)):
                numm = 0.0
                for(j in range(0, n_rxns)):
                    numm += prop-ders[x][j][i]
                prop_ders_sum[x][i] = numm
            asum = 0
            for(i in range(i, n_rxns)):
                asum += props[x][i]
            if(asum == 0):
                break
            for(i in range(0, n_rxns))
                numm = 0
                for(j in range(0, i+1)):
                    numm += props[x][j]/asum
                prop_cum[i] = numm
            rxn_to_fire_ind[x] = 0
            if(r_rxn_choose == 1):
                for(i in range(0, n_rxns)):
                    if(props[x][i] > 0):
                        rxn_to_fire_ind[x] = i
            else:
                while(prop_cum[rxn_to_fire_ind[x]] < r_rxn_choose):
                    rxn_to_fire_ind[x] = rxn_to_fire_ind[x] + 1

            dt[x] = math.log(1/r_timestep)/asum
            while(t[x] >= t_rec[ind_rec[x]]):
                ind_rec[x] = ind_rec[x] + 1
                for(i in range(0, n_specs)):
                    spec_profile[x][ind_rec[x]][i] = (double)N[x][i]
                for(i in range(0, n_params)):
                    traj_deriv_profile[x][ind_rec[x]][i] = W[x][i] -
                                                prop_ders_sum[x][i] * t_trunc[x]
            for(i in range(0, n_specs)):
                N[x][i] += stoich_mat[rxn_to_fire_ind[x]][i]
            t_prev[x] = t[x]
            t[x] = t[x] + dt[x]

            for(i in range(0, n_params)):
                w[x][i] += prop_ders[x][rxn_to_fire_ind[x]][i]/props[x][rxn_to_fire_ind[x]]
                w[x][i] - prop_ders_sum[x][i] * dt[x]

        while(ind_rec[x] < N_record):
            ind_rec[x] = ind_rec[x] + 1
            for(i in range(0, n_specs)):
                spec_profile[x][ind_rec[x]][i] = (double)N[x][i]
            for(i in range(0, n_params)):
                traj_deriv_profile[x][ind_rec[x]][i] =
                                W[x][i] - prop_ders_sum[x][i] * t_trunc[x]

        for(i in range(0, N_traj)):
            for(j in range(0, N_record)):
                for(k in range(0, n_specs)):
                    sim.spec_profiles_averages[j][k] += spec_profile[i][j][k]
                    for(l in range(0, n_params)):
                        sim.sensitivities[j][k][l] += spec_profile[i][j][k] *
                                                traj_deriv_profile[i][j][l] + 0
                for(k in range(0, n_params)):
                    sim.traj_deriv_avgs[j][k] += traj_deriv_profile[i][j][k]

        for(i in range(0, N_record)):
            for(j in range(0, n_specs)):
                sim.spec_profiles_averages[i][j] =
                        sim.spec_profiles_averages[i][j] / N_traj
                for(k in range(0, n_params)):
                    sim.sensitivities[i][j][k] = sim.sensitivities[i][j][k]/N_traj
            for(j in range(0, n_params)):
                sim.traj_deriv_avgs[i][j] = sim.traj_deriv_avgs[i][j]/N_traj

        for(i in range(0, N_record)):
            for(j in range(0, n_specs)):
                for(k in range(0, n_params)):
                    sim.sensitivities[i][j][k] -= sim.spec_profiles_averages[i][j] * sim.traj_deriv_avgs[i][k];
                    sim.sensitivities[i][j][k] = sim.sensitivities[i][j][k] * N_traj / (N_traj - 1);

        for(i in range(0, n_specs)):
            print('%n', spec_ranges[i])
        print(" ")
        for(i in range(0, N_record)):
            print("%n", t_rec[i])
            for(j in range(0, n_specs)):
                print(%n, sim.spec_profiles_averages[i][j])
            print(" ")
