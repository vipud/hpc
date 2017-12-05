#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <openacc.h>

char* species_out_flname;
char* traj_deriv_out_flname;
int** N;
double* t_trunc;
double* t;
double* t_prev;
double* dt;
int* rxn_to_fire_ind;
int* ind_rec;
double** props;

double** W;
double*** prop_ders;
double** prop_ders_sum;

double*** spec_profile;
double*** traj_deriv_profile;

struct Traj_stats{
  char* species_avgs_out_flname;
  char* SA_out_flname;
  double** spec_profiles_averages;
  double** traj_deriv_avgs;
  double*** sensitivities;
  double*** microscale_sensitivity_contributions;
};


int main(){
  int countiii = 0;

  //variables
  char* tmp;
  //input variables
  char *** readIn;
  int N_record;
  int N_traj;
  bool write_traj_files = false;
  int rand_seed = 1;
  bool two_time_scale = false;
  int n_fast_pairs = 0;
  int n_fast_rxns = 0;
  int n_slow_rxns = 0;
  int n_micro_steps = 1000;
  int n_rxns;
  int n_specs;
  int n_params;
  int t_final;
  double eps;
  char** spec_names;
  char** param_names;
  int* N_0;
  int** stoich_mat;
  double* rate_const;
  double* t_rec;
  int** fast_pairs;
  int* fast_rxns;
  int* slow_rxns;

  int fwd_rxn_ind;
  int rev_rxn_ind;

  int slow_ind = 0;

  char* fileName = "input.csv";

  
  struct Traj_stats sim;

  //file_reader
  char line[1024];
  FILE *stream;
  stream = fopen(fileName,"r");

  readIn = (char ***)malloc(sizeof(char**) * 13);
  for(int i=0; i<13; i++){
    readIn[i] = (char**)malloc(sizeof(char*)*10);
    for(int j = 0; j<10; j++){
      readIn[i][j] = (char*)malloc(sizeof(char*)*100);
    }
  }

  int count=0;
  while (fgets(line,1024,stream)){
    tmp = strdup(line);
    char* token;
    token = strtok(tmp,",");
    int j = 0;
    while(token!= NULL){
      readIn[count][j] = token;
      token = strtok(NULL,",");
      j++;
    }

    count++;
  }
  fclose(stream);

  rand_seed = atoi(readIn[0][1]);
  n_specs = atoi(readIn[1][1]);

  n_rxns = atoi(readIn[2][1]);
  n_params = n_rxns;

  spec_names = (char**)malloc(sizeof(char*)*n_specs);
  for(int i = 0; i<n_specs; i++){
    spec_names[i] = readIn[3][i+1];
  }

  N_0 = (int*) malloc(sizeof(int)*n_specs);
  for(int i = 0; i<n_specs; i++){
    N_0[i] = atoi(readIn[4][i+1]);
  }
  stoich_mat = (int**)malloc(sizeof(int*)*n_rxns);

  for(int i = 0; i<n_rxns; i++){
    stoich_mat[i] = (int*)malloc(sizeof(int)*n_specs);
    for(int j = 0; j<n_specs; j++){
      stoich_mat[i][j] =atoi(readIn[5][1+n_specs*i +j]);
    }
  }

  rate_const = (double*)malloc(sizeof(double)*n_rxns);
  for(int i=0; i<n_rxns; i++){
    rate_const[i] = atof(readIn[6][i+1]);
  }

  param_names =(char**) malloc(sizeof(char*) * n_params);
  for(int i=0; i<n_params; i++){
    param_names[i]=readIn[7][i+1];
  }
  t_final = atof(readIn[8][1]);

  N_traj = atoi(readIn[9][1]);

  N_record = atoi(readIn[10][1])+1;

  n_fast_pairs = atoi(readIn[11][1]);
  n_fast_rxns = n_fast_pairs*2;
  n_slow_rxns = n_rxns - n_fast_rxns;

  fast_rxns = (int*)malloc(sizeof(int) * n_fast_rxns);
  slow_rxns = (int*)malloc(sizeof(int) * n_slow_rxns);
  fast_pairs = (int**)malloc(sizeof(int*) *n_fast_rxns);
  for(int i=0; i<n_fast_rxns; i++){
    fast_pairs[i]=(int*)malloc(sizeof(int)*2);
  }

  for(int i =0; i<n_fast_pairs; i++){
    fwd_rxn_ind = atoi(readIn[12][1]);
    fast_pairs[2*i][0] = fwd_rxn_ind;
    fast_pairs[2*i +1][1] = fwd_rxn_ind;
    rev_rxn_ind = atoi(readIn[12][2]);
    fast_pairs[2*i][1]=rev_rxn_ind;
    fast_pairs[2*i +1][0] = rev_rxn_ind;

    fast_rxns[2*i]=rev_rxn_ind;
    fast_rxns[2*i +1] = rev_rxn_ind;
  }

  for(int i = 0; i<n_rxns; i++){
    for(int j = 0; j <n_fast_rxns; j++){
      if(fast_rxns[j] == i){
        slow_rxns[slow_ind] = i;
        slow_ind ++;
      }
    }
  }

  t_rec = (double*) malloc(sizeof(double)*N_record);
  for(int i = 0; i< N_record; i++){
    t_rec[i] = (double)i / (N_record -1) * t_final;
  }



  //initStats
  sim.species_avgs_out_flname = "species_avgs_out.txt";
  sim.SA_out_flname = "sensitivities_out.txt";


  sim.spec_profiles_averages =(double**) malloc(sizeof(double*)*N_record);
  for(int i =0; i< N_record; i++){
    sim.spec_profiles_averages[i] = (double*)malloc(sizeof(double)*n_specs);
  }

  sim.traj_deriv_avgs =(double**) malloc(sizeof(double*)*N_record);
  for(int i =0; i< N_record; i++){
    sim.traj_deriv_avgs[i] = (double*)malloc(sizeof(double)*n_params);
  }

  sim.sensitivities =(double***) malloc(sizeof(double**)*N_record);
  for(int i =0; i< N_record; i++){
    sim.sensitivities[i] = (double**)malloc(sizeof(double*)*n_specs);
    for(int j =0; j< n_specs; j++){
      sim.sensitivities[i][j] = (double*) malloc(sizeof(double)*n_params);
    }
  }


  for (int i = 0; i < N_record; i++){
    for(int j = 0; j < n_specs; j++){
      sim.spec_profiles_averages[i][j] = 0;
    }
  }

  for (int i = 0; i < N_record; i++){
    for(int j = 0; j < n_params; j++){
      sim.traj_deriv_avgs[i][j] = 0;
    }
  }

  for (int i = 0; i < N_record; i++){
    for(int j = 0; j < n_specs; j ++){
      for(int k = 0; k < n_params; k++){
        sim.sensitivities[i][j][k] = 0;
      }
    }
  }
//////////////////////////////////////////////////////////////////////////
  double randomNumbers[N_traj][2][2000];
  for(int i =0; i< N_traj; i++){
    srand(rand_seed + i);
    for(int j = 0; j<2; j++){
      for(int k =0; k< 500; k++){
        randomNumbers[i][j][k] = ((double) rand() / (double)(RAND_MAX));
      }
    }
  }
  species_out_flname = "species_out.txt";
  traj_deriv_out_flname = "traj_deriv_out.txt";
  N = (int**)malloc(sizeof(int*)*N_traj);
  t_trunc = (double*)malloc(sizeof(double)*N_traj);
  t = (double*)malloc(sizeof(double)*N_traj);
  t_prev = (double*)malloc(sizeof(double)*N_traj);
  dt = (double*)malloc(sizeof(double)*N_traj);
  rxn_to_fire_ind = (int*)malloc(sizeof(int)*N_traj);
  ind_rec = (int*)malloc(sizeof(int)*N_traj);

  props = (double**)malloc(sizeof(double*)*N_traj);
  prop_ders_sum = (double**)malloc(sizeof(double*)*N_traj);
  prop_ders = (double***)malloc(sizeof(double**)*N_traj);
  W = (double**)malloc(sizeof(double*)*N_traj);
  spec_profile = (double***)malloc(sizeof(double**)*N_traj);
  traj_deriv_profile = (double***)malloc(sizeof(double**)*N_traj);

  for (int i = 0; i < N_traj; i++){
    t[i] = 0;
    t_prev[i] = 0;
    ind_rec[i] = 0;
    N[i] = (int*)malloc(sizeof(int)*n_specs);
    for (int j = 0; j < n_specs; j++){
      N[i][j] = N_0[j];
    }
    
    props[i] = (double*)malloc(sizeof(double)* n_rxns);
    prop_ders_sum[i] = (double*)malloc(sizeof(double)*n_rxns);
    prop_ders[i] = (double**)malloc(sizeof(double*)*n_rxns);
    for(int j =0; j< n_rxns; j++){
      prop_ders[i][j] = (double*) malloc(sizeof(double)*n_params);
    }
    W[i] = (double*)malloc(sizeof(double)*n_params);
    for (int j = 0; j < n_params; j++){
    W[i][j] = 0;  
    }
    spec_profile[i] =(double**) malloc(sizeof(double*)*N_record);
    for(int j=0; j<N_record; j++){
      spec_profile[i][j] = (double*) malloc(sizeof(double)*n_specs);
      for(int k =0; k<n_specs;k++){
        spec_profile[i][j][k] = 0;
      }
    }
    traj_deriv_profile[i] = (double**)malloc(sizeof(double*)*N_record);
    for (int j = 0; j < N_record; j++){
      traj_deriv_profile[i][j] = (double*)malloc(sizeof(double)*N_record);
      for (int k = 0; k < n_params; k++){
        traj_deriv_profile[i][j][k] = 0;
      }
    }
  }

#pragma acc data copyin (N[:N_traj][:n_specs], t_trunc[:N_traj], t[:N_traj], t_prev[:N_traj], dt[:N_traj], rxn_to_fire_ind[:N_traj], ind_rec[:N_traj], props[:N_traj][:n_rxns], W[:N_traj][:n_params], prop_ders[:N_traj][:n_rxns][:n_params], prop_ders_sum[:N_traj][:n_params], spec_profile[:N_traj][:N_record][:n_specs], traj_deriv_profile[:N_traj][:N_record][:n_params], stoich_mat[:n_rxns][:n_specs], rate_const[:n_rxns], randomNumbers[:N_traj][:2][:2000], t_rec[:N_record]) copyout(spec_profile[:N_traj][:N_record][:n_specs], traj_deriv_profile[:N_traj][:N_record][:n_params])
{
#pragma acc parallel loop independent
  for (int x = 0; x < N_traj; x++){
    
    int r = 0;
    double r_rxn_choose;
    double r_timestep;
    double asum;
    double* prop_cum = (double*)malloc(sizeof(double) * n_rxns);
    for (; t[x] < t_final;){
   
      r_rxn_choose = randomNumbers[x][0][r];
      r_timestep = randomNumbers[x][1][r];
      r+=1;
      
      for (int i = 0; i < n_rxns; i++){
        
        props[x][i] = rate_const[i];
        for (int j = 0; j < n_specs; j++){
        
          if (stoich_mat[i][j] < 0){
            props[x][i] *= pow(N[x][j], -stoich_mat[i][j]);
          }
        }
        for (int k = 0; k < n_params; k++){
          if (i == k){         
            prop_ders[x][i][k] = 1.0;
            for (int l = 0; l < n_specs; l++){
              if (stoich_mat[i][l] < 0){
                prop_ders[x][i][k] *= pow(N[x][l], -stoich_mat[i][l]);
              }
            }
          }
          else{
            prop_ders[x][i][k] = 0;
          }
        }
      }
      for (int i = 0; i < n_params; i++){
        double numm = 0;
        for (int j = 0; j < n_rxns; j++){
          numm += prop_ders[x][j][i];
        }
        prop_ders_sum[x][i] = numm;
      }
      asum = 0;
      for (int i = 0; i < n_rxns; i++){
        asum += props[x][i];
      }

      if (asum == 0){ //must be changed for parallel structure
        break;
      }
      
      for (int i = 0; i < n_rxns; i++){
        double numm = 0;
        for (int j = 0; j <= i; j++){
          numm  += props[x][j]/asum;
        }
        prop_cum[i] = numm;
      }

      rxn_to_fire_ind[x] = 0;
      if (r_rxn_choose == 1){
        for (int i = 0; i < n_rxns; i++){
          if (props[x][i] > 0){
            rxn_to_fire_ind[x] = i;
          }
        }
      }
      else{
        for (; prop_cum[rxn_to_fire_ind[x]] < r_rxn_choose; rxn_to_fire_ind[x]++){}

      }

      dt[x] = log(1/r_timestep)/asum;
      for (; t[x] >= t_rec[ind_rec[x]]; ind_rec[x]++){
        t_trunc[x] = t_rec[ind_rec[x]] - t_prev[x];
        for (int i = 0; i < n_specs; i++){
          spec_profile[x][ind_rec[x]][i] = (double) N[x][i];
        }
        for (int i = 0; i < n_params; i++){
          traj_deriv_profile[x][ind_rec[x]][i] = W[x][i] - prop_ders_sum[x][i] * t_trunc[x];
        }
      }
      for (int i = 0; i < n_specs; i++){
        N[x][i] += stoich_mat[rxn_to_fire_ind[x]][i];
      }

      t_prev[x] = t[x];
      t[x] = t[x] + dt[x];

      for (int i = 0; i <n_params; i++){
        W[x][i] += prop_ders[x][rxn_to_fire_ind[x]][i] / props[x][rxn_to_fire_ind[x]];
        W[x][i] -= prop_ders_sum[x][i] * dt[x];
      }
    }
    
    for (; ind_rec[x] < N_record; ind_rec[x]++){
      t_trunc[x] = t_rec[ind_rec[x]] - t_prev[x];
#pragma acc loop seq
      for (int i = 0; i < n_specs; i++){
        spec_profile[x][ind_rec[x]][i] = (double) N[x][i];
      }
      for (int i = 0; i < n_params; i++){
        traj_deriv_profile[x][ind_rec[x]][i] = W[x][i] - prop_ders_sum[x][i] * t_trunc[x];
      }
    }
   // printf("trial %i \n ", countiii++);
  }
}

/////////////////////////////////////////////////////////////////////////////////
  for(int i = 0; i<N_traj; i++){
    for(int j = 0; j< N_record; j++){
      for(int k = 0; k < n_specs; k++){
        sim.spec_profiles_averages[j][k] += spec_profile[i][j][k];
        for(int l =0; l< n_params; l++){
          sim.sensitivities[j][k][l] += spec_profile[i][j][k] * traj_deriv_profile[i][j][l] + 0;
        }
      }

      for(int k =0; k < n_params; k++){
        sim.traj_deriv_avgs[j][k] += traj_deriv_profile[i][j][k];
      }
    }
  }
  for (int i = 0; i < N_record; i++){

      for(int j = 0; j < n_specs; j++){
          sim.spec_profiles_averages[i][j] = sim.spec_profiles_averages[i][j] / N_traj;
          for(int k = 0; k < n_params; k++){
              sim.sensitivities[i][j][k] = sim.sensitivities[i][j][k] / N_traj;
          }
      }
      for(int j = 0; j < n_params; j++){
          sim.traj_deriv_avgs[i][j] = sim.traj_deriv_avgs[i][j] / N_traj;
      }
  }
  for (int i = 0; i < N_record; i++){
      for(int j = 0; j < n_specs; j++){
          for(int k = 0; k < n_params; k++){
              sim.sensitivities[i][j][k] -= sim.spec_profiles_averages[i][j] * sim.traj_deriv_avgs[i][k];
              sim.sensitivities[i][j][k] = sim.sensitivities[i][j][k] * N_traj / (N_traj - 1);
          }
      }

  }

  for(int i = 0; i< n_specs;i++){
    printf("%s\t", spec_names[i]);
  }
  printf("\n");

  for (int i = 0; i < N_record; i++) {
    /* code */
    printf("%f\t",t_rec[i] );
    for (int j = 0; j < n_specs; j++) {
      /* code */
      printf("%f\t", sim.spec_profiles_averages[i][j] );
    }
    printf("\n" );
  }
  return 0;

}
