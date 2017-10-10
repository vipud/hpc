//C version of KMC
//header
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>



char* tmp;


//Structs

struct Traj_stats{
  char* species_avgs_out_flname;
  char* SA_out_flname;
  double** spec_profiles_averages;
  double** traj_deriv_avgs;
  double*** sensitivities;
  double*** microscale_sensitivity_contributions;
};

struct KMC_traj{
  char* species_out_flname;
  char* traj_deriv_out_flname;
  int* N;
  double t_trunc;
  double t;
  double t_prev;
  double dt;
  int rxn_to_file_ind;
  int ind_rec;
  double* props;

  double* W;
  double** prop_ders;
  double* prop_ders_sum;

  double** spec_profile;
  double** traj_deriv_profile;
};

struct KMC_traj_TTS{
  char* species_out_flname;
  char* traj_deriv_out_flname;

  int* N;
  double t_trunc;
  double t;
  double t_prev;
  double dt;
  int rxn_to_file_ind;
  int ind_rec;
  double* props;

  double* W;
  double** prop_ders;
  double* prop_ders_sum;

  double** spec_profile;
  double** traj_deriv_profile;

  double* N_micro_avg;
  double** micro_scale_sens;

  //Microscale vectors
  double* dEdth;
  double* dEdth_avg;
  double* rev_prop_ders;
  int** N_rec;
  double* q_cum;
  int* N_candidate;
  double* micro_props;
  double** micro_prop_ders;
  double** prop_ders_direct;           // Direct averaging of derviativeson the microscale
  double** prop_ders_indirect;         // Indirect averaging on the microscale
  double* slow_props;
  double* slow_props_cum;

};




//variables
//input variables
char *** readIn;
int N_record = 101;
int N_traj = 1000;
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

//functions


//read file
void file_reader(char* fileName){
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
  // printf("%s\n", spec_names[2]);

  N_0 = (int*) malloc(sizeof(int)*n_specs);
  for(int i = 0; i<n_specs; i++){
    N_0[i] = atoi(readIn[4][i+1]);
  }

  // rxn_stoich = (int*)malloc(sizeof(int)*n_specs);
  stoich_mat = (int**)malloc(sizeof(int*)*n_rxns);

  for(int i = 0; i<n_rxns; i++){
    stoich_mat[i] = (int*)malloc(sizeof(int)*n_specs);
    for(int j = 0; j<n_specs; j++){
      stoich_mat[i][j] =atoi(readIn[5][1+n_specs*i +j]);
    }
  }
  // for (size_t i = 0; i < 3; i++) {
  //   for (size_t j = 0; j < 3; j++) {
  //     printf("%d\n",stoich_mat[i][j] );
  //   }
  // }

  rate_const = (double*)malloc(sizeof(double)*n_rxns);
  for(int i=0; i<n_rxns; i++){
    rate_const[i] = atof(readIn[6][i+1]);
  }

  param_names =(char**) malloc(sizeof(char*) * n_params);
  for(int i=0; i<n_params; i++){
    param_names[i]=readIn[7][i+1];
  }

  // printf("%s\n", param_names[2]);

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


}

void initStats(struct Traj_stats *sim){
  //TODO
  sim->species_avgs_out_flname = "species_avgs_out.txt";
  sim->SA_out_flname = "sensitivities_out.txt";


  sim->spec_profiles_averages =(double**) malloc(sizeof(double*)*N_record);
  for(int i =0; i< N_record; i++){
    sim->spec_profiles_averages[i] = (double*)malloc(sizeof(double)*n_specs);
  }

  sim->traj_deriv_avgs =(double**) malloc(sizeof(double*)*N_record);
  for(int i =0; i< N_record; i++){
    sim->traj_deriv_avgs[i] = (double*)malloc(sizeof(double)*n_specs);
  }

  sim->sensitivities =(double***) malloc(sizeof(double**)*N_record);
  for(int i =0; i< N_record; i++){
    sim->sensitivities[i] = (double**)malloc(sizeof(double*)*n_specs);
    for(int j =0; j< n_specs; j++){
      sim->sensitivities[i][j] = (double*) malloc(sizeof(double)*n_params);
    }
  }


  for (int i = 0; i < N_record; ++i){
      for(int j = 0; j < n_specs; j++){
          sim->spec_profiles_averages[i][j] = 0;
      }
  }

  for (int i = 0; i < N_record; ++i){
      for(int j = 0; j < n_params; j++){
          sim->traj_deriv_avgs[i][j] = 0;
      }
  }

  for (int i = 0; i < N_record; ++i){
      for(int j = 0; j < n_specs; j ++){
          for(int k = 0; k < n_params; k++){
              sim->sensitivities[i][j][k] = 0;
          }
      }
  }
}

void finalize_stats(struct Traj_stats *sim){
  for (int i = 0; i < N_record; ++i){

      for(int j = 0; j < n_specs; j++){
          sim->spec_profiles_averages[i][j] = sim->spec_profiles_averages[i][j] / N_traj;

          for(int k = 0; k < n_params; k++){
              sim->sensitivities[i][j][k] = sim->sensitivities[i][j][k] / N_traj;
          }
      }

      for(int j = 0; j < n_params; j++){
          sim->traj_deriv_avgs[i][j] = sim->traj_deriv_avgs[i][j] / N_traj;
      }
  }

  for (int i = 0; i < N_record; ++i){

      for(int j = 0; j < n_specs; j++){

          for(int k = 0; k < n_params; k++){

              //cout << endl;
              //cout << sensitivities[i][j][k] << endl;
              sim->sensitivities[i][j][k] -= sim->spec_profiles_averages[i][j] * sim->traj_deriv_avgs[i][k];
              sim->sensitivities[i][j][k] = sim->sensitivities[i][j][k] * N_traj / (N_traj - 1);
              //cout << sensitivities[i][j][k] << endl;
          }
      }

  }
}

// void initRandomNumbers(int** randomNumbers){
//   //TODO
// }
//
// void initTraj(struct KMC_traj* trajs){
//   //TODO
// }
//
// void initTraj(struct KMC_traj_TTS* trajs_TTS){
//   //TODO
// }
//
// void simulate(struct KMC_traj traj, int random){
//   //TODO
// }
//
// void simulate(struct KMC_traj_TTS traj_TTS, int random){
//   //TODO
// }
//
// void run_simulations(struct KMC_traj* trajs){
//   //TODO
// }
//
// void run_simulations(struct KMC_traj_TTS * trajs_TTS){
//   //TODO
// }

int main(){

  // struct Traj_stats sim;
  // struct KMC_traj * trajs;
  // struct KMC_traj_TTS * trajs_TTS;
  // int** randomNumbers;
  //
  // file_reader("network.in");
  // initStats(sim);
  // initRandomNumbers(randomNumbers);
  // if(two_time_scale){
  //   initTraj(trajs);
  //   run_simulations(trajs_TTS);
  // }else{
  //   initTraj(trajs_TTS);
  //   run_simulations(trajs);
  // }

  file_reader("input.csv");
  struct Traj_stats* sim;
  initStats(sim);





  free(tmp);
  return 0;

}
