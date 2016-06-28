#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

void LinearGillespie_indiv(umat &nu, mat &a_mat, uvec &x0, vec &time_vec, umat& simulation_data)
{
  // Declare variables -----------------------------------------------------------
  // -----------------------------------------------------------------------------
//  umat        simulation_data=zeros<umat>(x0.n_elem, time_vec.n_elem);
  simulation_data.zeros();
  uvec                                                 xt;
  // .. General sampler parameters
  double                  cummul_a, ra0, r, a0, sim_time;
  uvec                                      indx_rates_a;
  vec   time_vec_aux(time_vec.n_elem+1), a_coeff, a_vec;
  int         t_i, next_reac;    


  // .. Set random Seed  
  //Rcpp::RNGScope scope;
//  GetRNGstate();

  indx_rates_a = conv_to<uvec>::from(a_mat.col(1));
  a_coeff      =                      a_mat.col(0);
  indx_rates_a = indx_rates_a-1; // start indices at 0
  // Init sim var
  xt       = x0;
  sim_time = 0;
  t_i      = 0;

  // Add dummy value at the end to avoid crash
  time_vec_aux.subvec(0,time_vec.n_elem-1) = time_vec;
  time_vec_aux(time_vec.n_elem) = pow(10,10);
  // Start the simulation
  do{   
      // Recalculate reaction rates
      a_vec = a_coeff%xt.elem(indx_rates_a);
  
      // Total reaction rate
      a0 = sum(a_vec);
      
      // If no more reactions possible fill simulation matrix and break
      if(a0 == 0){
        for(int i = t_i; i<simulation_data.n_cols; i++) simulation_data.col(i) = xt;
        break;      
      } 
  
      // Update sim time
      r         =  unif_rand();
      sim_time +=   -log(r)/a0;
  
      // Fill in necessary time points if appropriate (with old values)
      while(time_vec_aux(t_i)<sim_time & t_i<time_vec.n_elem){
        // Store result
        simulation_data.col(t_i) = xt;
        t_i++;
      }
  
      // Update concentrations
      r         =  unif_rand();
      ra0       =  r * a0;
      cummul_a = 0;
      for (next_reac = 0; next_reac < a_vec.n_elem-1; next_reac++) {
         cummul_a += a_vec(next_reac);
         if (ra0 <= cummul_a) break;
      }      
      xt += nu.col(next_reac); 
  }while(t_i<time_vec.n_elem);
  
//  PutRNGstate();
//  return(simulation_data);
}

// [[Rcpp::export]]
umat LinearGillespie(int numSim, umat nu, mat a_mat, uvec x0, vec time_vec)
{
  int sim_i;
  umat simulation_data_all=zeros<umat>(x0.n_elem, time_vec.n_elem);
  umat simulation_data=zeros<umat>(x0.n_elem, time_vec.n_elem);
  // .. Get random Seed  
  GetRNGstate();
  for(sim_i = 0; sim_i<numSim; sim_i++)
  {
    LinearGillespie_indiv(nu, a_mat, x0, time_vec, simulation_data);
    simulation_data_all = simulation_data_all + simulation_data;
  }
  // .. Return random seed    
  PutRNGstate();
  return(simulation_data_all);
}


