// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


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
arma::umat LinearGillespie(int numSim, arma::umat nu, arma::mat a_mat, arma::uvec x0, arma::vec time_vec)
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



// [[Rcpp::export]]
arma::mat Crypt_drift_c(double alpha, double beta, int Ns, arma::vec time_intervals, bool persisting, int splitNum = 0)
{
  double t;
  rowvec theory_timeseries_row0;

  vec    f_nm(Ns-1), g_nm(Ns-1), aux1(Ns-1), pow_gnm(Ns-1);
  double pi = 3.1415926535897932384626433832795028841971693993751058;

  vec              m  = linspace<vec>(1, Ns-1, Ns-1);
  double          cte = 2./Ns;
  double  al_bet_sqrt   = sqrt(alpha*beta);
  double  beta_ov_alpha = beta/alpha;
  double  aux_al_bet    = alpha + beta - 2. * al_bet_sqrt;

  mat theory_timeseries1 = zeros<mat>(Ns+1, time_intervals.n_elem);

  g_nm = 4*al_bet_sqrt*square(sin(0.5*pi*m/Ns)) + aux_al_bet;
  //f_nm = pow(cos(0.5*pi*m/Ns), 2);
  f_nm = pow(sin(pi*m/Ns), 2);
  f_nm = f_nm/g_nm;

  // .. For n = 0
  for(int i = 0; i<time_intervals.n_elem; i++)
  {
    t = time_intervals(i);
//    theory_timeseries1(0, i) = as_scalar(alpha*pow(beta_ov_alpha,0.5)*cte*sum(f_nm%(1-exp(-g_nm*t))));
    theory_timeseries1(0, i) = as_scalar(alpha*cte*sum(f_nm%(1-exp(-g_nm*t))));
  }  

  // .. For  1 <= n <= Ns-1
  for(int n = 1; n<Ns; n++)
  {
    f_nm = sin(pi*m/Ns)%sin(pi*m*n/Ns);
    for(int i = 0; i<time_intervals.n_elem; i++)
    {
      t = time_intervals(i);
      theory_timeseries1(n, i) = as_scalar(pow(beta_ov_alpha, (n-1)/2.)*cte*sum(f_nm%exp(-g_nm*t)));
      //theory_timeseries1(n, i) = as_scalar(pow(beta_ov_alpha, n/2.)*cte*sum(f_nm%exp(-g_nm*t)));
    }  
  }
  
  // .. For n = Ns
//  for(int nn = 0; nn<aux1.n_elem; nn++) aux1(nn) = pow(-1., nn);
//  f_nm = aux1%pow(cos(0.5*pi*m/Ns), 2);
  f_nm = sin(pi*m/Ns)%sin(pi*m*(Ns-1)/Ns);
  f_nm = f_nm/g_nm;

  for(int i = 0; i<time_intervals.n_elem; i++)
  {
    t = time_intervals(i);
    theory_timeseries1(Ns, i) = as_scalar(beta*pow(beta_ov_alpha, (Ns-2)/2.)*cte*sum(f_nm%(1-exp(-g_nm*t))));
//    theory_timeseries1(Ns, i) = as_scalar(beta*pow(beta_ov_alpha, (Ns-1)/2.)*cte*sum(f_nm%(1-exp(-g_nm*t))));
  }  

  // Store first row
  if(persisting)
  {
    theory_timeseries1.row(0) = 1-theory_timeseries1.row(0);
    for(int n = 1; n<(Ns+1); n++)
    {    
      theory_timeseries1.row(n) = theory_timeseries1.row(n)/theory_timeseries1.row(0);
    }
    theory_timeseries1.shed_row(0);
  }
  else // if not persisting remove first row to ease collapsing
  {
    theory_timeseries_row0 = theory_timeseries1.row(0);
    theory_timeseries1.shed_row(0);
  }

  if(splitNum!=0)
  {
    // .. Declare some vars
    int           old_row_num = theory_timeseries1.n_rows;
    int                          new_row_num = splitNum;
    int                                            indx_i;
    double                                 scaling_factor;
    mat new_x       = zeros<mat>(new_row_num, theory_timeseries1.n_cols);
    vec aux_vec     = zeros<vec>(new_row_num*old_row_num);
    vec vals_i      = zeros<vec>(old_row_num);
    // .. Change so that full is only for the last bin
    if(new_row_num > old_row_num)
    {
      new_x.row(new_row_num-1) = theory_timeseries1.row(old_row_num-1);
      new_row_num--;
      old_row_num--;
    }
    // .. Change !
    for(int j = 0; j < old_row_num; j++) aux_vec.subvec(j*new_row_num, (j+1)*new_row_num - 1) = j*ones<vec>(new_row_num);
    for(int i = 0; i < new_row_num; i++)
    {
      vals_i       = aux_vec.subvec(i*old_row_num,(i+1)*old_row_num-1);
      vec indx_use = unique(vals_i);
      for(int indx_i_loop = 0; indx_i_loop < indx_use.n_elem; indx_i_loop++)
      {
        indx_i         = indx_use(indx_i_loop);
      	scaling_factor = sum(vals_i == indx_i)/(double) new_row_num;
      	new_x.row(i)   = new_x.row(i) + scaling_factor*theory_timeseries1.row(indx_i);
      }
    }
    theory_timeseries1 = new_x;    
  }
  // If perisiting is false add first row N = 0
  if(!persisting)
  {
    theory_timeseries1.insert_rows( 0, theory_timeseries_row0);
  }
  return abs(theory_timeseries1);
}

