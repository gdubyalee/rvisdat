#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// From the vec locate the 1 (source element), -1 (displaced element) and 
// -2 (replaced element, may not exist)
void updateReact(uvec reactVec, uvec &xt)
{
  uvec nonzero;
  uvec pos_source(1), pos_disp(1), pos_repl(1);

  pos_source = find(reactVec == 1);
  pos_disp   = find(reactVec == -1);  
  pos_repl   = find(reactVec == -2);
  
  if(pos_repl.n_elem!=0) xt(pos_repl) = xt(pos_disp);
  xt(pos_disp) = xt(pos_source);
//  return(xt);
}

void replaceAndPushGillespie_indiv(umat &nu, mat &a_mat, uvec &x0, vec &time_vec, umat& simulation_data)
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
  
  a_vec = a_mat.col(0);
  // Start the simulation
  do{   
      // Recalculate reaction rates
      //a_vec = a_coeff%xt.elem(indx_rates_a);
  
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
      updateReact(nu.col(next_reac), xt);
      //xt += nu.col(next_reac); 
  }while(t_i<time_vec.n_elem);
  
//  PutRNGstate();
//  return(simulation_data);
}

umat transformMat(umat inMat)
{
  uvec sums_aux;
  umat inMat_c, inMat_b;
  umat outmat= zeros<umat>(inMat.n_rows+2, inMat.n_cols);
  inMat_c = inMat.rows(0,inMat.n_rows/2-1);
  inMat_b = inMat.rows(inMat.n_rows/2,inMat.n_rows-1);  
  
  sums_aux = sum(inMat_c.t(),1);
  // Do centre SC
  for(int i=0;i<inMat.n_cols;i++) outmat(sums_aux(i),i) = 1;
  
  sums_aux = sum(inMat_b.t(),1) + inMat.n_rows/2+1;
//  cout << sum(inMat_b.t(),1);
  // Do border SC
  for(int i=0;i<inMat.n_cols;i++) outmat(sums_aux(i),i) = 1;
  
  return(outmat);
}


// [[Rcpp::export]]
umat replaceAndPushGillespie(int numSim, umat nu, mat a_mat, uvec x0, vec time_vec)
{
  int sim_i;
  umat simulation_data_all = zeros<umat>(x0.n_elem+2, time_vec.n_elem);
  umat simulation_data_aux = zeros<umat>(x0.n_elem+2, time_vec.n_elem);
  umat simulation_data     = zeros<umat>(x0.n_elem, time_vec.n_elem);
  // .. Get random Seed  
  GetRNGstate();
  for(sim_i = 0; sim_i<numSim; sim_i++)
  {
    replaceAndPushGillespie_indiv(nu, a_mat, x0, time_vec, simulation_data); // output simulation_data with spatial data
    // Transform output to number of cells (border and center) labelled at each time point
    simulation_data_aux = transformMat(simulation_data);
    simulation_data_all = simulation_data_all + simulation_data_aux;
  }
  // .. Return random seed    
  PutRNGstate();
  return(simulation_data_all);
}


