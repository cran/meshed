#include "meshed.h"

using namespace std;

void Meshed::deal_with_tausq(MeshDataLMC& data, bool ref_pardata){
  // ref_pardata: set to true if this is called without calling deal_with_Lambda first
  if(forced_grid & arma::all(familyid == 0)){
    gibbs_sample_tausq_fgrid(data, ref_pardata);
  } else {
    gibbs_sample_tausq_std(ref_pardata);
  }
  
}


void Meshed::gibbs_sample_tausq_std(bool ref_pardata){
  if(verbose & debug){
    Rcpp::Rcout << "[gibbs_sample_tausq_std] start\n";
  }
  
  start = std::chrono::steady_clock::now();
  // note that at the available locations w already includes Lambda 
  
  double aprior = tausq_ab(0);
  double bprior = tausq_ab(1);
  
  arma::mat LHW = wU * Lambda.t();
  
  logpost = 0;
  for(int j=0; j<q; j++){
    if(familyid(j) == 0){
      // gibbs update
      arma::mat yrr = 
        y.submat(ix_by_q_a(j), oneuv*j) - 
        XB.submat(ix_by_q_a(j), oneuv*j) - 
        LHW.submat(ix_by_q_a(j), oneuv*j); //***
      
      double bcore = arma::conv_to<double>::from( yrr.t() * yrr );
      
      double aparam = aprior + ix_by_q_a(j).n_elem/2.0;
      double bparam = 1.0/( bprior + .5 * bcore );
      
      Rcpp::RNGScope scope;
      tausq_inv(j) = R::rgamma(aparam, bparam);
      logpost += 0.5 * (ix_by_q_a(j).n_elem + .0) * log(tausq_inv(j)) - 0.5*tausq_inv(j)*bcore;
      
      if(verbose & debug){
        Rcpp::Rcout << "[gibbs_sample_tausq] " << j << " | "
                    << aparam << " : " << bparam << " " << bcore << " --> " << 1.0/tausq_inv(j)
                    << "\n";
      }
    } else if(familyid(j) == 3){
      
      
      betareg_tausq_adapt.at(j).count_proposal();
      Rcpp::RNGScope scope;
      
      arma::vec one = arma::ones(1);
      arma::vec U_update = one * R::rnorm(0, 1);
      
      arma::vec new_tsqiv = 
        par_huvtransf_back(par_huvtransf_fwd(one*tausq_inv(j), tausq_unif_bounds.rows(oneuv * j)) + 
        betareg_tausq_adapt.at(j).paramsd * U_update, tausq_unif_bounds.rows(oneuv * j));
      
      double new_tsqi = new_tsqiv(0);
      //Rcpp::Rcout << arma::size(offsets) << " " << arma::size(XB) << " " << arma::size(LHW) << " " << arma::size(y) << endl;
      
      arma::vec start_logpost_vec = arma::zeros(ix_by_q_a(j).n_elem);
      arma::vec new_logpost_vec = arma::zeros(ix_by_q_a(j).n_elem);
      
#ifdef _OPENMP
#pragma omp parallel for 
#endif
      for(unsigned int ix=0; ix<ix_by_q_a(j).n_elem; ix++){
        int i = ix_by_q_a(j)(ix);
        
        double sigmoid = 1.0/(1.0 + exp(-offsets(i, j) - XB(i, j) - LHW(i, j)));
        
        start_logpost_vec(ix) = betareg_logdens(y(i, j), sigmoid, tausq_inv(j));
        new_logpost_vec(ix) = betareg_logdens(y(i, j), sigmoid, new_tsqi);
      }
      
      double start_logpost = arma::accu(start_logpost_vec);
      double new_logpost = arma::accu(new_logpost_vec);
      
      double prior_logratio = 0;
      
      if(aprior != 0){
        // for(int i=0; i<q; i++){
        //   prior_logratio += aprior * 
        //     (- log(new_tausq(i)) - log(tausq_inv(i)));
        // }
        
        prior_logratio = calc_prior_logratio(one * new_tsqi, one * tausq_inv(j), aprior, bprior);
      }
      
      if(std::isnan(prior_logratio)){
        Rcpp::Rcout << "NaN value from prior on tausq: a=" << aprior << " b=" << bprior << endl;
        Rcpp::stop("Terminated.");
      }
      
      double jacobian  = calc_jacobian(one * new_tsqi, one * tausq_inv(j), tausq_unif_bounds.rows(oneuv * j));
      double logaccept = new_logpost - start_logpost + 
        prior_logratio +
        jacobian;
      // 
      // Rcpp::Rcout << "new: " << new_logpost << " old " << start_logpost << endl;
      // 
      // // 
      // Rcpp::Rcout << "new " << new_tsqi << " vs " << tausq_inv(j) << endl
      //             << tausq_unif_bounds << endl
      //             << prior_logratio << endl
      //             << jacobian << endl;
      // 
      bool accepted = do_I_accept(logaccept);
      if(accepted){
        betareg_tausq_adapt.at(j).count_accepted();
        // make the move
        tausq_inv(j) = new_tsqi;
      } 
      
      betareg_tausq_adapt.at(j).update_ratios();
      betareg_tausq_adapt.at(j).adapt(U_update, exp(logaccept), brtausq_mcmc_counter(j)); 
      brtausq_mcmc_counter(j) ++;
    }
  }
  
  if(arma::any(familyid == 3)){
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(int i = 0; i<n_ref_blocks; i++){
      int r = reference_blocks(i);
      int u = block_names(r)-1;
      update_lly(u, alter_data, LambdaHw);
      //if(ref_pardata){
        update_lly(u, param_data, LambdaHw);
      //}
    }
    alter_data.loglik_w = arma::accu(alter_data.logdetCi_comps) + 
      arma::accu(alter_data.loglik_w_comps) + arma::accu(alter_data.ll_y); //***
    param_data.loglik_w = arma::accu(param_data.logdetCi_comps) + 
      arma::accu(param_data.loglik_w_comps) + arma::accu(param_data.ll_y); //***
  }
  
  
  if(verbose & debug){
    end = std::chrono::steady_clock::now();
    Rcpp::Rcout << "[gibbs_sample_tausq] "
                << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() 
                << "us.\n";
  }
  
}

void Meshed::gibbs_sample_tausq_fgrid(MeshDataLMC& data, bool ref_pardata){
  if(verbose & debug){
    Rcpp::Rcout << "[gibbs_sample_tausq_fgrid] start (sampling via Robust adaptive Metropolis)\n";
  }
  
  start = std::chrono::steady_clock::now();
  
  double aprior = tausq_ab(0);
  double bprior = tausq_ab(1);
  
  //arma::vec tau = sqrt(1.0/tausq_inv);
  
  tausq_adapt.count_proposal();
  Rcpp::RNGScope scope;
  
  arma::vec U_update = mrstdnorm(q, 1);
  
  arma::vec new_tausq = 
    par_huvtransf_back(par_huvtransf_fwd(1.0/tausq_inv, tausq_unif_bounds) + 
    tausq_adapt.paramsd * U_update, tausq_unif_bounds);
  
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(int i = 0; i<n_ref_blocks; i++){
    int r = reference_blocks(i);
    int u = block_names(r)-1;
    //for(int i=0; i<n_blocks; i++){
    //int u = block_names(i)-1;
    calc_DplusSi(u, alter_data, Lambda, 1.0/new_tausq);
    update_lly(u, alter_data, LambdaHw);
    
    //calc_DplusSi(u, param_data, tausq_inv);
    if(ref_pardata){
      update_lly(u, param_data, LambdaHw);
    }
  }
  
  double new_logpost = arma::accu(alter_data.ll_y);
  double start_logpost = arma::accu(param_data.ll_y);
  
  double prior_logratio = 0;
  
  if(aprior != 0){
    // for(int i=0; i<q; i++){
    //   prior_logratio += aprior * 
    //     (- log(new_tausq(i)) - log(tausq_inv(i)));
    // }
    
    prior_logratio = calc_prior_logratio(new_tausq, 1.0/tausq_inv, aprior, bprior);
  }
  
  double jacobian  = calc_jacobian(new_tausq, 1.0/tausq_inv, tausq_unif_bounds);
  double logaccept = new_logpost - start_logpost + 
    prior_logratio +
    jacobian;
  // 
  // Rcpp::Rcout << "new " << new_tausq(0) << " vs " << 1.0/tausq_inv << endl
  //             << tausq_unif_bounds << endl
  //             << prior_logratio << endl
  //             << jacobian << endl;
  // 
  bool accepted = do_I_accept(logaccept);
  if(accepted){
    tausq_adapt.count_accepted();
    // make the move
    tausq_inv = 1.0/new_tausq;
    param_data.DplusSi = alter_data.DplusSi;
    param_data.DplusSi_c = alter_data.DplusSi_c;
    param_data.DplusSi_ldet = alter_data.DplusSi_ldet;
    
    logpost = new_logpost;
  } else {
    logpost = start_logpost;
  }
  
  tausq_adapt.update_ratios();
  tausq_adapt.adapt(U_update, exp(logaccept), tausq_mcmc_counter); 
  tausq_mcmc_counter ++;
  
  if(verbose & debug){
    end = std::chrono::steady_clock::now();
    Rcpp::Rcout << "[gibbs_sample_tausq_fgrid] " << 
      tausq_adapt.accept_ratio << " average acceptance rate, "
                               << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() 
                               << "us.\n";
  }
}
