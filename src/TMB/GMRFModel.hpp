/// @file GMRFModel.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

//#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

/*
* Valid link functions
*/
enum valid_link {
  log_link                 = 0,
  logit_link               = 1,
  probit_link              = 2,
  inverse_link             = 3,
  identity_link            = 4
};
/*
* Valid distributions
*/
enum valid_family {
  poisson = 0,
  negative_binomial = 1,
  gaussian    = 2,
  gamma  = 3
};

/*
* Apply inverse link function for general cases
*/
template<class Type>
Type inverse_linkfun(Type& eta, int& link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = exp(eta);
    break;
  case identity_link:
    ans = eta;
    break;
  case logit_link:
    ans = invlogit(eta);
    break;
  case probit_link:
    ans = pnorm(eta);
    break;
  case inverse_link:
    ans = Type(1) / eta;
    break;
    // TODO: Implement remaining links
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}

// name of function below **MUST** match filename
// (here it's ModelA)
template <class Type>
Type GMRFModel(objective_function<Type>* obj) {
  using namespace R_inla;
  using namespace Eigen;
  using namespace density;
  // Observation level inputs
  DATA_VECTOR( y );                       // Response yariable length = n
  DATA_VECTOR_INDICATOR (keep , y);       // NOTICE " keep " this is for OSA residuals
  DATA_MATRIX( model_matrix ) ;           // model matrix To account for covariates dim[n,p]
  DATA_VECTOR( area );                    // Area for each observation sampling unit
  DATA_INTEGER( family );                 // response distribution 0 = Poisson, 1 = Negative Binomial, 2 = Gaussian, 3 = Gamma
  DATA_INTEGER( link ) ;                  // 0 = log, 1 = logit, 2 = probit, 3 = inverse, 4 = identity
  DATA_INTEGER( simulate_state ) ;        // not used in estimation, but if you build a model, you can use the simulate() call on it.
  // extrapolation inputs
  DATA_VECTOR( proj_area ) ;              // Area for each projection grid
  DATA_MATRIX( proj_model_matrix ) ;      // projection model matrix  dim[nrow(P),p]
  DATA_VECTOR( pred_indicator );          // for each row of Proj, whether it is observed or not
  // GMRF inputs
  DATA_SPARSE_MATRIX( Proj ) ;            // Sparse projection matrix maps nodes to each projection grid
  DATA_SPARSE_MATRIX( A ) ;               // Sparse matrix that maps, notes to each observation location
  DATA_STRUCT( spde, spde_t );            // spde object from INLA

  // Parameters
  PARAMETER_VECTOR( betas );              // Shoule be expanded for covariates, then need model matrix and projection model matrix
  PARAMETER( ln_kappa );                  // spatial decay/range parameter
  PARAMETER( ln_tau );                    // precision of GMRF covariance, note marginal variance of the GMRF is a function of kappa and tau 
  PARAMETER_VECTOR( omega );              // spatial random effect/latent variables
  PARAMETER( ln_phi );                    // overdispersion parameter, cv = gamma, sd = normal, phi = Neg Bin

  int i;
  vector<Type> nll(2);
  nll.setZero();
  // Transform parameters
  Type tau = exp(ln_tau);
  Type kappa = exp(ln_kappa);
  Type phi = exp(ln_phi);
  // Random effect contributions
  Eigen::SparseMatrix<Type> Q = Q_spde(spde, kappa);
  nll(0) = GMRF(Q)(omega);                              // Then the density is proportional to |Q|^.5*exp(-.5*x'*Q*x) mean = 0
  // The linear component of the model
  vector<Type> omega_i = (A * omega) / tau;             // maps mesh to observations usign A this is where the sccaling comes in
  vector<Type> eta = model_matrix * betas + omega_i;    // link X_beta
  // Apply link function
  vector<Type> mu(eta.size());
  for (i = 0; i < mu.size(); i++)
    mu(i) = area(i) * inverse_linkfun(eta(i), link);

  /////////////////
  // Extrapolation
  /////////////////
  vector<Type> omega_proj = (Proj * omega) / tau;                         // Project GF to all projected points account for tau
  vector<Type> linear_proj = proj_model_matrix * betas + omega_proj;     // this is in log space
  vector<Type> mu_proj(linear_proj.size());// = Proj_Area * exp(linear_proj);
  for (i = 0; i < linear_proj.size(); i++)
    mu_proj(i) = proj_area(i) * inverse_linkfun(linear_proj(i), link);

  Type fitted_non_sample_Total = (mu_proj * pred_indicator).sum();
  Type log_fitted_non_sample_Total = log(fitted_non_sample_Total);
  Type log_fitted_Total = log(fitted_non_sample_Total + y.sum());
  
  /////////////////
  // Observation likelihood
  /////////////////
  Type s1, s2;
  Type tmp_loglik;
  for (i = 0; i < y.size(); i++){
    if ( !isNA(y(i)) ) {
      switch (family) {
      case gaussian:
        tmp_loglik = dnorm(y(i), mu(i), phi, true);
        //SIMULATE{y(i) = rnorm(mu(i), sqrt(phi));}
        break;
      case poisson:
        tmp_loglik = dpois(y(i), mu(i), true);
        //SIMULATE{y(i) = rpois(mu(i));}
        break;
      case negative_binomial:
        s1 = log(mu(i));                       // log(mu)
        s2 = 2. * s1 - ln_phi;                         // log(var - mu)
        tmp_loglik = dnbinom_robust(y(i), s1 , s2, true);
        //SIMULATE{y(i) = rbinom(Type(1), mu(i));}
        break;
      case gamma:
        s1 = phi;           // shape
        s2 = mu(i) / phi;   // scale
        tmp_loglik = dgamma(y(i), s1, s2, true);
        //SIMULATE{y(i) = rgamma(s1, s2);}
        break;
      default:
        error("Family not implemented!");
      } // End switch
      // Add up
      nll(1) -= keep(i) * tmp_loglik;
    }
  }
  
  SIMULATE {
    // Simulate the population over the entire domain.
    vector<Type> y_sim(mu_proj.size());
    if(simulate_state == 1) {
      // Simulate a completley new GF
      GMRF(Q).simulate(omega);
    }
    omega_proj = (Proj * omega) / tau; // Project GF to all projected points account for tau
    linear_proj = proj_model_matrix * betas + omega_proj;     // this is in log space
    for (i = 0; i < linear_proj.size(); i++)
      mu_proj(i) = proj_area(i) * inverse_linkfun(linear_proj(i), link);
    
    for (i = 0; i < y_sim.size(); i++){
      if ( !isNA(y_sim(i)) ) {
        switch (family) {
        case gaussian:
          y_sim(i) = rnorm(mu_proj(i), phi);
          break;
        case poisson:
          y_sim(i) = rpois(mu_proj(i));
          break;
        case negative_binomial:
          s1 = mu_proj(i) * (Type(1.0) + mu_proj(i) / phi); // over_dispersion + lambda_sim;
          y_sim(i) = rnbinom2(mu_proj(i), s1);
          break;
        case gamma:
          s1 = phi;           // shape
          s2 = mu_proj(i) / phi;   // scale
          y_sim(i) = rgamma(s1, s2);
          break;
        default:
          error("Family not implemented!");
        } // End switch
      }
    }
    /////////////////
    // Reports
    /////////////////
    REPORT(y_sim);
  }
  
  /////////////////
  // Reports
  /////////////////
  REPORT( fitted_non_sample_Total );
  ADREPORT( fitted_non_sample_Total );   // We can use THorson epsilon method for ADREPORT
  REPORT( log_fitted_non_sample_Total );
  ADREPORT( log_fitted_non_sample_Total );
  REPORT( log_fitted_Total );
  ADREPORT( log_fitted_Total );
  //ADREPORT( mu_proj );
  REPORT( linear_proj );
  REPORT( omega_proj );
  REPORT( mu_proj );
  REPORT( omega ); // should be overwritten with simulated vals.
  REPORT( nll );
  REPORT( mu );
  REPORT( eta );
  double nu = 1.0;            // nu = alpha-d/2 = 2-1 by eqn (2) in Lindgren
  Type rho = sqrt(8*nu) / kappa;  // Distance at which correlation has dropped to 0.1 (p.4 in Lindgren)
  Type  marginal_variance = pow(exp(.5*log(Type(1.0)/(Type(4.0)*M_PI)) - ln_kappa - ln_tau),2);
  REPORT( marginal_variance );
  REPORT( rho );
  ADREPORT( tau );
  ADREPORT( marginal_variance );
  ADREPORT( kappa );
  REPORT( phi );
  ADREPORT( phi );
  REPORT( tau );
  REPORT( kappa );
  REPORT( marginal_variance );
  return nll.sum();
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
