#' SpatialModelEstimator
#' 
#' @details 
#' Return Estimates and objects from a model based GMRF approach
#'
#' @param spatial_df data.frame 
#' @param formula class(formula)  e.g. formula = y ~ 1  
#' @param mesh data.frame 
#' @param extrapolation_grid data.frame we expect columns names "x", "y", "area"
#' @param family 0 = Poisson, 1 = Negative Binomial, 2 = Gaussian, 3 = Gamma 
#' @param link link function 0 = log, 1 = logit, 2 = probit, 3 = inverse, 4 = identity
#' @param bias_correct <bool> whether bias correction is needed in deriving expectations and standard errors for 
#' @param control list of settings for nlminb
#' @param convergence list of convergence setting to check
#' @param start_par_list list of inital values for optimisaton, if you have conversion issues can be useful to look at alternative starting values
#' @export
#' @import sp
#' @import INLA
#' @import RANN
#' @import TMB
#' @return: list of estimated objects.


SpatialModelEstimator = function(spatial_df, formula, mesh, extrapolation_grid, family, link, bias_correct = F, control = list(eval.max = 10000, iter.max = 10000), convergence = list(grad_tol = 0.001), start_par_list = NULL) {
  if(class(spatial_df) != "SpatialPointsDataFrame")
    stop(paste0("spatial_df needs to be of class SpatialPointsDataFrame, see the example for more information"))
  if(!family %in% c(0:3))
    stop(paste0("family needs to be a value from 0 to 3, for a valid distribution"))
  if(!link %in% c(0:4))
    stop(paste0("link needs to be a value from 0 to 4, for a valid link function"))
  if(class(formula) != "formula") 
    stop(paste0("formula not of class formula, often due to wrapping as string, i.e. having quotation marks."))
  # get response variable
  repsons_var_label = all.vars(formula)[attr(terms(formula), "response")]
  if(!all(c("x","y","area") %in% colnames(extrapolation_grid))) 
    stop(paste0("extrapolation_grid: needs colnames 'x', 'y', 'area'"))
  if(!"area" %in% colnames(spatial_df@data))
    stop(paste0("spatial_df@data: needs colnames 'area'"))
  if(!repsons_var_label %in% colnames(extrapolation_grid))
    stop(paste0("extrapolation_grid: needs colname ", repsons_var_label, ", this is a dummy variable just for creating extrapolation matrix."))

  ## map mesh to observations
  A = inla.spde.make.A(mesh, loc = cbind(coordinates(spatial_df[, 1]), coordinates(spatial_df[, 2])))
  ## map mesh to extrapolation grid
  Proj <- inla.mesh.projector(mesh, loc = cbind(extrapolation_grid$x, extrapolation_grid$y))
  P = Proj$proj$A
  ## Create Sparse Matern objects to pass to TMB
  spde = inla.spde2.matern(mesh, alpha = 2)
  ## model formula
  model_matrix = model.matrix(object = formula, data = spatial_df@data)
  ## will need to create a dummy variable for the response variable in the extrapolation grid.
  proj_model_matrix = model.matrix(object = formula, data = extrapolation_grid)
  ## Index prediction sampled locations, these wont be included with 
  proj_indicator = rep(1,nrow(extrapolation_grid))
  proj_ndx = RANN::nn2(query = coordinates(spatial_df), data = cbind(extrapolation_grid$x, extrapolation_grid$y), k = 1)
  proj_indicator[proj_ndx$nn.idx] = 0
  
  ## create TMB data object
  data_tmb <- list(y = spatial_df@data[, repsons_var_label],
               area = spatial_df@data$area,
               spde = spde$param.inla[c("M0","M1","M2")],
               Proj = P,
               A = A,
               model_matrix = model_matrix,
               family = family,
               link = link,
               proj_area = extrapolation_grid$area,
               proj_model_matrix = proj_model_matrix,
               pred_indicator = proj_indicator
  )
  ## create TMB param object
  params_tmb = list()
  if(is.null(start_par_list)) {
    params_tmb$omega = rep(0, spde$n.spde)
    params_tmb$betas = c(mean(data_tmb$y / data_tmb$area), rep(0, ncol(model_matrix) - 1))
    ## get range of domain, based on mesh
    range_domain = diff(range(mesh$loc))
    params_tmb$ln_kappa = log(sqrt(8) / (range_domain * 0.25)) ## start by distance at which spatial autocorrelation (rho) = 0.1. based on 25% of domain
    ## variance terms half for phi, half for GMRF.
    params_tmb$ln_tau = log(1 / (0.5 * var(data_tmb$y / data_tmb$area)))
    params_tmb$ln_phi = log(0.5 * sd(data_tmb$y / data_tmb$area))
    if(family == 1) {
      params_tmb$ln_tau = 0 ## phi -> 0 = Poission
    } else if (family == 2) {
      params_tmb$ln_phi = log(0.5 * sd(data_tmb$y / data_tmb$area))
    } else if (family == 3) {
      params_tmb$ln_phi = log(0.3) # coeffecient of variation
    }
  } else {
    if(class(start_par_list) != "list")
      stop(paste0("If you supply start values it needs to be of class 'list'"))
    if(!all(names(start_par_list) %in% c("omega", "betas", "ln_kappa", "ln_tau", "ln_phi")))
      stop(paste0("you supplied start_par_list, needs to have all these parameters 'omega', 'betas', 'ln_kappa', 'ln_tau', 'ln_phi'"))
    params_tmb = start_par_list
  }
  ## optimise - model with no GMRF
  fixed_pars_no_gmrf = list(ln_kappa = as.factor(NA), ln_tau = as.factor(NA), omega = as.factor(rep(NA, spde$n.spde)))
  if(!family %in% c(1, 2, 3)) 
    fixed_pars_no_gmrf$ln_phi = as.factor(NA)
  
  ## Estiamte initial model without GMRF 
  data_tmb$model = "GMRFModel"
  
  obj_init = tryCatch (
    expr <- TMB::MakeADFun(data = data_tmb, parameters = params_tmb, map = fixed_pars_no_gmrf, DLL = "spatialsurvey_TMBExports", silent = T),
    error = function(e) {e}
  )
  if(inherits(obj_init, "error")) 
    stop(paste0("Failed to create obj_init, issue with MakeADFun() call error = ", obj_init$message))
  
  opt_init = tryCatch (
    expr = nlminb(obj_init$par, obj_init$fn, obj_init$gr, control = control),
    error = function(e) {e}
  )  
  if(inherits(opt_init, "error")) 
    stop(paste0("Failed to optimise opt_init, issue with nlminb() call. Error = ", opt_init$message))
  ## check convergence
  max_grad_this = max(abs(obj_init$gr(opt_init$par)));
  grad_tol = 0.001
  if(!is.null(convergence)) 
    grad_tol = convergence$grad_tol
  ## check hessian invertable and gradient above some satisfied 
  init_convergence = "There is no evidence that the model is not converged"
  if (is.character(try(chol(optimHess(opt_init$par, fn = obj_init$fn, gr = obj_init$gr)), silent = TRUE)) | (max_grad_this > grad_tol))
    init_convergence = "Failed conergence for non GMRF model";
  
  ## pull out estimates
  sd_report_init = TMB::sdreport(obj = obj_init)
  est_init = sd_report_init$value["fitted_non_sample_Total"]
  sd_init = sd_report_init$sd[names(sd_report_init$value) == "fitted_non_sample_Total"]
  est_ln_init = sd_report_init$value["log_fitted_non_sample_Total"]
  sd_ln_init = sd_report_init$sd[names(sd_report_init$value) == "log_fitted_non_sample_Total"]
  MLE_init_report = obj_init$report(obj_init$env$last.par.best)
  
  params_tmb$betas = as.numeric(opt_init$par[names(opt_init$par) %in% "betas"])
  
  if (any(names(opt_init$par) %in% "ln_phi"))
    params_tmb[["ln_phi"]] = opt_init$par[names(opt_init$par) %in% "ln_phi"]
  
  ## Now do the GMRF model and estimate
  obj_gmrf = NULL
  if(family %in% c(1, 2, 3)) {
    obj_gmrf = tryCatch(
      expr <- TMB::MakeADFun(data = data_tmb, random = c("omega"), parameters = params_tmb , DLL = "spatialsurvey_TMBExports", silent = T),
      error = function(e) {e}
    )
  } else {
    obj_gmrf = tryCatch(
      expr <- TMB::MakeADFun(data = data_tmb, random = c("omega"), parameters = params_tmb, map = list(ln_phi = as.factor(NA)), DLL = "spatialsurvey_TMBExports", silent = T),
      error = function(e) {e}
    )
  }
  if(inherits(obj_gmrf, "error")) 
    stop(paste0("Failed to create obj_gmrf, issue with MakeADFun() call error = ", obj_gmrf$message))
  
  ## Optimise
  opt_gmrf = tryCatch(
    expr = nlminb(obj_gmrf$par, obj_gmrf$fn, obj_gmrf$gr, control = control),
    error = function(e) {e}
  )
  if(inherits(opt_gmrf, "error")) 
    stop(paste0("Failed to optimise opt_gmrf, issue with nlminb() call. Error = ", opt_gmrf$message))
  ## Check convergence
  temp_pars = opt_gmrf$par
  max_grad_this = max(abs(obj_gmrf$gr(temp_pars)));
  
  ## start building the return object
  return_list = list(init_convergence = init_convergence, MLE_init_report = MLE_init_report)
  return_list$start_vals = params_tmb;
  return_list$sd_report_init = sd_report_init; ## no random effects so don't need to bias correct.
  return_list$init_opt = opt_init
  if (is.character(try(chol(optimHess(temp_pars, fn = obj_gmrf$fn, gr = obj_gmrf$gr)), silent = TRUE)) | (max_grad_this > grad_tol)) {
    ## try a different kappa
    temp_pars["ln_kappa"] = log(sqrt(8) / (range_domain * 0.1)) ## 10% of domain range
    temp_opt = tryCatch(
      expr = nlminb(temp_pars, obj_gmrf$fn, obj_gmrf$gr, control = control),
      error = function(e) {e}
    )
    temp_pars <<- temp_opt$par
    opt_gmrf$iterations = opt_gmrf$iterations + temp_opt$iterations
    if (is.character(try(chol(optimHess(temp_pars, fn = obj_gmrf$fn, gr = obj_gmrf$gr)), silent = TRUE))) {
      temp_pars["ln_kappa"] = log(sqrt(8) / (range_domain * 0.5)) ## 50% of domain range
      temp_opt = tryCatch(
        expr = nlminb(temp_pars, obj_gmrf$fn, obj_gmrf$gr, control = control),
        error = function(e) {e}
      )
      temp_pars <<- temp_opt$par
      opt_gmrf$iterations = opt_gmrf$iterations + temp_opt$iterations
    }
  }
  ## final convergence parameters
  max_grad_this = max(abs(obj_gmrf$gr(temp_pars)));
  hess = optimHess(temp_pars, fn = obj_gmrf$fn, gr = obj_gmrf$gr)
  gmrf_convergence = ""
  if (is.character(try(chol(hess), silent = TRUE))) {
    print("Hessian is not positive definite, so standard errors are not available for GMRF use init model");
    gmrf_convergence = "Hessian is not positive definite, so standard errors are not available. Look at eigen values from hession, to see problem params."
    return_list$hess_eigen_vals = eigen(hess)
    return_list$gmrf_failed_opt = opt_gmrf
  } else {
    if (max_grad_this > grad_tol) {
      gmrf_convergence = paste0("Gradient of MLE params > grad_tol = '", grad_tol, "standard errors are probably useable, but this does indicate non convergence")
    } else {
      gmrf_convergence = "There is no evidence that the model is not converged"
    }
    gmrf_sd_report = NULL
    if (bias_correct) {
      gmrf_sd_report = TMB::sdreport(obj_gmrf, bias.correct= TRUE, bias.correct.control = list(vars_to_correct = c("fitted_non_sample_Total", "log_fitted_non_sample_Total"), sd = TRUE),  getJointPrecision = F)
    } else {
      gmrf_sd_report = TMB::sdreport(obj_gmrf,  getJointPrecision = F)
    }
    return_list$gmrf_sd_report = gmrf_sd_report
    return_list$gmrf_max_grad = max_grad_this
    return_list$MLE_gmrf = obj_gmrf$report(obj_gmrf$env$last.par.best)
    return_list$MLE_par = obj_gmrf$env$last.par.best ## opt_gmrf won't have random effects, this will have bayes estiamtes for random effects.
    return_list$gmrf_opt = opt_gmrf
  }
  return_list$gmrf_convergence = gmrf_convergence
  return_list$tmb_data = data_tmb ## with this and MLE pars can generate thier on obj and play around.
  return(return_list)
}


