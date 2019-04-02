#' @docType package
#' @importFrom stats cov.wt pnorm qnorm quantile p.adjust as.formula
#' @importFrom survey svyvar svymean svyvar
#' @import Matrix
#' @import lavaan
#' @import parallel
#' @importFrom lavaan lavNames sem lav_matrix_bdiag summary coef
#' @importFrom parallel makePSOCKcluster clusterSetRNGStream stopCluster parLapply
#' @importFrom parallel mclapply

#' @name chisq.BRR
#' @title Adjust the model fit statistics
#' @description
#' This function is used to adjust model fit statistics for complex surveys with balanced repeated replications \cite{(Oberski, 2014; Satorra & Muthen, 1995)}. It saves time to only obtain the model fit statistics during the model selection stage.
#' @param model The model being fitted. It is written in lavaan model syntax \cite{(Rosseel, 2012)}.
#' @param lavaan.fit The model fit results using 'ML' estimator with sample main weights, but without adjusting the fit statistics or standard errors for complex surveys. Note that it is a lavaan object.
#' @param data The raw data including the variables of interest and the survey weights. It should be a dataset or dataframe.
#' @param mwgtname The variable name indicating the sample main weight in the dataset. See balanced repeated replications method \cite{(Wolter, 2007)} for more information about the main weight.
#' @param repwgtnames The variable names indicating the set of replicate weights in the dataset. See balanced repeated replications method \cite{(Wolter, 2007)} for more information about the replicate weight.
#' @param fayfactor The fayfactor used in the standard error calculation by fay's method \cite{(Fay & Train, 1995; Judkins, 1990)} for balanced repeated replications. Fayfactor is a value between 0 and 1. The default is 0.5.
#' @param estimator The method used to estimate the model. 'ML' is the default option and the only available option for current version. It is not required.
#' @param test The method used to generate adjusted standard errors. 'satorra.bentler' is the default option and the only available option for current version. It is not required.
#' @return The model fit results as a lavaan object \cite{(Rosseel, 2012)} with the adjusted model fit statistics.
#' @references Fay, R. E., & Train, G. F. (1995). Aspects of survey and model-based postcensal estimation of income and poverty characteristics for states and counties. In Proceedings of the Section on Government Statistics, American Statistical Association, Alexandria, VA (pp. 154-159).
#' @references Judkins, D. R. (1990). Fay’s method for variance estimation.Journal of Official Statistics,6(3), 223-239.
#' @references Oberski, D. (2014). lavaan. survey: An R package for complex survey analysis of structural equation models. Journal of Statistical Software, 57(1), 1-27. DOI:10.18637/jss.v057.i01
#' @references Rosseel, Y. (2012). Lavaan: An R package for structural equation modeling and more. Version 0.5–12 (BETA). Journal of statistical software, 48(2), 1-36. DOI:10.18637/jss.v048.i02
#' @references Satorra, A., & Muthen, B. (1995). Complex sample data in structural equation modeling. Sociological methodology, 25(1), 267-316.
#' @references Wolter, K. (2007). Introduction to variance estimation. New York, NY: Springer.
#' @examples
#' \dontshow{
#'  #Toy example for check:
#' R <- 20
#' wgtnames <- paste("repwgt", seq(0,R,by=1), sep="")
#' mwgtname=wgtnames[1]
#' repwgtnames=wgtnames[2:(R+1)]
#' fayfactor=0.5
#' model1 <- ' # outcome
#'              numcg ~ u2*1 + c*workban + b*sp_adltban
#'            # mediator
#'              sp_adltban ~ u1*1 + a*workban
#'            # indirect effect (a*b)
#'              ab := a*b
#'            # total effect
#'              total := c + (a*b)
#'           '
#' fit <- lavaan::sem(model=model1, data=MedData, estimator='ML', test='standard')
#' chisq.BRR(model1,fit,MedData,mwgtname, repwgtnames)
#' }
#' \donttest{
#' R <- 160
#' wgtnames <- paste("repwgt", seq(0,R,by=1), sep="")
#' mwgtname=wgtnames[1]
#' repwgtnames=wgtnames[2:(R+1)]
#' fayfactor=0.5
#'
#' model3 <- ' # outcome
#'             numcg ~ u0*1 + c*workban + b1*sp_adltban + b2*sp_kidsban
#'            # mediator
#'               sp_adltban ~ u1*1 + a1*workban
#'               sp_kidsban ~ u2*1 + a2*workban
#'            # indirect effect (a*b)
#'               a1b1 := a1*b1
#'              a2b2 := a2*b2
#'            # total effect
#'              total := c + (a1*b1) + (a2*b2)
#'           '
#'
#' fit <- lavaan::sem(model=model3, data=MedData, estimator='ML', test='standard')
#' chisq.BRR(model3,fit,MedData,mwgtname, repwgtnames)
#' #
#' # MedSurvey 1.1.0 Adjusted Model Fit Statistics using BRR
#' #
#' # chisq   df   pvalue    CFI      RMSEA      SRMR         AIC       BIC
#' #
#' # 305.25   1  0.00000   0.40561  0.27852   0.07416   88699.43   88768.45
#'
#' }
#'
##' @export
chisq.BRR <- function(model, lavaan.fit, data, mwgtname, repwgtnames, fayfactor=0.5, estimator=c('ML'), test=c('satorra.bentler')){
  if (missing(model)) stop("A model is needed.")
  if (missing(data)) stop("data are needed.")
  if (missing(mwgtname)) stop("varibale name of the main weight are needed.")
  if (missing(repwgtnames)) stop("varibale names of replicate weights are needed.")
  if (missing(lavaan.fit)) stop("a lavaan object of model fit results is needed.")
  if (!is.null(fayfactor)) {if(fayfactor >= 0 & fayfactor <1) {} else stop("fayfactor can only be a value between 0 and 1.") } else fayfactor <- 0.5
  if (!is.null(estimator)) {estimator<-match.arg(estimator)} else {estimator<-'ML'}

   #Correct the model fit statistic:
  des.rep <- survey::svrepdesign(ids=~1, weights=data[,mwgtname], data=as.data.frame(data),
                                 repweights=as.data.frame(data)[,repwgtnames], type="Fay", rho=fayfactor)
  wn <- as.integer(sum(data[,mwgtname]))
  n <- nrow(data)
  #This function Gamma.svy() is used to calculate the covariance matrix of sample covariances for complex surveys.
  #It was modified from an embedded function in Oberski, D. (2014)'s lavaan.survey.
  Gamma.svy <- function (lavaan.fit, survey.design){
    ov.names <- lavaan::lavaanNames(lavaan.fit, type = "ov", group = 1)
    ov.formula <- as.formula(paste("~", paste(ov.names, collapse = "+")))
    Dplus <- lavaan::lav_matrix_duplication_ginv(length(ov.names))
    get.stats.design <- function(survey.design, sample.nobs) {
      sample.cov <- as.matrix(survey::svyvar(ov.formula, design = survey.design, na.rm = TRUE))
      Gamma.cov <- attr(sample.cov, "var")
      Gamma.cov <- Dplus %*% Gamma.cov %*% t(Dplus)
      sample.mean <- survey::svymean(ov.formula, design = survey.design, na.rm = TRUE)
      Gamma.mean <- attr(sample.mean, "var")
      Gamma <- lavaan::lav_matrix_bdiag(Gamma.mean, Gamma.cov)
      Gamma<- Gamma * sample.nobs
      attr(sample.cov, "var") <- NULL
      tmp <- as.vector(sample.mean)
      names(tmp) <- names(sample.mean)
      sample.mean <- tmp
      stats <- list(Gamma = Gamma, sample.cov = sample.cov,
           sample.mean = sample.mean)
      stats
    }
    if (!any(class(survey.design) == "svyimputationList")) {
      sample.nobs <- lavaan::lavInspect(lavaan.fit, "nobs")
      stats <- get.stats.design(survey.design, sample.nobs)
    } else {stop("The survey design type 'svyimputationList' is currently not supported.")}

    stats
  }
  gamma.svy <- Gamma.svy(lavaan.fit, survey.design=des.rep)
  tcall <- lavaan::lavInspect(lavaan.fit, "call")
  tcall$data <- NULL
  tcall$sample.cov <- gamma.svy$sample.cov
  tcall$sample.mean <- gamma.svy$sample.mean
  tcall$sample.nobs <- n
  tcall$estimator <- estimator
  tcall$test <- test
  if (substr(estimator, 1, 2) == "ML") {  tcall$NACOV <- 	gamma.svy$Gamma }
  else if (substr(estimator, 1, 2) == "MLM") {tcall$NACOV <- 	gamma.svy$Gamma }
  main.res.new <- eval(as.call(tcall))
  cat('\tMedSurvey 1.1.0 Adjusted Model Fit Statistics using BRR  \n\n');
  fitmeasures <- as.vector(lavaan::fitMeasures(main.res.new, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr", "aic", "bic")))
  cat('chisq   df   pvalue    CFI      RMSEA      SRMR         AIC       BIC \n\n')
  cat(sprintf("%8.2f %3.0f  %6.5f   %6.5f  %6.5f   %6.5f %10.2f %10.2f \n",
              fitmeasures[1], fitmeasures[2], fitmeasures[3],fitmeasures[4], fitmeasures[5],fitmeasures[6]
              , fitmeasures[7], fitmeasures[8]))
  #return(main.res.new)
  invisible(main.res.new)
}

#' @name med.fit.BRR
#' @title Estimate the mediation effects and standard errors adjusting for complex surveys with BRR
#' @description
#' This function is used to estimate the mediation effects adjusted for complex surveys with balanced repeated replications (BRR) \cite{ (Mai, Ha, Soulakova, 2019)}.
#' @param model The model being fitted. It is written in lavaan model syntax \cite{(Rosseel, 2012)}.
#' @param data The raw data including the variables of interest and the survey weights. It should be a dataset or dataframe.
#' @param mwgtname The variable name indicating the sample main weight in the dataset. See balanced repeated replications method \cite{(Wolter, 2007)} for more information about the main weight.
#' @param repwgtnames The variable names indicating the set of replicate weights in the dataset. See balanced repeated replications method \cite{(Wolter, 2007)} for more information about the replicate weight.
#' @param fayfactor The fayfactor used in the standard error calculation by fay's method \cite{(Fay & Train, 1995; Judkins, 1990)} for balanced repeated replications. Fayfactor is a value between 0 and 1. The default is 0.5.
#' @param estimator The method used to estimate the model. 'ML' is the default option and the only available option for current version. It is not required.
#' @param test The method used to generate adjusted standard errors. 'satorra.bentler' is the default option and the only available option for current version. It is not required.
#' @param parallel Parallel computing (\code{"no"} or \code{"parallel"} or \code{"snow"}). It is "no" by default, which means it will not use parallel computing.
#' The option "parallel" is to use multiple cores in a computer for parallel computing. It is used with the number of cores (\cite{ncore}).
#' The option "snow" is to use clusters for parallel computing. It is used with the number of clusters (\cite{cl}).
#' @param ncore Number of processors used for parallel computing. By default, ncore = Sys.getenv ('NUMBER_OF_PROCESSORS').
#' @param cl Number of clusters. It is NULL by default. When it is NULL, the program will detect the number of clusters automatically.
#' @param ... Extra arguments. For example, ordered=c('z1','z2') is an argument to tell 'z1' and 'z2' are ordinal variables. It is not required.
#' @return The model fit results as a lavaan object with the adjusted estimates, standard errors, and model fit statistics. It is a lavaan object \cite{(Rosseel, 2012)}.
#' @references Fay, R. E., & Train, G. F. (1995). Aspects of survey and model-based postcensal estimation of income and poverty characteristics for states and counties. In Proceedings of the Section on Government Statistics, American Statistical Association, Alexandria, VA (pp. 154-159).
#' @references Judkins, D. R. (1990). Fay’s method for variance estimation.Journal of Official Statistics,6(3), 223-239.
#' @references Mai, Y., Ha, T., & Soulakova, J. N. (2019). Multimediation Method With Balanced Repeated Replications For Analysis Of Complex Surveys. Structural Equation Modeling: A Multidisciplinary Journal. DOI:10.1080/10705511.2018.1559065
#' @references Rosseel, Y. (2012). Lavaan: An R package for structural equation modeling and more. Version 0.5–12 (BETA). Journal of statistical software, 48(2), 1-36. DOI:10.18637/jss.v048.i02
#' @references Wolter, K. (2007). Introduction to variance estimation. New York, NY: Springer.
#' @examples
#' \dontshow{
#'  #Toy example for check:
#' R <- 20
#' wgtnames <- paste("repwgt", seq(0,R,by=1), sep="")
#' mwgtname=wgtnames[1]
#' repwgtnames=wgtnames[2:(R+1)]
#' model1 <- ' # outcome
#'              numcg ~ u2*1 + c*workban + b*sp_adltban
#'            # mediator
#'              sp_adltban ~ u1*1 + a*workban
#'            # indirect effect (a*b)
#'              ab := a*b
#'            # total effect
#'              total := c + (a*b)
#'           '
#' fit.BRR <- med.fit.BRR(model=model1, data=MedData, mwgtname=mwgtname,
#'          repwgtnames=repwgtnames, fayfactor=0.5)
#' lavaan::summary(fit.BRR)
#' }
#' \donttest{
#' R <- 160
#' wgtnames <- paste("repwgt", seq(0,R,by=1), sep="")
#' mwgtname=wgtnames[1]
#' repwgtnames=wgtnames[2:(R+1)]
#'
#' model2 <- ' # outcome
#'               numcg ~ u0*1 + c*workban + b1*sp_adltban + b2*sp_kidsban
#'             # mediator
#'               sp_adltban ~ u1*1 + a1*workban
#'               sp_kidsban ~ u2*1 + a2*workban
#'             #covariance of residuals
#'               sp_adltban ~~ sp_kidsban
#'             # indirect effect (a*b)
#'               a1b1 := a1*b1
#'               a2b2 := a2*b2
#'             # total effect
#'               total := c + (a1*b1) + (a2*b2)
#'            '
#' fit.BRR <- med.fit.BRR(model=model2, data=MedData, mwgtname=mwgtname,
#'          repwgtnames=repwgtnames, fayfactor=0.5, parallel='parallel', ncore=2)
#' lavaan::summary(fit.BRR)
#' #
#' # lavaan 0.6-3 ended normally after 41 iterations
#' #
#' # Optimization method                           NLMINB
#' # Number of free parameters                         12
#' #
#' # Number of observations                          3922
#' #
#' # Estimator                                         ML      Robust
#' # Model Fit Test Statistic                       0.000       0.000
#' # Degrees of freedom                                 0           0
#' # Minimum Function Value               0.0000000000000
#' # Scaling correction factor                                     NA
#' # for the Satorra-Bentler correction
#' #
#' # Parameter Estimates:
#' #
#' #   Information                                 Expected
#' # Information saturated (h1) model          Structured
#' # Standard Errors                                  BRR
#' #
#' # Regressions:
#' #                    Estimate  Std.Err  z-value  P(>|z|)
#' # numcg ~
#' #    workban    (c)   -0.101    0.039   -2.572    0.010
#' #    sp_adltbn (b1)   -0.253    0.048   -5.270    0.000
#' #    sp_kidsbn (b2)   -0.361    0.051   -7.006    0.000
#' # sp_adltban ~
#' #    workban   (a1)    0.069    0.018    3.753    0.000
#' # sp_kidsban ~
#' #    workban   (a2)    0.020    0.016    1.250    0.211
#' #
#' # Covariances:
#' #                    Estimate  Std.Err  z-value  P(>|z|)
#' # .sp_adltban ~~
#' #    .sp_kidsban        2.784    0.195   14.300    0.000
#' #
#' # Intercepts:
#' #                    Estimate  Std.Err  z-value  P(>|z|)
#' #   .numcg     (u0)   18.485    0.566   32.668    0.000
#' #   .sp_adltbn (u1)    4.221    0.167   25.281    0.000
#' #   .sp_kidsbn (u2)    7.926    0.143   55.272    0.000
#' #
#' # Variances:
#' #                    Estimate  Std.Err  z-value  P(>|z|)
#' #   .numcg            54.283    1.716   31.628    0.000
#' #   .sp_adltban       11.011    0.239   46.140    0.000
#' #   .sp_kidsban        9.402    0.209   44.998    0.000
#' #
#' # Defined Parameters:
#' #                    Estimate  Std.Err  z-value  P(>|z|)
#' #    a1b1             -0.017    0.006   -2.905    0.004
#' #    a2b2             -0.007    0.006   -1.234    0.217
#' #    total            -0.125    0.040   -3.169    0.002
#' }
#'
##' @export
med.fit.BRR<-function(model=NULL, data=NULL, mwgtname=NULL, repwgtnames=NULL, fayfactor=0.5,
                      estimator=c('ML'), test=c('satorra.bentler'), parallel=c('no','parallel','snow'), ncore=Sys.getenv('NUMBER_OF_PROCESSORS'), cl=NULL, ...){
  if (missing(model)) stop("A model is needed.")
  if (missing(data)) stop("data are needed.")
  if (missing(mwgtname)) stop("varibale name of the main weight are needed.")
  if (missing(repwgtnames)) stop("varibale names of replicate weights are needed.")
  if (!is.null(fayfactor)) {if(fayfactor >= 0 & fayfactor <1){} else stop("fayfactor can only be a value between 0 and 1.") } else fayfactor <- 0.5
  if (!is.null(estimator)) {estimator<-match.arg(estimator)} else {estimator<-'ML'}
  if (!is.null(test)) {test<-match.arg(test)} else {test<-'satorra.bentler'}
  if (!is.null(parallel)) {parallel <- match.arg(parallel) } else parallel<-'no'
  if (!is.null(ncore)) {if(ncore >= 1){ncore <- as.integer(ncore)} else {ncore<-1} } else ncore<-4

   #test the model and read the ov names:
  temp.res<-lavaan::sem(model=model, data=data, estimator=estimator, test="standard", ...)
  idx <- 1:length( temp.res@ParTable$lhs )
  cnames<-paste(temp.res@ParTable$lhs[idx], temp.res@ParTable$op[idx], temp.res@ParTable$rhs[idx])
  ovnames<-lavaan::lavNames(temp.res,'ov')
  i <- 1
  for (i in 1: length(ovnames) ){
    if (!is.numeric(data[1,ovnames[i]])) stop(paste('value of ',ovnames[i], 'is not numerical', sep=''))
  }

  #calculate the weighted sample mean and cov
  wn <- as.integer(sum(data[,mwgtname]))
  n <- nrow(data)
  wsmpcov <- stats::cov.wt(x=data[,ovnames], wt = data[,mwgtname], cor = FALSE, center = TRUE, method = "unbiased")
  main.res<- lavaan::sem(model=model, sample.cov = wsmpcov$cov, sample.mean = wsmpcov$center, sample.nobs = n, test="standard", ...) # lavaan::summary(main.res)


  #Correct the model fit statistic:
  lavaan.fit=main.res
  des.rep <- survey::svrepdesign(ids=~1, weights=data[,mwgtname], data=as.data.frame(data),
                                 repweights=as.data.frame(data)[,repwgtnames], type="Fay", rho=fayfactor)
  wn <- as.integer(sum(data[,mwgtname]))
  n <- nrow(data)
  #This function is used to calculate the covariance matrix of sample covariances for complex surveys.
  Gamma.svy <- function (lavaan.fit, survey.design){
    ov.names <- lavaan::lavaanNames(lavaan.fit, type = "ov", group = 1)
    ov.formula <- as.formula(paste("~", paste(ov.names, collapse = "+")))
    Dplus <- lavaan::lav_matrix_duplication_ginv(length(ov.names))
    get.stats.design <- function(survey.design, sample.nobs) {
      sample.cov <- as.matrix(survey::svyvar(ov.formula, design = survey.design, na.rm = TRUE))
      Gamma.cov <- attr(sample.cov, "var")
      Gamma.cov <- Dplus %*% Gamma.cov %*% t(Dplus)
      sample.mean <- survey::svymean(ov.formula, design = survey.design, na.rm = TRUE)
      Gamma.mean <- attr(sample.mean, "var")
      Gamma <- lavaan::lav_matrix_bdiag(Gamma.mean, Gamma.cov)
      Gamma<- Gamma * sample.nobs
      attr(sample.cov, "var") <- NULL
      tmp <- as.vector(sample.mean)
      names(tmp) <- names(sample.mean)
      sample.mean <- tmp
      stats <- list(Gamma = Gamma, sample.cov = sample.cov,
                    sample.mean = sample.mean)
      stats
    }
    if (!any(class(survey.design) == "svyimputationList")) {
      sample.nobs <- lavaan::lavInspect(lavaan.fit, "nobs")
      stats <- get.stats.design(survey.design, sample.nobs)
    } else {stop("The survey design type 'svyimputationList' is currently not supported.")}

    stats
  }
  gamma.svy <- Gamma.svy(lavaan.fit, survey.design=des.rep)
  tcall <- lavaan::lavInspect(lavaan.fit, "call")
  tcall$data <- NULL
  tcall$sample.cov <- gamma.svy$sample.cov
  tcall$sample.mean <- gamma.svy$sample.mean
  tcall$sample.nobs <- n
  tcall$estimator <- estimator
  tcall$test <- test
  if (substr(estimator, 1, 2) == "ML") {  tcall$NACOV <- 	gamma.svy$Gamma }
  else if (substr(estimator, 1, 2) == "MLM") {tcall$NACOV <- 	gamma.svy$Gamma }
  main.res.new <- eval(as.call(tcall))

  #Calculate the BRR standard errors
  RR <- length(repwgtnames)
  runonce<-function(i){
    ## Step 1: generate weighted covariance matrix
    wsmpcov <- stats::cov.wt(x=data[,c(ovnames)], wt = data[,repwgtnames[i]], cor = FALSE, center = TRUE, method = "unbiased")
    ## Step 2: fit the model
    rep.res<- lavaan::sem(model=model, sample.cov = wsmpcov$cov, sample.mean = wsmpcov$center, sample.nobs = n, ...)
    return(rep.res@ParTable)
  }
  ## run parallel or not
  # this is from the boot function in package boot
  old_options <- options(); options(warn = -1)
  have_mc <- have_snow <- FALSE
  ncore <- ncore
  if (parallel != "no" && ncore > 1L) {
    if (parallel == "parallel") have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") have_snow <- TRUE
    if (!have_mc && !have_snow) ncore <- 1L
  }

  res <- if (ncore > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      #require(parallel)
      mclapply(seq_len(RR), runonce, mc.cores = ncore)
    } else if (have_snow) {
      #require(snow)
      #list(...) # evaluate any promises
      if (is.null(cl)) {
        cl <- makePSOCKcluster(rep("localhost", ncore))
        if(RNGkind()[1L] == "L'Ecuyer-CMRG")
          clusterSetRNGStream(cl)
        res <- parLapply(cl, seq_len(RR), runonce)
        stopCluster(cl)
        res
      } else parLapply(cl, seq_len(RR), runonce)
    }
  } else lapply(seq_len(RR), runonce)

  rep.est <-do.call(rbind, lapply(res, "[[", 'est'))
  rep.se <-do.call(rbind, lapply(res, "[[", 'se'))
  est.main <- main.res@ParTable$est
  tmat <- matrix(1, nrow = RR, ncol = length(est.main))%*%diag(est.main)
  se.BRR <- apply( (rep.est - tmat)^2, 2, function(x){
    sqrt(mean(x)/(1-fayfactor)^2 )} )
  z.BRR <- est.main/se.BRR
  p.BRR <- stats::pnorm(abs(z.BRR),mean=0,sd=1,lower.tail=F, log.p=F)*2.0 #two-sided p-values
  fit.BRR <- main.res.new
  fit.BRR@ParTable$se <- se.BRR
  fit.BRR@ParTable$z <- z.BRR
  fit.BRR@ParTable$p <- p.BRR
  fit.BRR@Options$se <- "BRR"
  fit.BRR@Options$population.num <- wn
  options(old_options)
  return(fit.BRR)
}

#' @name med.p.adjust
#' @title To adjust the p values for multimediation tests
#' @description
#' This function is used to adjust the p values when there are multiple mediators \cite{(Mai et al., 2019)}.
#' @param fit The model fit results of a model with multiple mediators. Note that it is a lavaan object.
#' @param med.eff A vector of labels. The labels should be of the mediation effects in the estimated model.
#' @param p.adj.method The method used to adjust for multiplicity (\code{'holm'} or \code{'hochberg'} or \code{'hommel'} or \code{'bonferroni'} or \code{'BH'} or \code{'BY'} or \code{'fdr'}).
#' Conservative method includes the Bonferroni correction ('bonferroni') in which the p-values are multiplied by the number of comparisons.
#' Less conservative corrections are also included by Holm (1979) ('holm'), Hochberg (1988) ('hochberg'), Hommel (1988) ('hommel'), Benjamini & Hochberg (1995) ('BH' or its alias 'fdr'), and Benjamini & Yekutieli (2001) ('BY'), respectively.
#' It is 'holm' by default. It is not required.
#' @return The adjusted p values along with the effect labels and original p values. It is a list.
#' @references Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289–300. DOI:10.2307/2346101
#' @references Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165–1188. DOI:10.1214/aos/1013699998
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6, 65–70.
#' @references Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika, 75, 383–386. DOI:10.1093/biomet/75.2.383
#' @references Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika, 75, 800–803. DOI:10.1093/biomet/75.4.800
#' @references Rosseel, Y. (2012). Lavaan: An R package for structural equation modeling and more. Version 0.5–12 (BETA). Journal of statistical software, 48(2), 1-36. DOI:10.18637/jss.v048.i02
#' @references Shaffer, J. P. (1995). Multiple hypothesis testing. Annual Review of Psychology, 46, 561–576.
#' @references Sarkar, S. (1998). Some probability inequalities for ordered MTP2 random variables: a proof of Simes conjecture. Annals of Statistics, 26, 494–504. DOI:10.1214/aos/1028144846
#' @examples
#' \dontshow{
#'  #Toy example for check:
#' R <- 20
#' wgtnames <- paste("repwgt", seq(0,R,by=1), sep="")
#' mwgtname=wgtnames[1]
#' repwgtnames=wgtnames[2:(R+1)]
#' model2 <- ' # outcome
#'               numcg ~ u0*1 + c*workban + b1*sp_adltban + b2*sp_kidsban
#'             # mediator
#'               sp_adltban ~ u1*1 + a1*workban
#'               sp_kidsban ~ u2*1 + a2*workban
#'             #covariance of residuals
#'               sp_adltban ~~ sp_kidsban
#'             # indirect effect (a*b)
#'               a1b1 := a1*b1
#'               a2b2 := a2*b2
#'             # total effect
#'               total := c + (a1*b1) + (a2*b2)
#'            '
#' fit.BRR2 <- med.fit.BRR(model=model2, data=MedData, mwgtname=mwgtname,
#'              repwgtnames=repwgtnames)
#' med.p.adjust(fit=fit.BRR2, med.eff=c('a1b1' , 'a2b2'))
#' }
#' \donttest{
#' R <- 160
#' wgtnames <- paste("repwgt", seq(0,R,by=1), sep="")
#' mwgtname=wgtnames[1]
#' repwgtnames=wgtnames[2:(R+1)]
#' fayfactor=0.5
#'
#' model2 <- ' # outcome
#'               numcg ~ u0*1 + c*workban + b1*sp_adltban + b2*sp_kidsban
#'             # mediator
#'               sp_adltban ~ u1*1 + a1*workban
#'               sp_kidsban ~ u2*1 + a2*workban
#'             #covariance of residuals
#'               sp_adltban ~~ sp_kidsban
#'             # indirect effect (a*b)
#'               a1b1 := a1*b1
#'               a2b2 := a2*b2
#'             # total effect
#'               total := c + (a1*b1) + (a2*b2)
#'            '
#' fit.BRR2 <- med.fit.BRR(model=model2, data=MedData, mwgtname=mwgtname,
#'              repwgtnames=repwgtnames, fayfactor, parallel='parallel', ncore=4)
#' temp <- med.p.adjust(fit=fit.BRR2, med.eff=c('a1b1' , 'a2b2'))
#' #
#' # Adjustment for multi mediation tests:
#' #
#' #      Effect          p Value      adj.p Value
#' #       a1b1      0.003667674      0.007335347
#' #       a2b2      0.217228711      0.217228711
#' #
#' # NOTE: 	 p Value adjustment method is holm
#' #
#' #########################################
#' # To catch the unformatted results:
#' temp
#' #
#' # $med.eff
#' # [1] "a1b1" "a2b2"
#' #
#' # $org.p.value
#' # [1] 0.003667674 0.217228711
#' #
#' # $adj.p.value
#' # [1] 0.007335347 0.217228711
#' }
##' @export
med.p.adjust <- function(fit=NULL, med.eff=NULL, p.adj.method=c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr')){
  if (missing(fit)) stop("A lavaan object/class is needed.")
  if (missing(med.eff)) stop("the labels of mediation effects are required.")
  if (!is.null(p.adj.method)) {p.adj.method <- match.arg(p.adj.method)} else {p.adj.method <- 'holm'}
  label <- fit@ParTable$label
  p <- fit@ParTable$p
  idx <-  match(label, med.eff)
  med.label <- label[!is.na(idx)]
  med.p <- p[!is.na(idx)]
  adj.p <- stats::p.adjust(med.p, method = p.adj.method)
  cat("\tAdjustment for multi mediation tests:", "\n\n", sep="")
  tabtitle <- sprintf("%15s  %15s  %15s \n", "Effect", "p Value", "adj.p Value" )
  cat(tabtitle)
  i=1
  for (i in 1:length(med.eff)) {
    tab0 <- sprintf("%15s  %15.9f  %15.9f \n", med.label[i], med.p[i], adj.p[i])
    cat(tab0)
  }

  if(!is.null(p.adj.method)) note <- paste("\t p Value adjustment method is ", p.adj.method, "\n", sep="")
  if (!is.null(note))
    cat("\n\tNOTE: ", note, "\n", sep = "")
  invisible(list(med.eff= med.label, org.p.value=med.p, adj.p.value=adj.p))
}

#' @name med.summary
#' @title To print the summary results of the mediation analysis
#' @description
#' This function is used to print the summary results of the mediation analysis with adjustment for multiplicity.
#' @param fit The model fit results of a mediation model. Note that it is a lavaan object.
#' @param med.eff A vector of labels. The labels should be of the mediation effects in the estimated model.
#' @param p.adj.method The method used to adjust for multiplicity (\code{'holm'} or \code{'hochberg'} or \code{'hommel'} or \code{'bonferroni'} or \code{'BH'} or \code{'BY'} or \code{'fdr'}).
#' Conservative method includes the Bonferroni correction ('bonferroni') in which the p-values are multiplied by the number of comparisons.
#' Less conservative corrections are also included by Holm (1979) ('holm'), Hochberg (1988) ('hochberg'), Hommel (1988) ('hommel'), Benjamini & Hochberg (1995) ('BH' or its alias 'fdr'), and Benjamini & Yekutieli (2001) ('BY'), respectively.
#' It is 'holm' by default. It is not required.
#' @return A list including the effect labels, estimates, standard errors, p values, and adjusted p values if there are more than one mediation effects.
#' @references Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289–300. DOI:10.2307/2346101
#' @references Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165–1188. DOI:10.1214/aos/1013699998
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6, 65–70.
#' @references Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika, 75, 383–386. DOI:10.1093/biomet/75.2.383
#' @references Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika, 75, 800–803. DOI:10.1093/biomet/75.4.800
#' @references Mai, Y., Ha, T., & Soulakova, J. N. (2019). Multimediation Method With Balanced Repeated Replications For Analysis Of Complex Surveys. Structural Equation Modeling: A Multidisciplinary Journal. DOI:10.1080/10705511.2018.1559065
#' @references Rosseel, Y. (2012). Lavaan: An R package for structural equation modeling and more. Version 0.5–12 (BETA). Journal of statistical software, 48(2), 1-36. DOI:10.18637/jss.v048.i02
#' @references Shaffer, J. P. (1995). Multiple hypothesis testing. Annual Review of Psychology, 46, 561–576.
#' @references Sarkar, S. (1998). Some probability inequalities for ordered MTP2 random variables: a proof of Simes conjecture. Annals of Statistics, 26, 494–504. DOI:10.1214/aos/1028144846
#' @examples
#' \dontshow{
#'  #Toy example for check:
#' R <- 20
#' wgtnames <- paste("repwgt", seq(0,R,by=1), sep="")
#' mwgtname=wgtnames[1]
#' repwgtnames=wgtnames[2:(R+1)]
#' model2 <- ' # outcome
#'               numcg ~ u0*1 + c*workban + b1*sp_adltban + b2*sp_kidsban
#'             # mediator
#'               sp_adltban ~ u1*1 + a1*workban
#'               sp_kidsban ~ u2*1 + a2*workban
#'             #covariance of residuals
#'               sp_adltban ~~ sp_kidsban
#'             # indirect effect (a*b)
#'               a1b1 := a1*b1
#'               a2b2 := a2*b2
#'             # total effect
#'               total := c + (a1*b1) + (a2*b2)
#'            '
#' fit.BRR2 <- med.fit.BRR(model=model2, data=MedData, mwgtname=mwgtname,
#'              repwgtnames=repwgtnames)
#' med.summary(fit=fit.BRR2, med.eff=c('a1b1' , 'a2b2'))
#' }
#' \donttest{
#' R <- 160
#' wgtnames <- paste("repwgt", seq(0,R,by=1), sep="")
#' mwgtname=wgtnames[1]
#' repwgtnames=wgtnames[2:(R+1)]
#' fayfactor=0.5
#'
#' model2 <- ' # outcome
#'               numcg ~ u0*1 + c*workban + b1*sp_adltban + b2*sp_kidsban
#'             # mediator
#'               sp_adltban ~ u1*1 + a1*workban
#'               sp_kidsban ~ u2*1 + a2*workban
#'             #covariance of residuals
#'               sp_adltban ~~ sp_kidsban
#'             # indirect effect (a*b)
#'               a1b1 := a1*b1
#'               a2b2 := a2*b2
#'             # total effect
#'               total := c + (a1*b1) + (a2*b2)
#'            '
#' fit.BRR2 <- med.fit.BRR(model=model2, data=MedData, mwgtname=mwgtname,
#'              repwgtnames=repwgtnames, fayfactor, parallel='parallel')
#' temp <- med.summary(fit=fit.BRR2, med.eff=c('a1b1' , 'a2b2'))
#' #
#' # MedSurvey 1.1.0
#' #
#' # Multimediation with Complex Survey Data:
#' #
#' #   Effect             Est.          BRR SE.          p Value      adj.p Value
#' #
#' #   a1b1     -0.017475544      0.006014820      0.003667674      0.007335347
#' #   a2b2     -0.007244189      0.005870823      0.217228711      0.217228711
#' #
#' # NOTE:
#' #   p Value adjustment method is holm
#' #   Standard errors type is BRR SE.
#' #
#' #
#' ######################################
#' # To catch the unformatted results:
#' temp
#' #
#' # $med.label
#' # [1] "a1b1" "a2b2"
#' #
#' # $med.est
#' # [1] -0.017475544 -0.007244189
#' #
#' # $med.se
#' # [1] 0.006014820 0.005870823
#' #
#' # $org.p.value
#' # [1] 0.003667674 0.217228711
#' #
#' # $adj.p.value
#' # [1] 0.007335347 0.217228711
#' #
#' # $se.type
#' # [1] "BRR SE."
#' #
#' # $p.adj.method
#' # [1] "holm"
#' #
#' }
#'
##' @export
med.summary <- function(fit=NULL, med.eff=NULL, p.adj.method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
  if (missing(fit)) stop("A lavaan object/class is needed.")
  if (missing(med.eff)) stop("the labels of mediation effects are required.")
  if (!is.null(p.adj.method)) {p.adj.method<-match.arg(p.adj.method)} else {p.adj.method<-'holm'}
  note1 <- paste("\t p Value adjustment method is ", p.adj.method, sep="")
  se.type <- paste(fit@Options$se," SE.",sep="")
  note2 <- paste("\t Standard errors type is ", se.type, '\n', sep="")
  note <- NULL
  cat('\tMedSurvey 1.1.0 \n\n');
  cat("\tMultimediation with Complex Survey Data:","\n\n", sep="")
  x <- fit@ParTable
  idx <- match(x$label, med.eff)
  med.label <- x$label[!is.na(idx)]
  med.est <- x$est[!is.na(idx)]
  med.se <- x$se[!is.na(idx)]
  med.p <- x$p[!is.na(idx)]
  if (length(med.eff)>=2){
    adj.p <- stats::p.adjust(med.p, method = p.adj.method)
    tabtitle <- sprintf("%15s  %15s  %15s  %15s  %15s \n", "Effect", "Est.", se.type, "p Value", "adj.p Value" )
    cat(tabtitle)
    i=1
    for (i in 1:length(med.eff)) {
      tab0 <- sprintf("%15s  %15.9f  %15.9f  %15.9f  %15.9f \n", med.label[i],med.est[i], med.se[i], med.p[i], adj.p[i])
      cat(tab0)
    }
    note <- paste(note, note1 ,note2, sep="\n")
    if (!is.null(note)) cat("\n\tNOTE: ", note, sep = "")
  } else {
    p.adj.method <- NULL
    adj.p <- NULL
    tabtitle <- sprintf("%15s  %15s  %15s  %15s \n", "Effect", "Est.", se.type, "p Value")
    cat(tabtitle)
    i=1
    for (i in 1:length(med.eff)) {
      tab0 <- sprintf("%15s  %15.9f  %15.9f  %15.9f \n", med.label[i],med.est[i], med.se[i], med.p[i])
      cat(tab0)
    }
    note <- paste(note, note2, sep="\n")
    if (!is.null(note))
      cat("\n\tNOTE: ", note, '\n', sep = "")
  }
  invisible(x)
  res <- (list(med.label = med.label, med.est=med.est, med.se=med.se, org.p.value=med.p, adj.p.value=adj.p, se.type=se.type, p.adj.method=p.adj.method))
  invisible(res)
}
