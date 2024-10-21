#' Estimate ICC for nominal or ordinal categorical response data
#'
#' @param cid Cluster id variable.
#' @param y Categorical response variable.
#' @param data Dataframe containing 'cid' and 'y'.
#' @param alpha Significance level for confidence interval computation. Default is 0.05.
#' @param method Method used to estimate categorical ICC. A single method or multiple methods can be specified. Default is both resampling and moments estimators. See iccmult::iccmulti for more details.
#' @param binmethod Method used to estimate binary ICC. A single or multiple methods can be specified. By default all 16 methods are returned. See full details in ICCbin::iccbin().
#' @param ci.type Type of confidence interval to be computed for binary ICC. By default, all 5 types will be returned See full details in ICCbin::iccbin() for more.
#' @param kappa Value of Kappa to be used in computing Stabilized ICC when the binary response method 'stab' is chosen. Default value is 0.45.
#' @param nAGQ An integer scaler, as in lme4::glmer(), denoting the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood. Used when the binary response method 'lin' is chosen. Default value is 1.
#' @param M Number of Monte Carlo replicates used in binary ICC computation method 'sim'. Default is 1000.
#' @param nowarnings Flag to turn off estimation warnings. Default is False.
#'
#' @returns Data frame or list of data frames with single column estimate of ICC, se(ICC), and lower and upper CI bounds.
#'
#' @importFrom stats na.omit qnorm rnorm runif
#' @importFrom gtools permutations
#' @importFrom dirmult weirMoM
#' @import ICCbin
#' @import lme4
#' @import dirmult
#' @export
#'
#' @examples
#' iccdat4 <- rccat(rho=0.15, prop=c(0.15,0.25,0.20,0.40), noc=10, csize=25)
#' iccmulti(cid=cid, y=y, data=iccdat4)
#' iccdat3 <- rccat(rho=0.10, prop=c(0.30,0.25,0.45), noc=15, csize=50)
#' iccmulti(cid=cid, y=y, data=iccdat3)
iccmulti = function(cid, y, data, alpha=0.05, method=c("rm","mom"),
                    binmethod = c("aov", "aovs", "keq", "kpr", "keqs",
                                  "kprs", "stab", "ub", "fc", "mak", "peq",
                                  "pgp", "ppr", "rm", "lin", "sim"),
                    ci.type = c("aov", "wal", "fc", "peq", "rm"),
                    kappa = 0.45, nAGQ = 1, M = 1000, nowarnings=FALSE){
  # ICC estimation function that runs on (any) number of categories (2+)
  # cid       = Column name indicating cluster id in dataframe data
  # y         = Column name indicating categorical response in data frame `data'
  # data      = Dataframe containing `cid' and `y'
  # alpha     = Significance level for confidence interval computation. Default
  #             value is 0.05.
  # method    = Method used to compute categorical ICC. A single method or
  #             multiple methods can be specified. Default is both resampling
  #             & moments methods. See Details of iccmult::iccmulti() for more.
  # binmethod = Method used to compute binary ICC. A single or multiple methods
  #             can be specified. By default, all 16 methods will be used. See
  #             details of ICCbin::iccbin() for more.
  # ci.type   = Type of confidence interval to be computed for binary ICC. By
  #             default, all 5 types will be reported. See details of
  #             ICCbin::iccbin() for more.
  # kappa     = Value of Kappa to be used in computing Stabilized ICC when the
  #             binary response method stab is chosen. Default value is 0.45.
  # nAGQ      = An integer scaler, as in glmer function of package lme4,
  #             denoting the number of points per axis for evaluating the
  #             adaptive Gauss-Hermite approximation to the log-likelihood.
  #             Used when the binary response method lin is chosen. Default
  #             value is 1.
  # M         = Number of Monte Carlo replicates used in ICC computation method
  #             sim. Default is 1000.
  # nowarnings = Flag to turn off estimation warnings. Default is False.

  ##~~~ NOTE: ICCbin::iccbin may return an error between lmer and Matrix package
  ##~~~ when method=lin or sim. Run the following install call to avoid this
  ##~~~ error:
  ##~~~ > install.packages("lme4", type = "source")

  requireNamespace("ICCbin","lme4","dirmult","gtools")

  CALL <- match.call()
  ic <- list(cid = substitute(cid), y = substitute(y))
  if (is.character(ic$y)) {
    if (missing(data))
      stop("Supply either the unqouted name of an object containing 'y' or supply both 'data' and then 'y' as an unquoted column name to 'data'")
    ic$y <- eval(as.name(y), data, parent.frame())
  }
  if (is.name(ic$y))
    ic$y <- eval(ic$y, data, parent.frame())
  if (is.call(ic$y))
    ic$y <- eval(ic$y, data, parent.frame())
  if (is.character(ic$y))
    ic$y <- eval(as.name(ic$y), data, parent.frame())
  if (is.character(ic$cid)) {
    if (missing(data))
      stop("Supply either the unquoted name of an object containing 'cid' or supply both 'data' and then 'cid' as an unquoted column name to 'data'")
    ic$cid <- eval(as.name(cid), data, parent.frame())
  }
  if (is.name(ic$cid))
    ic$cid <- eval(ic$cid, data, parent.frame())
  if (is.call(ic$cid))
    ic$cid <- eval(ic$cid, data, parent.frame())
  if (is.character(ic$cid) && length(ic$cid) == 1)
    ic$cid <- eval(as.name(ic$cid), data, parent.frame())
  dt <- data.frame(ic)
  dt <- stats::na.omit(dt)
  k <- length(unique(dt$cid))
  if (!is.null(attributes(dt)$na.action)) {
    warning("NAs removed from data rows:\n", paste0(unclass(attributes(dt)$na.action),collapse=","), "\n")
  }
  if (!is.factor(dt$cid)) {
    warning("'cid' has been coerced to a factor")
    dt$cid <- as.factor(dt$cid)
  }
  else {
    if (length(levels(dt$cid)) > k) {
      dt$x <- factor(as.character(dt$cid), levels = unique(dt$cid))
      warning("Missing levels of 'cid' have been removed")
    }
  }

  ### Function to capture unique warnings
  saveWarnings <- function(expr) {
    warnings <- NULL
    wHandler <- function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
    result <- withCallingHandlers(expr, warning = wHandler)
    list(result = result, warnings = warnings)
  }

  # Order records by cluster
  cordata = dt[order(dt$cid),]
  # convert factor outcome
  if(is.factor(cordata$y)) cordata$y = as.numeric(cordata$y)
  # Pull unique variable values
  yvals = sort(unique(cordata$y))

  ncat = length(yvals)

  rho = se = NULL

  if(ncat>=3){
    if("rm" %in% method){
      warns = NULL

      for(i in 1:ncat){
        ysub = NULL

        # ICC for component i vs others
        cordatap = cordata;
        cordatap$ysub = ifelse(cordatap$y==yvals[i],1,0)

        # Save iccbin results and warnings
        iball <- saveWarnings({
          ICCbin::iccbin(cid, ysub, data=cordatap, method="rm", ci.type="rm", alpha=alpha)
        })
        # iccbin warnings
        warns=c(warns,iball$warnings)
        # iccbin estimates
        ib = iball$result

        # Check if ICC estimable:
        if(ib$estimates[,"ICC"]=="-"){ rho = c(rho,NA) }else{ rho = c(rho,ib$estimates[,"ICC"])}

        # Check if CI estimable:
        sumci = sum(ib$ci[,c("LowerCI","UpperCI")]==rep("-",2))
        if(sumci>0){
          ci=c(0,0);
          if(ib$ci[,"LowerCI"]=="-"){ ci[1] = NA }else{ ci[1] = ib$ci[,"LowerCI"] }
          if(ib$ci[,"UpperCI"]=="-"){ ci[2] = NA }else{ ci[2] = ib$ci[,"UpperCI"] }
          se = NA
        }
        else{
          ci = as.numeric(ib$ci[,c("LowerCI","UpperCI")])
          se = c(se,(ci[2]-ci[1])/(2*qnorm(1-alpha/2)))
        }
      }

      if(ncat>3){
        for(j in 2:floor(ncat/2)){
          combos = gtools::permutations(2,ncat,c(0,1),rep=T)
          combos = combos[which(rowSums(combos)==j),]
          if(j==floor(ncat/2) & (ncat%%2)==0) combos = combos[(nrow(combos)/2+1):nrow(combos),]

          for(k in 1:nrow(combos)){
            pick = c(1:ncat)[which(combos[k,]>0)]
            # ICC for (2+)-categories vs others
            cordatap = cordata;
            cordatap$ysub = ifelse(cordatap$y%in%yvals[pick],1,0)

            # Save iccbin results and warnings
            iball <- saveWarnings({
              ICCbin::iccbin(cid, ysub, data=cordatap, method="rm", ci.type="rm", alpha=alpha)
            })
            # iccbin warnings
            warns=c(warns,iball$warnings)
            # iccbin estimates
            ib = iball$result

            # Check if ICC estimable:
            if(ib$estimates[,"ICC"]=="-"){ rho = c(rho,NA) }else{ rho = c(rho,ib$estimates[,"ICC"])}

            # Check if CI estimable:
            sumci = sum(ib$ci[,c("LowerCI","UpperCI")]==rep("-",2))

            if(sumci>0){
              ci=c(0,0);
              ci[which(ib$ci[,c("LowerCI","UpperCI")]=="-")]=NA;
              ci[which(ib$ci[,c("LowerCI","UpperCI")]!="-")] = ib$ci[,c("LowerCI","UpperCI")][which(ib$ci[,c("LowerCI","UpperCI")]!="-")]
              se = NA
            }
            else{
              ci = as.numeric(ib$ci[,c("LowerCI","UpperCI")])
              se = c(se,(ci[2]-ci[1])/(2*qnorm(1-alpha/2)))
            }
          }
        }
      }

      # # # # Change capitalization of warnings # # # #
      make_lower <- function(text, exception) {
        words <- unlist(strsplit(text, " "))
        words <- ifelse(words %in% exception | words == words[1],
                        words, tolower(words))
        paste(words, collapse = " ")
      }
      multiwarns = unname( sapply(unique(warns),function(x) make_lower(x,
                                                                       c("ICC","'Resampling'","'Resampling"))) )
      multiwarns = paste0(multiwarns," for at least one categorical calculation")
      if(nowarnings==FALSE & is.null(warns)==FALSE){
        # Print only one of ICCbin calculation warnings
        for(w in 1:length(multiwarns)){
          warning(multiwarns[w])
        }
      }


      ################################################
      ## Overall X-category resampling ICC estimate ##
      rhom = mean(rho)

      # Variance & 95% CI of overall resampling estimate
      vrhom = 1/length(rho)*sum(se^2)
      rhomci = rhom+c(-1,1)*qnorm(1-alpha/2)*sqrt(vrhom)

      # Check if CIs <0 or >1, reset to 0 or 1
      if(is.na(rhomci[1])==F & is.na(rhomci[2])==F){
        if(rhomci[1] < 0 | rhomci[2] > 1) {
          if(rhomci[1]<0) rhomci[1] = 0
          if(rhomci[2]>1) rhomci[2] = 1
          warning("One or both of 'Resampling based' method confidence limits fell outside of [0, 1]")
        }
      }
      else{
        warning("Resampling based confidence interval for categorical ICC is not estimable")
      }

      # Check if rhom is missing
      if(is.na(rhom)==TRUE){
        warning("Categorical ICC not estimable by 'Resampling' method")
      }
      else{
        # If estimate is <0 or >1 set to missing
        if (rhom < 0 | rhom > 1) {
          rhom = "-"
          warning("Categorical ICC not estimable by 'Resampling' method")
          # warning("Resampling based confidence interval for ICC is not estimable")
        }
      }

      outrm = data.frame("Avg.Resampling.Rho"=c(rhom,sqrt(vrhom),rhomci))
      rownames(outrm) = c("ICC","se(ICC)","LCL","UCL")

    }
    if("mom" %in% method){
        if (!requireNamespace("dirmult", quietly = TRUE)) {
          stop("Package 'dirmult' is needed for method 'mom'. Please install it.",
               call. = FALSE)
        }
      #################################
      ## X-category MoM ICC estimate ##
      moms = dirmult::weirMoM(as.matrix(table(cordata$cid,cordata$y)),se=T)
      rhomom = moms$theta
      rhomomse = moms$se

      # 95% CI of overall resampling estimate
      rhomomci = rhomom+c(-1,1)*qnorm(1-alpha/2)*rhomomse

      # Check if CIs <0 or >1, reset to 0 or 1
      if (rhomomci[1] < 0 | rhomomci[2] > 1) {
        if(rhomomci[1]<0) rhomomci[1] = 0
        if(rhomomci[2]>1) rhomomci[2] = 1
        warning("One or both of 'Moments' method confidence limits fell outside of [0, 1]")
      }

      # If estimate is <0 or >1 set to missing
      if (rhomom < 0 | rhomom > 1) {
        rhomom = "-"
        outmom = data.frame("MomentEst.Rho"=c(rhomom,rhomomse,rhomomci))
        warning("ICC not estimable by 'Moments' Method")
      }
      else{
        outmom = data.frame("MomentEst.Rho"=c(rhomom,rhomomse,rhomomci))
      }
      rownames(outmom) = c("ICC","se(ICC)","LCL","UCL")
    }

    ###########################################################################
    ## Output final results ##
    ###########################################################################
    if(length(method)==1){
      if(method=="rm") out=outrm
      else out=outmom
    }
    else{
      out = list(rm=outrm,mom=outmom)
    }
  }
  else{
    ## Redefine binary outcome ##
    dt$ysub = ifelse(dt$y==yvals[2],1,0)

    ## Compute binary ICC estimate(s) ##
    ib = ICCbin::iccbin(cid=cid, y=ysub, data=dt, method=binmethod, ci.type=ci.type,
                        alpha=alpha, kappa=kappa, nAGQ=nAGQ, M=M)

    ## Output final results ##
    out=ib
  }

  return(out)
}
