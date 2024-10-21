#' Generate Correlated Clustered Categorical Data
#'
#' @param rho Numeric value between 0 and 1 of the desired ICC value.
#' @param prop Numeric vector of each response category's probability, each taking value between 0 and 1.
#' @param prvar Numeric value or vector of values between 0 and 1 denoting percent variation in each assumed event rate. Default is 0.
#' @param noc Numeric value of number of clusters to be generated.
#' @param csize Numeric value of desired cluster size.
#' @param csvar Numeric value between 0 and 1 denoting percent variation in cluster sizes. Default is 0.
#' @param allevtcl Logical value specifying whether all clusters must have all categories. Default is True.
#' @param drawn Maximum number of attempts to apply variation to event probabilities.
#' @param nowarnings Flag to turn off warnings. Default is False.
#'
#' @returns Dataframe with two columns, a column identifier 'cid' and categorical response 'y', and one row for each observation within each cluster
#' @export
#'
#' @examples
#' rccat(rho=0.2, prop=c(0.2, 0.3, 0.5), prvar=0, noc=5, csize=20, csvar=0.2)
#' rccat(rho=0.1, prop=c(0.2, 0.4, 0.3, 0.1), prvar=0.10, noc=30, csize=40, csvar=0)
rccat = function(rho, prop, prvar=0, noc, csize, csvar=0, allevtcl=TRUE,
                         drawn=10, nowarnings=FALSE){
  ### Code to generate clustered (correlated) multicategory data
  ### based on: Biswas A (2004). Generating Correlated Ordinal Random Samples. Stat Probs & Letters
  ### Code to generate clustered (correlated) multicategory data ###
  # rho        = numeric value between 0 and 1 of the assumed ICC value
  # prop       = numeric vector of each response category's probability, each
  #               taking value between 0 and 1
  # prvar      = numeric value or vector of values between 0 and 1 denoting
  #               percent variation in each assumed event rate, default is 0.
  # noc        = numeric value of number of clusters to be generated
  # csize      = numeric value of desired cluster size
  # csvar      = numeric value between 0 and 1 denoting percent variation in
  #               cluster sizes, default is 0
  # allevtcl   = logical value specifying whether all clusters must have all
  #               categories, default is TRUE
  # drawn      = maximum number of attempts to apply variation to event
  #               probabilities
  # nowarnings = flag to turn off warnings. Default is False.

  if(noc<1){
    stop("Number of clusters must be at least 2.")
  }
  if(csvar<0 | csvar>1){
    stop("Cluster variation must be within [0,1].")
  }
  if(sum(sapply(prvar, function(x) (x<0 | x>1)))>0){
    stop("Event rate variation must be within [0,1].")
  }
  if(length(prop)>1 & sum(prop)!=1){
    stop("Categorical event rate probabilities must sum to 1.")
  }

  ## If prop is a single number, assume event rate for binary response
  if(length(prop)==1) prop = c(1-prop,prop)


  ## set correlation rho = (b/(sum(a_u) + b))^2 --> b = sqrt(rho)/(1-sqrt(rho))
  b = sqrt(rho)/(1-sqrt(rho))
  out = NULL

  ## Initialize vector of cluster ids
  cluster <- c()

  ## Initialize categorical outcome vector
  y <- c()


  ## Save warnings
  rcwarn = NULL


  ################################################
  ### Iterate following steps for each cluster ###
  ################################################
  for (i in 1:noc) {
    ###### Induce variation in cluster size ######

    ## Define min cluster size as (1-csvar)*csize or at least # categories
    min_csize <- ifelse((csize - round(csize * csvar)) >=
                          length(prop), csize - round(csize * csvar), length(prop))

    ## Define cluster size of current [i] cluster
    csizen <- abs(round(csize + (csize * csvar) * rnorm(1)))

    ## If too small, draw again
    while (csizen < min_csize) {
      csizen <- abs(round(csize + (csize * csvar) * rnorm(1)))
    }


    ###### Induce variation in event probabilities ######

    ##~~ If all event probs have same variation percentage ~~##
    if(length(prvar)==1){

      # Flag for probability (re)draws
      drawp=drawcount=1

      ## Define min event proportion: (1-prvar)*prop or at least 0
      min_pr <- sapply(prop, function(a)
        ifelse((a - a * prvar) >= 0, a - a * prvar, 0))

      ## Define max event proportion: (1+prvar)*prop or no more than 1
      max_pr <- sapply(prop, function(a)
        ifelse((a + a * prvar) <= 1, a + a * prvar, 1))

      ## Create new event probability vector
      prn = prop

      while(drawp==1 & drawcount<drawn){
        # Update all but one probability to ensure all sum to 1
        updatep = sample(1:length(prop),length(prop)-1)
        for(k in updatep){
          ## Define event proportion within current [i] cluster
          prn[k] <- abs(prop[k] + (prop[k] * prvar) * rnorm(1))
          ## If too small or too big, draw again
          while (prn[k] < min_pr[k] | prn[k] > max_pr[k]) {
            prn[k] <- abs(prop[k] + (prop[k] * prvar) * rnorm(1))
          }
        }
        # Update final category
        prn[-updatep] = 1-sum(prn[updatep])
        if(min_pr[-updatep]<=prn[-updatep] & prn[-updatep]<=max_pr[-updatep]){
          drawp=0
        }
        drawcount = drawcount+1
      }

    }else{
      # Flag for probability (re)draws
      drawp=drawcount=1

      # Update all but one probability to ensure all sum to 1
      updatep = sample(1:length(prop),length(prop)-1)
      ##~~ event probs with distinct variation percentages ~~##
      prn = prop
      while(drawp==1 & drawcount<drawn){
        for(k in updatep){
          ## Define min event proportion: (1-prvar)*prop or at least 0
          min_pr <- ifelse((prop[k] - prop[k] * prvar[k]) >= 0, prop[k] -
                             prop[k] * prvar[k], 0)
          ## Define max event proportion: (1+prvar)*prop or no more than 1
          max_pr <- ifelse((prop[k] + prop[k] * prvar[k]) <= 1, prop[k] +
                             prop[k] * prvar[k], 1)
          ## Define event proportion within current [i] cluster
          prn[k] <- abs(prop[k] + (prop[k] * prvar[k]) * rnorm(1))
          ## If too small or too big, draw again
          while (prn[k] < min_pr | prn[k] > max_pr) {
            prn[k] <- abs(prop[k] + (prop[k] * prvar[k]) * rnorm(1))
          }
        }

        ## Update final category
        ## Define min event proportion: (1-prvar)*prop or at least 0
        min_pr <- ifelse((prop[-updatep] - prop[-updatep] * prvar[-updatep])>=0,
                         prop[-updatep] - prop[-updatep] * prvar[-updatep], 0)
        ## Define max event proportion: (1+prvar)*prop or no more than 1
        max_pr <- ifelse((prop[-updatep] + prop[-updatep] * prvar[-updatep])<=1,
                         prop[-updatep] + prop[-updatep] * prvar[-updatep], 1)
        # Define final category probability
        prn[-updatep] = 1-sum(prn[updatep])
        # Check it falls within its min & max pr
        if(min_pr<= prn[-updatep] & prn[-updatep]<= max_pr){ drawp=0}

        drawcount = drawcount+1
      }
    }
    if(drawp==1 & prvar!=0){
      rcwarn = c(rcwarn,"Limit reached when varying event probabilities. Probabilities reset to original values with 0% variation.")
      # Reset probabilities
      prn = prop
    }


    #################################################################
    ######                Biswas A3 algorithm                  ######
    ###### (see pseudocode in Chakraborty Comm. in Stats 2023) ######
    #################################################################

    #######################################
    ### Step 1: Set prob distribution of random effect Y0 to desired prob distn
    # Calculate cumulative probability distribution
    cdf0 = c(0,cumsum(prn))

    #######################################
    ### Step 2: draw random effect Y0 ###
    # Draw a random number between [0,1]
    r0 = runif(1)
    # For k'=0,1,...,length(prop)-1 set y0 = k' if cdf0[k-1] < r0 <= cdf0[k], k' = k-1
    #   where 0,...length(prop)-1 are (the desired values of final outcome y - 1)
    y0 = seq(0,length(prn)-1,1)[sapply(1:(length(cdf0)-1),
                                       function(x) (cdf0[x]<r0 & r0<=cdf0[x+1]))]

    #######################################
    ### Step 3: Update probability distribution and draw outcome values ###
    # Initialize update prob distn to current prob
    prnew = prn

    ## Update the probability corresponding to value y0:
    # The prob corresponding to category [k] that y0 equals is now =
    #     (Pr(k)+b)/(sum(prop)+b)
    prnew[which(seq(0,length(prn)-1,1)==y0)] =
      ( prn[which(seq(0,length(prn)-1,1)==y0)] + b ) / ( sum(prn)+b )

    # Update the remaining probabilities
    # Those prop not corresponding to category [k] = y0 are now =
    #     Pr([-k])/(sum(prop)+b)
    prnew[which(seq(0,length(prn)-1,1)!=y0)] =
      sapply(prn[which(seq(0,length(prn)-1,1)!=y0)],
             function(x) x / ( sum(prn)+b ) )

    # Define cumulative prob distn of multicategory outcome for cluster [i]
    cdfi = c(0,cumsum(prnew))


    #################################
    ### Draw categorical outcomes ###
    #################################
    yi = clusteri = NULL
    for(j in 1:csizen){
      rj = runif(1)
      # Set individual j to category [k] if cdfi[k-1] < rj <= cdfi[k]
      yij = seq(0,length(prn)-1,1)[sapply(1:(length(cdfi)-1),
                                          function(x) (cdfi[x]<rj & rj<=cdfi[x+1]))]
      yi = c(yi,yij)
      # Repeat cluster ID for all individual values
      clusteri=c(clusteri,i)
    }
    ### If must have all 3 categories in all clusters
    if(allevtcl==TRUE & length(unique(yi))<length(prn)){
      while(length(unique(yi))<length(prn)){
        # Drop previous pulls of y
        yi = NULL

        for(j in 1:csizen){
          rj = runif(1)
          # Set individual j to category [k] if cdfi[k-1] < rj <= cdfi[k]
          yij = seq(0,length(prn)-1,1)[sapply(1:(length(cdfi)-1),
                                              function(x) (cdfi[x]<rj & rj<=cdfi[x+1]))]
          yi = c(yi,yij)
        }
      }
    }
    y = c(y,yi)
    cluster = c(cluster,clusteri)
  }

  # Print only one of ICCbin calculation warnings
  if(is.null(rcwarn)==FALSE & nowarnings==FALSE){
    rcmultiwarns = unique(rcwarn)
    for(w in 1:length(rcmultiwarns)){
      warning(rcmultiwarns[w])
    }
  }

  out = data.frame(cid=as.factor(cluster),y=y)
  if(length(prop)>2) out$y = out$y+1
  return(out)
}
