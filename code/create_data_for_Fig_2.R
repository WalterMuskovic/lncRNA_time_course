## Figure 2 - The effects of transcript stability
# We look at 7 genes with rapid induction dynamics; FOS, HES1, FOSB, CTGF, SRF, DSTN and TPM4.
# For these genes we will use a simple model of transcription:

# $\frac{dM}{dt} = \beta P(t) - \alpha M(t)$

# to infer the half-lives of these transcripts.
  
# To assist in visualizing the pre-mRNA expression data we will also fit impulse models to the
# expression estimates obtained from RNA-seq. The impulse model is described by
# [Chechik and Koller](https://doi.org/10.1089/cmb.2008.13TT). They describe a a parametric 
# model that is a product of two sigmoid functions. Their original model has six free parameters;
# h0, h1, h2, t1, t2, and lambda. Three amplitude (h) parameters determine the initial amplitude
# (h0), the peak amplitude (h1), and the steady state amplitude (h2) of the response. The onset
# time t1 is the time of first transition (where rise or fall is maximal) and the offset t2 is
# the time of second transition. Lambda is the slope parameter. As the authors point out, this
# model is easily generalizable and we modify it slightly here to include lambda1 and lambda2,
# such that the two transitions may have different slopes. So we get:

#$\frac{1}{h_{1}} \left[ h_{0} + (h_{1} -  h_{0})(\frac{1}{1 + e^{-\lambda_{1}(t-t1)}}) \right]\left[ h_{2} + (h_{1} -  h_{2})(\frac{1}{1 + e^{\lambda2(t-t2)}}) \right]$

# To estimate the seven free parameters we'll use the `nlsLM()` function from the `minpack.lm`
# R package. Note that we'll also use the six-parameter version and a nine-parameter model that
# allows up to three transitions.

# Load R packages
library(DESeq2)
library(tidyverse)
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
library(minpack.lm)
library(deSolve)

# If not already done, get data for figure 2
if(!file.exists("data_track/fig_2_data.rds")){
  
  # Create a data frame to store the data needed for figure 2
  fig_2_data <- data.frame(gene_name = c("FOS", "HES1", "FOSB", "CTGF", "SRF", "DSTN", "TPM4"),
                           gene_id = c("ENSG00000170345.9", "ENSG00000114315.3", "ENSG00000125740.13", "ENSG00000118523.5", "ENSG00000112658.7", "ENSG00000125868.15", "ENSG00000167460.15"),
                           transcript_id = c("ENST00000303562.8", "ENST00000232424.3", "ENST00000353609.7", "ENST00000367976.3", "ENST00000265354.5", "ENST00000246069.11", "ENST00000643579.1"),
                           impulse_model_p = c("2", "2", "2b", "3", "2b", "3", "3"),
                           mRNA_expression = NA, 
                           last_10Kb = NA,
                           impulse_model_fit = NA,
                           transcription_model_fit = NA,
                           alpha = NA)
  
  # Import gene counts
  gene_counts <- readRDS("data/t98g_filtered_coords_clust.rds")
  
  # Restrict to genes of interest
  gene_counts <- gene_counts[match(fig_2_data$gene_name,  gene_counts$gene_name),]
  
  # Add to fig_2_data
  fig_2_data$mRNA_expression <- gene_counts$exon
  fig_2_data$last_10Kb <- gene_counts$last_10Kb
  
  ## We will also get impulse model fits to assist in the visualisation of the pre-mRNA data
  # We will be using the  nlsLM function from the minpack.lm package to do the heavy lifting here
  
  # Define  impulse models 
  # Define 6-parameter impulse model function that allows two transitions, where both transitions have the same slope, lambda.
  impulse2 <- function(h0, h1, h2, lambda, t, t1, t2) { (1/h1)*((h0+(h1-h0)*(1/(1+exp(-lambda*(t-t1)))))*(h2+(h1-h2)*(1/(1+exp(lambda*(t-t2)))))) }
  # Define 7-parameter impulse model function that allows two transitions, where the transitions can have different slopes, defined by lambda1 and lambda2.
  impulse2b <- function(h0, h1, h2, lambda1, lambda2, t, t1, t2) { (1/h1)*(h0+(h1-h0)*(1/(1+exp(-lambda1*(t-t1)))))*(h2+(h1-h2)*(1/(1+exp(lambda2*(t-t2))))) }
  # Define a 9 parameter impulse model that allows up to three transitions
  impulse3 <- function(h0, h1, h2, h3, lambda1, lambda2, t, t1, t2, t3) { (1/h1)*((h0+(h1-h0)*(1/(1+exp(-lambda1*(t-t1)))))*(h2+(h1-h2)*(1/(1+exp(-lambda1*(t-t2)))))*(h3+(h2-h3)*(1/(1+exp(-lambda2*(t-t3)))))) }
  
  # Define time
  t<-seq(0,400,10)
  
  # Define function that will fit the impulse2 model and return the best fit after 100 iterations
  fit_impulse2 <- function(input_time_series, t){
    temp <- data.frame(t=t, M=input_time_series)
    best_NRMSD <- 1E6
    best_fit <- NA
    for(i in 1:100){
      print(str_glue('nlsLM fit - attempt {i} of 100'))
      # Set fit to 0 to start with
      impulse_fit <- 0
      # Continue trying to the model to the data until a fit is found
      while(is.numeric(impulse_fit)){
        impulse_fit <- tryCatch({
          nlsLM(M ~ impulse2(h0, h1, h2, lambda, t, t1, t2), 
                data = temp, start = c(h0=runif(1, min=min(input_time_series), max = max(input_time_series)), # the initial expression value
                                       h1=runif(1, min=min(input_time_series), max = max(input_time_series)), # the peak expression value
                                       h2=runif(1, min=min(input_time_series), max = max(input_time_series)), # the steady-state expression value
                                       lambda=runif(n=1, min=0, max=1), # Slope of the transitions
                                       t1=sample(1:max(t),1), # onset time is the time of the first transition (where rise or fall is maximal)
                                       t2=sample(1:max(t),1)), # offset time is the time of the second transition
                control = list(maxiter = 1024,warnOnly=TRUE), trace = F)
        }, error = function(e) { # if nlsLM can't converge, ouput zero
          c(0)
        }, finally = {
          
        })
      }
      # calculate the NRMSD
      current_NRMSD <- sqrt(sum((residuals(impulse_fit))^2)/length(residuals(impulse_fit)))/diff(range(temp$M))
      # Check if current fit is more optimal
      if(current_NRMSD < best_NRMSD){ 
        best_NRMSD <- current_NRMSD
        best_fit <- impulse_fit
      }
    }
    return(best_fit)
  }
  
  # Define function that will fit the impulse2b model and return the best fit after 100 iterations
  fit_impulse2b <- function(input_time_series, t){
    temp <- data.frame(t=t, M=input_time_series)
    best_NRMSD <- 1E6
    best_fit <- NA
    for(i in 1:100){
      print(str_glue('nlsLM fit - attempt {i} of 100'))
      # Set fit to 0 to start with
      impulse_fit <- 0
      # Continue trying to the model to the data until a fit is found
      while(is.numeric(impulse_fit)){
        impulse_fit <- tryCatch({
          nlsLM(M ~ impulse2b(h0, h1, h2, lambda1,lambda2, t, t1, t2), 
                data = temp, start = c(h0=runif(1, min=min(input_time_series), max = max(input_time_series)), # the initial expression value
                                       h1=runif(1, min=min(input_time_series), max = max(input_time_series)), # the peak expression value
                                       h2=runif(1, min=min(input_time_series), max = max(input_time_series)), # the steady-state expression value
                                       lambda1=runif(n=1, min=0, max=1), # Slope of the first transition
                                       lambda2=runif(n=1, min=0, max=1), # Slope of the second transition
                                       t1=sample(1:max(t),1), # onset time is the time of the first transition (where rise or fall is maximal)
                                       t2=sample(1:max(t),1)), # offset time is the time of the second transition
                control = list(maxiter = 1024,warnOnly=TRUE), trace = F)
        }, error = function(e) { # if nlsLM can't converge, ouput zero
          c(0)
        }, finally = {
          
        })
      }
      # calculate the NRMSD
      current_NRMSD <- sqrt(sum((residuals(impulse_fit))^2)/length(residuals(impulse_fit)))/diff(range(temp$M))
      # Check if current fit is more optimal
      if(current_NRMSD < best_NRMSD){ 
        best_NRMSD <- current_NRMSD
        best_fit <- impulse_fit
      }
    }
    return(best_fit)
  }
  
  # Define function that will fit the impulse3 model and return the best fit after 100 iterations
  fit_impulse3 <- function(input_time_series, t){
    temp <- data.frame(t=t, M=input_time_series)
    best_NRMSD <- 1E6
    best_fit <- NA
    for(i in 1:100){
      print(str_glue('nlsLM fit - attempt {i} of 100'))
      # Set fit to 0 to start with
      impulse_fit <- 0
      # Continue trying to the model to the data until a fit is found
      while(is.numeric(impulse_fit)){
        impulse_fit <- tryCatch({
          nlsLM(M ~ impulse3(h0, h1, h2, h3, lambda1,lambda2, t, t1, t2, t3), 
                data = temp, start = c(h0=runif(1, min=min(input_time_series), max = max(input_time_series)),
                                       h1=runif(1, min=min(input_time_series), max = max(input_time_series)),
                                       h2=runif(1, min=min(input_time_series), max = max(input_time_series)),
                                       h3=runif(1, min=min(input_time_series), max = max(input_time_series)),
                                       lambda1=runif(n=1, min=0, max=1),
                                       lambda2=runif(n=1, min=0, max=1),
                                       t1=sample(1:max(t),1),
                                       t2=sample(1:max(t),1),
                                       t3=sample(1:max(t),1)),
                control = list(maxiter = 1024,warnOnly=TRUE), trace = F)
        }, error = function(e) { # if nlsLM can't converge, ouput zero
          c(0)
        }, finally = {
          
        })
      }
      # calculate the NRMSD
      current_NRMSD <- sqrt(sum((residuals(impulse_fit))^2)/length(residuals(impulse_fit)))/diff(range(temp$M))
      # Check if current fit is more optimal
      if(current_NRMSD < best_NRMSD){ 
        best_NRMSD <- current_NRMSD
        best_fit <- impulse_fit
      }
    }
    return(best_fit)
  }
  
  # Get:
  # Impulse model fit to P(t) - for plotting purposes only
  # Transcription model fit to M(t)
  # The degradation rate (alpha) given by the transcription model fit
  
  set.seed(1234)
  for (i in 1:nrow(fig_2_data)){
    
    # Get impulse model fit to the pre-mRNA profile
    print(str_glue('Fitting impulse model for gene: {fig_2_data$gene_name[i]}'))
    if(fig_2_data$impulse_model_p[i]=="2"){ fig_2_data[i, "impulse_model_fit"][[1]] <- list(predict(fit_impulse2(input_time_series = unlist(fig_2_data$last_10Kb[[i]]), t=t), data.frame(t=t))) }
    if(fig_2_data$impulse_model_p[i]=="2b"){ fig_2_data[i, "impulse_model_fit"][[1]] <- list(predict(fit_impulse2b(input_time_series = unlist(fig_2_data$last_10Kb[[i]]), t=t), data.frame(t=t))) }
    if(fig_2_data$impulse_model_p[i]=="3"){ fig_2_data[i, "impulse_model_fit"][[1]] <- list(predict(fit_impulse3(input_time_series = unlist(fig_2_data$last_10Kb[[i]]), t=t), data.frame(t=t))) }
    
    # Fit transcription model
    # get expression
    pre_mRNA <- as.numeric(fig_2_data$last_10Kb[[i]])
    mRNA <- as.numeric(fig_2_data$mRNA_expression[[i]])
    # state variable
    state <- c(m=mRNA[1])
    # Time specification.
    times <- seq(0,400,10)
    # Define P(t) as an external variable
    p <- approxfun(as.data.frame(list(times = times, p = pre_mRNA)), rule = 2)
    # Model equations
    get_m <- function(t, state, parameters) {
      with(as.list(c(state, parameters)), {
        # Define rates of change
        P <- p(t)
        dm <- beta*P - alpha*m
        # Return the rate of change
        list(dm)
      })
    }
    # Define function that accepts values for alpha & beta then returns the values for M(t) this gives
    get_fit_m <- function(alpha, beta){
      parameters <- c(alpha=alpha, beta=beta)
      out <- data.frame(ode(y = state, times = times, func = get_m, parms = parameters, method = "euler"))
      out$m
    }
    # Fit the model 
    NRMSD <- 1E10; best_model <- NA
    for(n in 1:10){
      transcription_model_fit <- 0
      while(is.numeric(transcription_model_fit)){
        transcription_model_fit <- tryCatch({
          nlsLM(m ~ get_fit_m(alpha, beta), data = data.frame(m=mRNA), start = list(alpha=runif(1,0,1), beta=runif(1,0,1)), control = list(maxiter = 1024, warnOnly=TRUE), trace = F)
        }, error = function(e) { # if nlsLM can't converge, ouput zero
          c(0)
        }, finally = {
          
        })
      }
      current_NRMSD <- sqrt(sum((residuals(transcription_model_fit))^2)/length(residuals(transcription_model_fit)))/diff(range(mRNA))
      if(current_NRMSD < NRMSD){
        NRMSD <- current_NRMSD
        best_model <- transcription_model_fit
      }
    }
    
    # Save the transcription model fit
    fig_2_data[i, "transcription_model_fit"][[1]] <- list(predict(best_model))
    # Save the inferred degradation rate
    fig_2_data$alpha[i] <- coefficients(best_model)["alpha"]
  }
  
  # Save out plot data
  saveRDS(fig_2_data, "data_track/fig_2_data2.rds")
}

print("Finished creating data required for Fig_2")
