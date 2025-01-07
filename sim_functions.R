### DATA GENERATION ##
# create data from multivariate distribution
MVN <- function(n, d, mu = 0, sigma = 1, rho = 0) {
  CovMat <- function(d, rho = 0, sigma = 1) {
    if (rho > 1 | rho < -1){stop("rho must be between -1 and 1")}
    A <- matrix(data = sigma*rho, nrow = d, ncol = d)
    diag(A) <- sigma
    return(A)
  }
  
  # Set mean vector
  Mu <- rep(mu, d)
  # Set covariance matrix
  Sigma <- CovMat(d, rho, sigma)
  # Generate multivariate normal samples
  X <- rmvnorm(n, Mu, Sigma)
  return(data.frame(X))
}

ContaminateDF <- function(df, size = .10, type = c("MCAR", "MAR", "MNAR"), MAR_rule = NULL) {
  
  if (type == "MCAR") {
    new_df <- as.matrix(df)
    l <- prod(dim(df))
    na_id <- sample(1:l, size = l * size, replace = F)
    new_df[na_id] <- NA
    new_df <- as.data.frame(new_df)
  }
  
  if (type %in% c("MAR", "MNAR")) {
    if (is.null(MAR_rule)) {
      stop("Submit a MAR_rule value with strings inside: e.g. 'c(\"new_df$X1[df$X2 > 5] <- NA\", ...)'.")
    }
    new_df <- df
    # Loop over MAR_rule and evaluate each rule
    for (rule in MAR_rule) {
      eval(parse(text = rule))
    }
  }
  return(new_df)
}

### SIMULATION ###
# function to convert correlation/covariance matrix into df
cor_to_df <- function(mat, type = "Cov") {
  cor_df <- as.data.frame(as.table(mat))
  colnames(cor_df) <- c("V1", "V2", type)
  if (type == "Cov") {
    cor_df <- cor_df[upper.tri(mat, diag = TRUE), ]
  }
  else {
    cor_df <- cor_df[upper.tri(mat, diag = FALSE), ] # don't need 1s
  }
  
  # create column "V1.V2" and transpose for data stacking
  cor_df$Combination <- paste(cor_df$V1, cor_df$V2, sep = ".")
  cor_df_t <- as.data.frame(t(cor_df[, type, drop = FALSE]))
  colnames(cor_df_t) <- cor_df$Combination
  return(cor_df_t)
}

impute_sim <- function(gen_params, na_params) {
  results <- list()
  
  ### Generate Data ###
  n <- gen_params$n
  d <- gen_params$d
  mu <- gen_params$mu # column means
  sigma <- gen_params$sigma # variance of each column
  rho = gen_params$rho # collinearity
  
  # Ground Truth #
  data <- MVN(n, d, mu, sigma, rho)
  mu.t <- colMeans(data)
  sigma.t <- cor_to_df(cov(data), "Cov")
  cor.t <- cor_to_df(cor(data), "Cor")
  results$Truth <- list(
    Mean = mu.t,
    Covariance = sigma.t,
    Correlation = cor.t
  )
  
  ### Contaminate Data with Missing Values ### 
  size <- na_params$size 
  type <- na_params$type 
  MAR_rule <- na_params$MAR_rule
  na_data <- ContaminateDF(data, size, type, MAR_rule)
  mcar_t <- mcar_test(na_data)
  results$MCARp <- mcar_t$p.value
  results$MCARstat <- mcar_t$statistic
  
  ### Methodology ###
  ## Method 1: Listwise Deletion ##
  na_data1 <- drop_na(na_data)
  mu1 <- colMeans(na_data1)
  sigma1 <- cor_to_df(cov(na_data1), "Cov")
  cor1 <- cor_to_df(cor(na_data1), "Cor")
  results$LWD <- list(
    Mean = mu1,
    Covariance = sigma1,
    Correlation = cor1
  )
  
  ## Method 2: Expectation Maximization Algorithm
  mat <- as.matrix(na_data)
  s <- prelim.norm(mat)
  thetahat <- em.norm(s, criterion = 1e-5, showits = F)
  rngseed(sample(1:B, 1)) # needs to be run first for impute to work
  na_data2 <- imp.norm(s, thetahat, mat)
  mu2 <- colMeans(na_data2)
  sigma2 <- cor_to_df(cov(na_data2), "Cov")
  cor2 <- cor_to_df(cor(na_data2), "Cor")
  results$EM <- list(
    Mean = mu2,
    Covariance = sigma2,
    Correlation = cor2
  )
  
  ## Method 3: Sampling Missing Values
  imputed_data <- mice(na_data, method = "norm", m = 5, printFlag = F) 
  na_data3 <- complete(imputed_data, "all") 
  na_data3 <- Reduce("+", na_data3) / length(na_data3)  
  mu3 <- colMeans(na_data3)
  sigma3 <- cor_to_df(cov(na_data3), "Cov")
  cor3 <- cor_to_df(cor(na_data3), "Cor")
  results$CS <- list(
    Mean = mu3,
    Covariance = sigma3,
    Correlation = cor3
  )
  
  return(results)
  
}

### EXTRACTION ###

extract_mcar_values <- function(sim_results) {
  # Extract MCARp and MCARstat from each simulation result
  mcar_values <- lapply(sim_results, function(res) {
    data.frame(
      MCARp = res$MCARp,
      MCARstat = res$MCARstat
    )
  })
  
  # Combine the extracted values into a single dataframe
  mcar_df <- do.call(rbind, mcar_values)
  return(mcar_df)
}

methods <- c("LWD", "EM", "CS")
statistics <- c("Mean", "Covariance", "Correlation")

process_results <- function(bs, methods, statistics, truth_method = "Truth") {
  results_list <- list()
  
  # Extract and combine results for each statistic and method
  for (stat in statistics) {
    truth_stat <- do.call(rbind, lapply(bs, function(x) x[[truth_method]][[stat]]))
    method_results <- lapply(methods, function(method) {
      method_stat <- do.call(rbind, lapply(bs, function(x) x[[method]][[stat]]))
      method_stat - truth_stat  # Calculate deviation from truth
    })
    
    # Combine results and label by method
    combined_results <- do.call(rbind.data.frame, method_results) %>%
      mutate(Method = rep(methods, each = nrow(truth_stat)))
    
    results_list[[stat]] <- combined_results
  }
  
  return(results_list)
}

combine_results <- function(results) {
  # Initialize an empty data frame
  res_combined_mean <- data.frame()
  res_combined_cov <- data.frame()
  res_combined_cor <- data.frame()
  tests_combined <- data.frame()
  
  # Iterate over each row of the results data frame
  for (i in seq_len(nrow(results))) {
    # Extract current parameters
    rho <- results$rho[i]
    type <- results$type[i]
    
    # Expand the res and test lists
    res_mean <- cbind(results$res[[i]]$Mean, 'rho' = rho, 'type' = type)
    res_cov <- cbind(results$res[[i]]$Covariance, 'rho' = rho, 'type' = type)
    res_cor <- cbind(results$res[[i]]$Correlation, 'rho' = rho, 'type' = type)
    tests_subset <- cbind(results$test[[i]], 'rho' = rho, 'type' = type)
    
    # Combine res and test for this iteration
    res_combined_mean <- rbind(res_combined_mean, res_mean)
    res_combined_cov <- rbind(res_combined_cov, res_cov)
    res_combined_cor <- rbind(res_combined_cor, res_cor)
    tests_combined <- rbind(tests_combined, tests_subset)
    
  }
  
  full_list <- list('means' = res_combined_mean, 
                    'cov' = res_combined_cov, 
                    'cor' = res_combined_cor, 
                    'tests' = tests_combined)
  
  return(full_list)
}


## SIMULATING OVER PARAM_GRID ## 

get_mar_rule <- function(type, size) {
  if (type == "MAR") {
    return(c(
      "new_df$X1[df$X2 > qnorm(1-size)] <- NA",
      "new_df$X2[df$X3 > qnorm(1-size)] <- NA",
      "new_df$X3[df$X1 > qnorm(1-size)] <- NA"
    ))
  } else if (type == "MNAR") {
    return(c(
      "new_df$X1[df$X1 > qnorm(1-size)] <- NA",
      "new_df$X2[df$X2 > qnorm(1-size)] <- NA",
      "new_df$X3[df$X3 > qnorm(1-size)] <- NA"
    ))
  } else {
    return(NULL)  # MCAR has no MAR rules
  }
}

# Function to calculate statistics for a dataframe
calculate_statistics <- function(df) {
  # Group by Method
  df %>%
    group_by(Method, rho, type) %>%
    dplyr::summarise(
      # For each numeric column, calculate multiple statistics
      across(
        where(is.numeric),  # Select only numeric columns
        list(
          Mean = ~ mean(.x, na.rm = TRUE),
          Variance = ~ var(.x, na.rm = TRUE)
          # P2.5 = ~ quantile(.x, 0.025, na.rm = TRUE),
          # Median = ~ median(.x, na.rm = TRUE),
          # P97.5 = ~ quantile(.x, 0.975, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"  # Dynamically generate names for the output columns
      ),
      .groups = "drop"  # Remove grouping after summarizing
    )
}

